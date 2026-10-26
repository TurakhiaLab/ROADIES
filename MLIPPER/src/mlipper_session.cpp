#include "mlipper_session.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "io/input_validation.hpp"
#include "io/jplace.hpp"
#include "io/tree_newick.hpp"
#include "gpu/gpu_admission.hpp"
#include "util/model_utils.hpp"
#include "optimize/optimization_types.hpp"
#include "optimize/model_optimization_backend.hpp"
#include "optimize/model_parameter_optimizer.hpp"
#include "placement/placement.cuh"
#include "pmatrix/pmat_gpu.cuh"
#include "spr/local_spr.hpp"
#include "spr/local_spr_internal.hpp"
#include "tree/divide_and_conquer.hpp"
#include "util/checked_size.hpp"
#include "util/mlipper_util.h"

namespace {

constexpr int kSmallTipSPRTopKPerUnit = 8;
constexpr int kSectorNNIRadius = 2;
constexpr int kSectorNNIClusterThreshold = 3;
constexpr int kSectorNNITopKPerUnit = 2048;
constexpr int kSectorNNIRounds = 1;

struct SplitFingerprint {
    size_t side_size = 0;
    std::uint64_t hash_a = 0;
    std::uint64_t hash_b = 0;

    bool operator==(const SplitFingerprint& other) const {
        return side_size == other.side_size &&
            hash_a == other.hash_a && hash_b == other.hash_b;
    }

    bool operator<(const SplitFingerprint& other) const {
        return std::tie(side_size, hash_a, hash_b) <
            std::tie(other.side_size, other.hash_a, other.hash_b);
    }
};

struct SplitFingerprintHash {
    size_t operator()(const SplitFingerprint& value) const {
        size_t seed = std::hash<size_t>{}(value.side_size);
        seed ^= std::hash<std::uint64_t>{}(value.hash_a) +
            0x9e3779b9U + (seed << 6) + (seed >> 2);
        seed ^= std::hash<std::uint64_t>{}(value.hash_b) +
            0x9e3779b9U + (seed << 6) + (seed >> 2);
        return seed;
    }
};

std::uint64_t split_tip_hash(const std::string& name, std::uint64_t seed) {
    std::uint64_t value = seed;
    for (unsigned char byte : name) {
        value ^= static_cast<std::uint64_t>(byte);
        value *= 1099511628211ULL;
    }
    return value;
}

std::unordered_map<int, SplitFingerprint> build_full_tree_split_fingerprints(
    const TreeBuildResult& tree)
{
    std::vector<size_t> descendant_count(tree.nodes.size(), 0);
    std::vector<std::uint64_t> descendant_a(tree.nodes.size(), 0);
    std::vector<std::uint64_t> descendant_b(tree.nodes.size(), 0);
    for (int node_id : tree.postorder) {
        const TreeNode& node = tree.nodes.at(static_cast<size_t>(node_id));
        if (node.is_tip) {
            descendant_count[static_cast<size_t>(node_id)] = 1;
            descendant_a[static_cast<size_t>(node_id)] =
                split_tip_hash(node.name, 1469598103934665603ULL);
            descendant_b[static_cast<size_t>(node_id)] =
                split_tip_hash(node.name, 1099511628211ULL);
        } else {
            descendant_count[static_cast<size_t>(node_id)] =
                mlipper::util::checked_add_size(
                    descendant_count.at(static_cast<size_t>(node.left)),
                    descendant_count.at(static_cast<size_t>(node.right)),
                    "tree split descendant count");
            descendant_a[static_cast<size_t>(node_id)] =
                descendant_a.at(static_cast<size_t>(node.left)) ^
                descendant_a.at(static_cast<size_t>(node.right));
            descendant_b[static_cast<size_t>(node_id)] =
                descendant_b.at(static_cast<size_t>(node.left)) ^
                descendant_b.at(static_cast<size_t>(node.right));
        }
    }
    const size_t total = descendant_count.at(
        static_cast<size_t>(tree.root_id));
    const std::uint64_t total_a = descendant_a.at(
        static_cast<size_t>(tree.root_id));
    const std::uint64_t total_b = descendant_b.at(
        static_cast<size_t>(tree.root_id));
    std::unordered_map<int, SplitFingerprint> result;
    for (const TreeNode& node : tree.nodes) {
        if (node.id == tree.root_id || node.parent < 0 || node.is_tip ||
            tree.nodes.at(static_cast<size_t>(node.parent)).is_tip) {
            continue;
        }
        size_t side_size = descendant_count.at(static_cast<size_t>(node.id));
        std::uint64_t side_a = descendant_a.at(static_cast<size_t>(node.id));
        std::uint64_t side_b = descendant_b.at(static_cast<size_t>(node.id));
        const size_t complement_size = total - side_size;
        const std::uint64_t complement_a = total_a ^ side_a;
        const std::uint64_t complement_b = total_b ^ side_b;
        if (complement_size < side_size ||
            (complement_size == side_size &&
             std::pair{complement_a, complement_b} <
                 std::pair{side_a, side_b})) {
            side_size = complement_size;
            side_a = complement_a;
            side_b = complement_b;
        }
        if (side_size >= 2 && total - side_size >= 2) {
            result.emplace(
                node.id,
                SplitFingerprint{side_size, side_a, side_b});
        }
    }
    return result;
}

using SplitSet = std::unordered_set<SplitFingerprint, SplitFingerprintHash>;
using BranchEdgeSet = std::unordered_set<int>;

std::vector<BranchEdgeSet> partition_branch_edge_cores(
    const TreeBuildResult& tree,
    const BranchEdgeSet& branches,
    int core_edges)
{
    if (core_edges <= 0) {
        throw std::invalid_argument(
            "partition_branch_edge_cores: core_edges must be positive");
    }
    std::unordered_map<int, std::vector<int>> incident;
    for (int child_id : branches) {
        if (child_id < 0 ||
            child_id >= static_cast<int>(tree.nodes.size())) {
            throw std::runtime_error(
                "partition_branch_edge_cores: branch child is out of range");
        }
        const int parent_id =
            tree.nodes[static_cast<size_t>(child_id)].parent;
        if (parent_id < 0) {
            throw std::runtime_error(
                "partition_branch_edge_cores: branch has no parent");
        }
        incident[child_id].push_back(child_id);
        incident[parent_id].push_back(child_id);
    }

    std::unordered_map<int, std::vector<int>> adjacency;
    for (int child_id : branches) {
        adjacency.emplace(child_id, std::vector<int>{});
    }
    for (auto& [node_id, edges] : incident) {
        (void)node_id;
        std::sort(edges.begin(), edges.end());
        for (int edge : edges) {
            for (int neighbor : edges) {
                if (edge != neighbor) {
                    adjacency[edge].push_back(neighbor);
                }
            }
        }
    }
    for (auto& [edge, neighbors] : adjacency) {
        (void)edge;
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(
            std::unique(neighbors.begin(), neighbors.end()),
            neighbors.end());
    }

    BranchEdgeSet remaining = branches;
    std::vector<BranchEdgeSet> cores;
    while (!remaining.empty()) {
        const int seed = *std::min_element(
            remaining.begin(), remaining.end());
        std::deque<int> queue{seed};
        remaining.erase(seed);
        BranchEdgeSet core;
        while (!queue.empty() &&
               core.size() < static_cast<size_t>(core_edges)) {
            const int edge = queue.front();
            queue.pop_front();
            core.insert(edge);
            for (int neighbor : adjacency.at(edge)) {
                if (core.size() >= static_cast<size_t>(core_edges) ||
                    queue.size() >=
                        static_cast<size_t>(core_edges) - core.size()) {
                    break;
                }
                if (remaining.erase(neighbor) != 0) {
                    queue.push_back(neighbor);
                }
            }
        }
        // Return queued but unprocessed edges to the pool. They will seed a
        // later connected sector instead of being silently dropped.
        while (!queue.empty()) {
            remaining.insert(queue.front());
            queue.pop_front();
        }
        cores.push_back(std::move(core));
    }
    return cores;
}

struct WholeTreeSplitIndex {
    std::unordered_map<SplitFingerprint, std::pair<int, int>, SplitFingerprintHash>
        endpoints;
    std::vector<std::vector<int>> graph;
};

WholeTreeSplitIndex build_whole_tree_split_index(const TreeBuildResult& tree)
{
    WholeTreeSplitIndex index;
    index.graph.resize(tree.nodes.size());
    for (const TreeNode& node : tree.nodes) {
        if (node.parent < 0) continue;
        index.graph[static_cast<size_t>(node.id)].push_back(node.parent);
        index.graph[static_cast<size_t>(node.parent)].push_back(node.id);
    }
    for (auto& neighbors : index.graph) {
        std::sort(neighbors.begin(), neighbors.end());
    }
    const auto fingerprints = build_full_tree_split_fingerprints(tree);
    for (const auto& [child, fingerprint] : fingerprints) {
        const int parent = tree.nodes.at(static_cast<size_t>(child)).parent;
        auto [it, inserted] = index.endpoints.emplace(
            fingerprint, std::pair{child, parent});
        if (!inserted && std::pair{child, parent} < it->second) {
            it->second = {child, parent};
        }
    }
    return index;
}

std::unordered_map<SplitFingerprint, SplitSet, SplitFingerprintHash>
whole_tree_edge_adjacency(
    const WholeTreeSplitIndex& index,
    const SplitSet& keys)
{
    std::unordered_map<int, std::vector<SplitFingerprint>> incident;
    for (const SplitFingerprint& key : keys) {
        const auto found = index.endpoints.find(key);
        if (found == index.endpoints.end()) continue;
        incident[found->second.first].push_back(key);
        incident[found->second.second].push_back(key);
    }
    std::unordered_map<SplitFingerprint, SplitSet, SplitFingerprintHash> adjacency;
    for (const SplitFingerprint& key : keys) adjacency.emplace(key, SplitSet{});
    for (const auto& [node, edges] : incident) {
        (void)node;
        for (const SplitFingerprint& edge : edges) {
            for (const SplitFingerprint& other : edges) {
                if (!(edge == other)) adjacency[edge].insert(other);
            }
        }
    }
    return adjacency;
}

std::vector<SplitSet> partition_owned_edge_cores(
    const WholeTreeSplitIndex& index,
    SplitSet remaining,
    int core_edges)
{
    const auto adjacency = whole_tree_edge_adjacency(index, remaining);
    std::vector<SplitSet> cores;
    while (!remaining.empty()) {
        const SplitFingerprint seed = *std::min_element(
            remaining.begin(), remaining.end());
        std::unordered_map<
            SplitFingerprint, SplitFingerprint, SplitFingerprintHash> parent;
        SplitSet has_parent;
        std::vector<SplitFingerprint> order{seed};
        has_parent.insert(seed);
        for (size_t cursor = 0; cursor < order.size(); ++cursor) {
            std::vector<SplitFingerprint> neighbors(
                adjacency.at(order[cursor]).begin(),
                adjacency.at(order[cursor]).end());
            std::sort(neighbors.begin(), neighbors.end());
            for (const SplitFingerprint& neighbor : neighbors) {
                if (remaining.count(neighbor) &&
                    !has_parent.count(neighbor)) {
                    has_parent.insert(neighbor);
                    parent.emplace(neighbor, order[cursor]);
                    order.push_back(neighbor);
                }
            }
        }
        SplitSet core;
        if (order.size() <= static_cast<size_t>(core_edges)) {
            core.insert(order.begin(), order.end());
        } else {
            std::unordered_map<SplitFingerprint, int, SplitFingerprintHash>
                subtree_size;
            for (const SplitFingerprint& key : order) subtree_size[key] = 1;
            for (size_t idx = order.size(); idx-- > 1;) {
                subtree_size[parent.at(order[idx])] += subtree_size[order[idx]];
            }
            bool found_root = false;
            SplitFingerprint root{};
            int root_size = -1;
            for (const SplitFingerprint& key : order) {
                const int size = subtree_size[key];
                if (size <= core_edges &&
                    (!found_root || size > root_size ||
                     (size == root_size && root < key))) {
                    root = key;
                    root_size = size;
                    found_root = true;
                }
            }
            std::deque<SplitFingerprint> queue{root};
            while (!queue.empty()) {
                const SplitFingerprint key = queue.front();
                queue.pop_front();
                core.insert(key);
                for (const SplitFingerprint& neighbor : adjacency.at(key)) {
                    const auto parent_it = parent.find(neighbor);
                    if (parent_it != parent.end() &&
                        parent_it->second == key) {
                        queue.push_back(neighbor);
                    }
                }
            }
        }
        for (const SplitFingerprint& key : core) remaining.erase(key);
        cores.push_back(std::move(core));
    }
    return cores;
}

mlipper::model_optimization::SiteBatchSelection
choose_global_optimization_site_batch(
    size_t nodes,
    size_t tips,
    size_t sites,
    size_t rate_categories,
    size_t states,
    bool per_rate_scaling,
    size_t free_gpu_bytes,
    mlipper::model_optimization::OptimizationWorkspaceKind workspace_kind)
{
    auto selection =
        mlipper::model_optimization::chooseOptimizationSiteBatchSize(
            nodes, tips, sites, rate_categories, states, per_rate_scaling,
            workspace_kind,
            free_gpu_bytes);
    return selection;
}

void copy_host_to_device(void* dst, const void* src, size_t bytes);

int count_live_tip_nodes(const TreeBuildResult& tree)
{
    int tip_count = 0;
    for (const TreeNode& node : tree.nodes) {
        if (node.is_tip) {
            ++tip_count;
        }
    }
    return std::max(1, tip_count);
}

struct FullTreeRegraftCoordinate {
    int target_child_global_id = -1;
    double proximal_length = 0.0;
};

FullTreeRegraftCoordinate map_local_regraft_to_full_tree(
    const TreeBuildResult& local_tree,
    const std::vector<int>& local_to_global,
    const TreeBuildResult& full_tree,
    int local_target_child_id,
    double local_proximal_length)
{
    if (local_target_child_id < 0 ||
        local_target_child_id >= static_cast<int>(local_tree.nodes.size())) {
        throw std::runtime_error(
            "map_local_regraft_to_full_tree: local target child id is out of range.");
    }
    if (local_to_global.size() != local_tree.nodes.size()) {
        throw std::runtime_error(
            "map_local_regraft_to_full_tree: local-to-global mapping size does not match the local subtree.");
    }

    const TreeNode& local_target_child =
        local_tree.nodes[static_cast<size_t>(local_target_child_id)];
    if (local_target_child.parent < 0) {
        throw std::runtime_error(
            "map_local_regraft_to_full_tree: local regraft target is the root edge.");
    }

    const int local_parent_id = local_target_child.parent;
    const int global_parent_id =
        local_to_global[static_cast<size_t>(local_parent_id)];
    const int global_child_id =
        local_to_global[static_cast<size_t>(local_target_child_id)];
    if (global_parent_id < 0 ||
        global_parent_id >= static_cast<int>(full_tree.nodes.size()) ||
        global_child_id < 0 ||
        global_child_id >= static_cast<int>(full_tree.nodes.size())) {
        throw std::runtime_error(
            "map_local_regraft_to_full_tree: local edge endpoint does not map into the full tree.");
    }

    const std::vector<int> full_path =
        mlipper::divide_and_conquer::build_tree_path_nodes(
            full_tree,
            global_parent_id,
            global_child_id);
    if (full_path.size() < 2 ||
        full_path.front() != global_parent_id ||
        full_path.back() != global_child_id) {
        throw std::runtime_error(
            "map_local_regraft_to_full_tree: full-tree path does not match the local edge orientation.");
    }

    double remaining_proximal =
        std::max(0.0, local_proximal_length);
    double total_path_length = 0.0;
    for (size_t path_idx = 1; path_idx < full_path.size(); ++path_idx) {
        const int previous_id = full_path[path_idx - 1];
        const int current_id = full_path[path_idx];
        const TreeNode& previous =
            full_tree.nodes[static_cast<size_t>(previous_id)];
        const TreeNode& current =
            full_tree.nodes[static_cast<size_t>(current_id)];
        if (current.parent == previous_id) {
            total_path_length +=
                static_cast<double>(current.branch_length_to_parent);
        } else if (previous.parent == current_id) {
            total_path_length +=
                static_cast<double>(previous.branch_length_to_parent);
        } else {
            throw std::runtime_error(
                "map_local_regraft_to_full_tree: local edge endpoints are not connected in the full tree.");
        }
    }
    remaining_proximal =
        std::min(remaining_proximal, total_path_length);

    for (size_t path_idx = 1; path_idx < full_path.size(); ++path_idx) {
        const int previous_id = full_path[path_idx - 1];
        const int current_id = full_path[path_idx];
        const TreeNode& previous =
            full_tree.nodes[static_cast<size_t>(previous_id)];
        const TreeNode& current =
            full_tree.nodes[static_cast<size_t>(current_id)];
        const bool walking_down = current.parent == previous_id;
        const int target_child_id = walking_down ? current_id : previous_id;
        const double edge_length = walking_down
            ? static_cast<double>(current.branch_length_to_parent)
            : static_cast<double>(previous.branch_length_to_parent);
        if (remaining_proximal <= edge_length ||
            path_idx + 1 == full_path.size()) {
            return FullTreeRegraftCoordinate{
                target_child_id,
                walking_down
                    ? std::clamp(remaining_proximal, 0.0, edge_length)
                    : std::clamp(
                          edge_length - remaining_proximal,
                          0.0,
                          edge_length),
            };
        }
        remaining_proximal -= edge_length;
    }

    throw std::runtime_error(
        "map_local_regraft_to_full_tree: failed to resolve a full-tree regraft coordinate.");
}

void replay_local_spr_move_on_full_tree(
    TreeBuildResult& full_tree,
    const TreeBuildResult& local_tree,
    const std::vector<int>& local_to_global,
    int local_prune_root_id,
    int local_target_child_id,
    double pendant_length,
    double proximal_length)
{
    if (local_prune_root_id < 0 ||
        local_prune_root_id >= static_cast<int>(local_to_global.size())) {
        throw std::runtime_error(
            "replay_local_spr_move_on_full_tree: local prune root id is out of range.");
    }
    const int full_prune_root_id =
        local_to_global[static_cast<size_t>(local_prune_root_id)];
    if (full_prune_root_id < 0 ||
        full_prune_root_id >= static_cast<int>(full_tree.nodes.size())) {
        throw std::runtime_error(
            "replay_local_spr_move_on_full_tree: local prune root does not map into the full tree.");
    }

    TreeBuildResult local_pruned = local_tree;
    PruneInfo local_prune_info;
    if (!prune_subtree_for_spr(
            local_pruned,
            local_prune_root_id,
            local_prune_info)) {
        throw std::runtime_error(
            "replay_local_spr_move_on_full_tree: failed to prune the local subtree.");
    }

    // A local move can fail while mapping a suppressed edge back into the full
    // topology. Build it off to the side so the authoritative tree is unchanged
    // unless every prune/map/regraft step succeeds.
    TreeBuildResult candidate_full_tree = full_tree;
    PruneInfo full_prune_info;
    if (!prune_subtree_for_spr(
            candidate_full_tree, full_prune_root_id, full_prune_info)) {
        throw std::runtime_error(
            "replay_local_spr_move_on_full_tree: failed to prune the mapped full-tree subtree.");
    }
    // Candidate edges and proximal lengths are defined on the topology after
    // prune suppression. Map that topology, not the stale pre-prune edge.
    const FullTreeRegraftCoordinate full_regraft =
        map_local_regraft_to_full_tree(
            local_pruned,
            local_to_global,
            candidate_full_tree,
            local_target_child_id,
            proximal_length);
    regraft_subtree_for_spr(
        candidate_full_tree,
        full_prune_info,
        full_regraft.target_child_global_id,
        pendant_length,
        full_regraft.proximal_length,
        candidate_full_tree.nodes[
            static_cast<size_t>(full_prune_info.pruned_id)].is_tip
            ? OPT_BRANCH_LEN_MIN
            : TOPOLOGY_INTERNAL_BRANCH_LEN_MIN,
        candidate_full_tree.nodes[static_cast<size_t>(
            full_regraft.target_child_global_id)].is_tip
            ? OPT_BRANCH_LEN_MIN
            : TOPOLOGY_INTERNAL_BRANCH_LEN_MIN);
    full_tree = std::move(candidate_full_tree);
}

void copy_local_branch_lengths_to_full_tree(
    const TreeBuildResult& local_tree,
    const std::vector<int>& local_to_global,
    TreeBuildResult& full_tree)
{
    for (size_t local_id = 0;
         local_id < local_tree.nodes.size() &&
         local_id < local_to_global.size();
         ++local_id) {
        const int global_id = local_to_global[local_id];
        const int local_parent = local_tree.nodes[local_id].parent;
        if (global_id < 0 || local_parent < 0 ||
            local_parent >= static_cast<int>(local_to_global.size())) {
            continue;
        }
        const int global_parent =
            local_to_global[static_cast<size_t>(local_parent)];
        if (global_parent < 0 ||
            global_id >= static_cast<int>(full_tree.nodes.size()) ||
            full_tree.nodes[static_cast<size_t>(global_id)].parent !=
                global_parent) {
            continue;
        }
        full_tree.nodes[static_cast<size_t>(global_id)]
            .branch_length_to_parent =
            local_tree.nodes[local_id].branch_length_to_parent;
    }
}

void copy_host_to_device(void* dst, const void* src, size_t bytes) {
    CUDA_CHECK(cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice));
}

template <typename T>
void cuda_free_and_clear(T*& pointer)
{
    if (pointer == nullptr) return;
    CUDA_CHECK(cudaFree(pointer));
    pointer = nullptr;
}

void require_matching_device(
    int candidate_device,
    int& expected_device,
    const char* context,
    const char* detail)
{
    if (candidate_device < 0) {
        return;
    }
    if (expected_device >= 0 && expected_device != candidate_device) {
        throw std::runtime_error(std::string(context) + ": " + detail);
    }
    expected_device = candidate_device;
}

bool query_device_pointer_device(const void* ptr, int* device_out) {
    if (ptr == nullptr) {
        return false;
    }

    cudaPointerAttributes attrs{};
    const cudaError_t err = cudaPointerGetAttributes(&attrs, ptr);
    if (err != cudaSuccess) {
        cudaGetLastError();
        return false;
    }

#if CUDART_VERSION >= 10000
    const auto memory_type = attrs.type;
#else
    const auto memory_type = attrs.memoryType;
#endif

    if (memory_type != cudaMemoryTypeDevice &&
        memory_type != cudaMemoryTypeManaged) {
        return false;
    }

    if (device_out != nullptr) {
        *device_out = attrs.device;
    }
    return true;
}

void require_alignment_matches_device_tree(
    const parse::Alignment& alignment,
    const DeviceTree& device_tree,
    const char* context)
{
    if (alignment.names.size() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(
            std::string(context) +
            ": alignment tip count exceeds the DeviceTree integer limit.");
    }
    if (device_tree.sites != alignment.sites) {
        throw std::runtime_error(
            std::string(context) +
            ": alignment site count does not match DeviceTree.");
    }
    if (device_tree.tips != static_cast<int>(alignment.names.size())) {
        throw std::runtime_error(
            std::string(context) +
            ": alignment tip count does not match DeviceTree.");
    }
}

void require_model_matches_device_tree(
    const parse::ModelConfig& model_config,
    const DeviceTree& device_tree,
    const char* context)
{
    if (device_tree.states > 0 && device_tree.states != model_config.states) {
        throw std::runtime_error(
            std::string(context) +
            ": DeviceTree state count does not match model.");
    }
}

parse::ModelConfig resolved_model_config(
    parse::ModelConfig model_config,
    bool use_empirical_freqs,
    const parse::Alignment* tree_alignment,
    const std::vector<unsigned>* pattern_weights = nullptr)
{
    if (use_empirical_freqs) {
        if (tree_alignment != nullptr &&
            !tree_alignment->names.empty() &&
            tree_alignment->sites > 0) {
            model_config.freqs =
                mlipper::model::estimate_empirical_pi(
                    *tree_alignment,
                    model_config.states,
                    pattern_weights ? *pattern_weights
                                    : std::vector<unsigned>{});
        } else {
            model_config.freqs.clear();
        }
    }

    mlipper::input::validate_model_inputs(model_config);
    return model_config;
}

} // namespace

namespace mlipper {

MlipperSession::MlipperSession() = default;

MlipperSession::~MlipperSession()
{
    try {
        releaseOwnedDeviceTreeMemory();
    } catch (...) {
        // Destruction cannot propagate CUDA cleanup failures. Clear the
        // independently owned workspaces so each destructor remains idempotent.
        releasePlacementOps();
        releaseSubtreeWorkspace();
        releaseTopologyRefinementWorkspaces();
    }
    releaseGpuReservation();
}

void MlipperSession::releaseGpuReservation() noexcept {
    gpu_reservation_.reset();
    active_gpu_device_ = -1;
}

cudaStream_t MlipperSession::currentStream() const noexcept {
    return cudaStream_t{};
}

LocalSPRSessionWorkspaceSet& MlipperSession::ensureTopologyRefinementWorkspaces()
{
    if (!topology_refinement_session_workspaces_) {
        topology_refinement_session_workspaces_ =
            std::make_unique<LocalSPRSessionWorkspaceSet>();
    }
    return *topology_refinement_session_workspaces_;
}

void MlipperSession::clearPendingPlacement() {
    pending_queries_ = PlacementQueryBatch{};
    pending_query_names_.clear();
}

void MlipperSession::releasePlacementOps() noexcept {
    try {
        int device_to_use = active_gpu_device_;
        if (device_to_use < 0 && owned_device_tree_.device_id >= 0) {
            device_to_use = owned_device_tree_.device_id;
        }
        if (placement_ops_.d_ops != nullptr && device_to_use >= 0) {
            gpu::set_device_or_throw(device_to_use);
            free_placement_op_buffer(placement_ops_, currentStream());
        }
    } catch (...) {
    }

    placement_ops_ = PlacementOpBuffer{};
    clearPendingPlacement();
}

void MlipperSession::releaseSubtreeWorkspace() noexcept {
    try {
        int device_to_use = active_gpu_device_;
        if (device_to_use < 0 &&
            subtree_workspace_.dev.device_id >= 0) {
            device_to_use = subtree_workspace_.dev.device_id;
        }
        if (device_to_use >= 0) {
            gpu::set_device_or_throw(device_to_use);
        }
        if (subtree_placement_ops_.d_ops != nullptr) {
            free_placement_op_buffer(
                subtree_placement_ops_,
                currentStream());
        }
        release_subtree_workspace(subtree_workspace_);
    } catch (...) {
    }

    subtree_placement_ops_ = PlacementOpBuffer{};
    subtree_workspace_ = OwnedSubtreeWorkspace{};
}

void MlipperSession::releaseDivideAndConquerSubtreeCache()
{
    CUDA_CHECK(cudaStreamSynchronize(currentStream()));
    releaseSubtreeWorkspace();
    // Sector NNI retains ranking/evaluation scratch sized for the largest
    // sector seen so far. A D&C batch boundary is a safe point to release it
    // because the following batch does not reuse its CLVs or candidates.
    releaseTopologyRefinementWorkspaces();
}

void MlipperSession::installStreamingDirectionalBoundaryMessages(
    const std::vector<divide_and_conquer::DirectionalBoundaryPort>& ports)
{
    DeviceTree& destination = subtree_workspace_.dev;
    if (destination.sites != tree_alignment_.sites) {
        throw std::runtime_error(
            "streaming boundary destination does not span the full alignment.");
    }
    launch_init_tip_clv(destination, currentStream());
    OwnedSubtreeWorkspace component_workspace;
    PlacementOpBuffer component_ops;
    component_ops.tuning = subtree_placement_ops_.tuning;
    size_t free_bytes = 0;
    size_t total_bytes = 0;
    CUDA_CHECK(cudaMemGetInfo(&free_bytes, &total_bytes));
    const auto boundary_selection =
        model_optimization::chooseOptimizationSiteBatchSize(
            cpu_tree_.nodes.size(),
            cpu_host_packing_.tip_node_ids.size(),
            tree_alignment_.sites,
            static_cast<size_t>(model_config_.ncat),
            static_cast<size_t>(model_config_.states),
            model_config_.per_rate_scaling,
            model_optimization::OptimizationWorkspaceKind::GlobalParameters,
            free_bytes);
        const size_t boundary_site_batch = boundary_selection.batch_sites;
        if (boundary_site_batch == 0) {
            throw std::runtime_error(
                "streaming boundary batching selected a zero-sized batch");
        }
    {
        const size_t clv_site_span = mlipper::util::checked_product(
            "streaming boundary CLV site span",
            static_cast<size_t>(model_config_.ncat),
            static_cast<size_t>(model_config_.states));
        const size_t scaler_site_span = model_config_.per_rate_scaling
            ? static_cast<size_t>(model_config_.ncat)
            : 1;
        try {
            for (size_t site_begin = 0;
                 site_begin < tree_alignment_.sites;
                 site_begin += boundary_site_batch) {
                const size_t site_count = std::min(
                    boundary_site_batch, tree_alignment_.sites - site_begin);
                const size_t site_end = mlipper::util::checked_add_size(
                    site_begin, site_count, "streaming boundary site range");
                if (site_end > static_cast<size_t>(
                        std::numeric_limits<std::ptrdiff_t>::max())) {
                    throw std::overflow_error(
                        "streaming boundary site range exceeds iterator limits");
                }
                std::vector<std::string> rows;
                rows.reserve(tree_alignment_.sequences.size());
                for (const std::string& sequence : tree_alignment_.sequences) {
                    rows.push_back(sequence.substr(site_begin, site_count));
                }
                HostPacking host = pack_host_arrays_from_tree_and_msa(
                    cpu_tree_, tree_alignment_.names, rows, site_count,
                    model_config_.states);
                host.pattern_weights.assign(
                    cpu_host_packing_.pattern_weights.begin() +
                        static_cast<std::ptrdiff_t>(site_begin),
                    cpu_host_packing_.pattern_weights.begin() +
                        static_cast<std::ptrdiff_t>(site_end));
                fill_pmats_in_host_packing(
                    cpu_tree_, host, cpu_eig_, cpu_rate_multipliers_,
                    model_config_.states, model_config_.ncat);
                load_subtree_workspace(
                    component_workspace,
                    SubtreeWorkspaceLoadConfig{
                        cpu_tree_, host, cpu_eig_, cpu_rate_weights_,
                        cpu_rate_multipliers_, cpu_pi_, site_count,
                        model_config_.states, model_config_.ncat,
                        model_config_.per_rate_scaling, nullptr, false, 0, 0,
                    },
                    currentStream(),
                    "installStreamingDirectionalBoundaryMessages(shared)");
                UpdateTreeClvsUpwardOnly(
                    component_workspace.dev, cpu_tree_, host,
                    component_ops, currentStream());
                for (const auto& port : ports) {
                    if (port.outside_global_node_id < 0 ||
                        port.outside_global_node_id >= component_workspace.dev.N ||
                        port.local_node_id < 0 ||
                        port.local_node_id >= destination.N) {
                        throw std::runtime_error(
                            "shared streaming boundary port is invalid");
                    }
                    CUDA_CHECK(cudaMemcpyAsync(
                        destination.d_clv_up +
                            static_cast<size_t>(port.local_node_id) *
                                destination.per_node_elems() +
                            site_begin * clv_site_span,
                        component_workspace.dev.d_clv_up +
                            static_cast<size_t>(port.outside_global_node_id) *
                                component_workspace.dev.per_node_elems(),
                        mlipper::util::checked_allocation_bytes<fp_t>(
                            mlipper::util::checked_product(
                                "streaming boundary CLV copy",
                                site_count, clv_site_span),
                            "streaming boundary CLV copy"),
                        cudaMemcpyDeviceToDevice, currentStream()));
                    CUDA_CHECK(cudaMemcpyAsync(
                        destination.d_site_scaler_up +
                            static_cast<size_t>(port.local_node_id) *
                                destination.scaler_elems() +
                            site_begin * scaler_site_span,
                        component_workspace.dev.d_site_scaler_up +
                            static_cast<size_t>(port.outside_global_node_id) *
                                component_workspace.dev.scaler_elems(),
                        mlipper::util::checked_allocation_bytes<unsigned>(
                            mlipper::util::checked_product(
                                "streaming boundary scaler copy",
                                site_count, scaler_site_span),
                            "streaming boundary scaler copy"),
                        cudaMemcpyDeviceToDevice, currentStream()));
                }
                CUDA_CHECK(cudaStreamSynchronize(currentStream()));
                if (component_ops.d_ops != nullptr) {
                    free_placement_op_buffer(component_ops, currentStream());
                }
                release_subtree_workspace(component_workspace);
            }
        } catch (...) {
            if (component_ops.d_ops != nullptr) {
                free_placement_op_buffer(component_ops, currentStream());
            }
            release_subtree_workspace(component_workspace);
            throw;
        }
    }
}

void MlipperSession::releaseTopologyRefinementWorkspaces() noexcept
{
    if (!topology_refinement_session_workspaces_) {
        return;
    }
    try {
        release_local_spr_session_workspace_set(
            *topology_refinement_session_workspaces_,
            currentStream());
    } catch (...) {
    }
    topology_refinement_session_workspaces_.reset();
}

void MlipperSession::loadSubtreeWorkspace(
    TreeBuildResult& tree,
    HostPacking& host_packing,
    const PlacementQueryBatch& queries,
    const char* context,
    bool commit_to_tree,
    int insertion_capacity)
{
    if (queries.size() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(
            std::string(context) +
            ": query count exceeds the subtree workspace integer limit");
    }
    const bool needs_rebuild =
        subtree_workspace_requires_rebuild(
            subtree_workspace_,
            tree,
            host_packing,
            tree_alignment_.sites,
            model_config_.states,
            model_config_.ncat,
            model_config_.per_rate_scaling);
    if (needs_rebuild) {
        if (active_gpu_device_ < 0) {
            throw std::runtime_error(
                std::string(context) +
                ": no active CUDA device is attached to the session.");
        }

        gpu::set_device_or_throw(active_gpu_device_);

        PlacementTuningConfig tuning = placement_ops_.tuning;
        if (subtree_placement_ops_.d_ops != nullptr) {
            tuning = subtree_placement_ops_.tuning;
        }

        releaseSubtreeWorkspace();
        subtree_placement_ops_.tuning = tuning;
    }

    load_subtree_workspace(
        subtree_workspace_,
        SubtreeWorkspaceLoadConfig{
            tree,
            host_packing,
            cpu_eig_,
            cpu_rate_weights_,
            cpu_rate_multipliers_,
            cpu_pi_,
            tree_alignment_.sites,
            model_config_.states,
            model_config_.ncat,
            model_config_.per_rate_scaling,
            &queries,
            commit_to_tree,
            static_cast<int>(queries.size()),
            std::max(0, insertion_capacity),
        },
        currentStream(),
        context);
    UpdateTreeClvs(
        subtree_workspace_.dev,
        tree,
        host_packing,
        subtree_placement_ops_,
        currentStream());
}

void MlipperSession::ensurePlacementReady(const char* context) const {
    if (!has_initialized_cpu_state_) {
        throw std::runtime_error(
            std::string(context) +
            ": initializeCPU must complete before placement.");
    }
    if (gpu_mode_ != GpuMode::DivideAndConquer &&
        owned_device_tree_.d_tipchars == nullptr) {
        throw std::runtime_error(
            std::string(context) +
            ": initializeGPU must bind a DeviceTree before placement.");
    }
    if (active_gpu_device_ < 0) {
        throw std::runtime_error(
            std::string(context) +
            ": no active CUDA device is attached to the session.");
    }
    if (gpu_mode_ != GpuMode::DivideAndConquer &&
        (placement_ops_.d_ops == nullptr || placement_ops_.num_ops <= 0)) {
        throw std::runtime_error(
            std::string(context) +
            ": placement ops are not initialized. initializeGPU must finish successfully first.");
    }
}

SequenceRecord MlipperSession::buildLoadedQueryRecord(int query_index) const
{
    if (!has_query_alignment_) {
        throw std::runtime_error(
            "buildLoadedQueryRecord: no query alignment has been loaded into the session.");
    }
    if (query_index < 0) {
        throw std::runtime_error(
            "buildLoadedQueryRecord: query_index must be >= 0.");
    }
    if (query_index >= static_cast<int>(query_alignment_.names.size())) {
        throw std::runtime_error(
            "buildLoadedQueryRecord: query_index is out of range for the loaded query alignment.");
    }

    SequenceRecord query;
    query.name = query_alignment_.names[static_cast<size_t>(query_index)];
    query.sequence = query_alignment_.sequences[static_cast<size_t>(query_index)];
    return query;
}

PlacementResult MlipperSession::makePlacementResultForTarget(
    int target_id,
    double loglikelihood,
    double distal_length,
    double pendant_length,
    const std::string& tip_name) const
{
    if (target_id < 0 || target_id >= static_cast<int>(cpu_tree_.nodes.size())) {
        throw std::runtime_error(
            "makePlacementResultForTarget: target node id is out of range.");
    }

    const TreeNode& target = cpu_tree_.nodes[static_cast<size_t>(target_id)];
    if (target.parent < 0) {
        throw std::runtime_error(
            "makePlacementResultForTarget: target edge must not be the root edge.");
    }

    PlacementResult out;
    out.parent_node_label =
        cpu_tree_.nodes[static_cast<size_t>(target.parent)].stable_node_label;
    out.child_node_label = target.stable_node_label;
    out.parent_side_branch_length =
        static_cast<double>(target.branch_length_to_parent) - distal_length;
    out.child_side_branch_length = distal_length;
    out.pendant_branch_length = pendant_length;
    out.log_likelihood = loglikelihood;
    out.query_name = tip_name;
    return out;
}

int MlipperSession::prepareGpuDevice(
    const MlipperGpuConfig& gpu_config,
    int estimated_process_memory_mb,
    gpu::DeviceReservation* transferred_reservation)
{
    if (transferred_reservation != nullptr &&
        transferred_reservation->device >= 0) {
        releaseGpuReservation();
        gpu_reservation_ = std::move(*transferred_reservation);
        gpu::set_device_or_throw(gpu_reservation_.device);
        active_gpu_device_ = gpu::current_device_or_throw();
        return active_gpu_device_;
    }

    if (gpu_config.acquire_mode !=
        MlipperGpuAcquireMode::UseCurrentDevice) {
        releaseGpuReservation();
    }
    gpu::DeviceReservation selected_reservation =
        gpu::select_device_or_wait_or_throw(
        gpu_config,
        estimated_process_memory_mb);

    active_gpu_device_ = gpu::current_device_or_throw();
    if (selected_reservation.device >= 0) {
        gpu_reservation_ = std::move(selected_reservation);
    } else if (gpu_reservation_.device >= 0 &&
               gpu_reservation_.device != active_gpu_device_) {
        releaseGpuReservation();
        active_gpu_device_ = gpu::current_device_or_throw();
    } else {
        gpu::ensure_reservation_capacity_or_wait(
            gpu_reservation_,
            estimated_process_memory_mb);
    }
    return active_gpu_device_;
}

void MlipperSession::clearInitializedCpuState()
{
    // Every field below is derived from the loaded tree, alignment, or model.
    // Clearing them together prevents a later GPU initialization from mixing
    // representations built from different input generations.
    cpu_tree_ = TreeBuildResult{};
    cpu_host_packing_ = HostPacking{};
    cpu_eig_ = EigResult{};
    cpu_rate_weights_.clear();
    cpu_rate_multipliers_.clear();
    cpu_pi_.clear();
    next_stable_node_label_ = 0;
    clearPendingPlacement();
    has_initialized_cpu_state_ = false;
}

void MlipperSession::validateLoadedInputsAgainstDeviceTree(
    const DeviceTree& device_tree,
    const char* context) const
{
    if (has_tree_alignment_) {
        require_alignment_matches_device_tree(
            tree_alignment_, device_tree, context);
    }
    if (has_model_config_) {
        require_model_matches_device_tree(
            model_config_, device_tree, context);
    }
}

std::string MlipperSession::loadBackboneTree(const std::string& newick_text)
{
    if (newick_text.empty()) {
        throw std::runtime_error(
            "loadBackboneTree: tree input is empty.");
    }

    // Parsing and rooted/unrooted normalization are performed once by the
    // libpll-backed tree builder during initializeCPU().
    backbone_tree_newick_ = newick_text;
    has_backbone_tree_ = true;
    clearInitializedCpuState();
    return backbone_tree_newick_;
}

const parse::Alignment& MlipperSession::loadAlignment(
    const parse::Alignment& tree_alignment)
{
    if (tree_alignment.names.empty() ||
        tree_alignment.sequences.empty() ||
        tree_alignment.sites == 0) {
        throw std::runtime_error(
            "loadAlignment: tree alignment must contain at least one sequence and one site.");
    }
    // Public session APIs may be used without the CLI validation layer. Validate
    // the complete candidate before replacing the currently loaded generation.
    input::validate_alignment_names(tree_alignment, "loadAlignment");
    input::validate_alignment_symbols(tree_alignment, 4, "loadAlignment");

    tree_alignment_ = tree_alignment;
    has_tree_alignment_ = true;
    query_alignment_ = parse::Alignment{};
    has_query_alignment_ = false;
    has_custom_pattern_weights_ = false;
    loaded_pattern_weights_.clear();
    clearInitializedCpuState();
    return tree_alignment_;
}

const parse::Alignment& MlipperSession::loadAlignment(
    const parse::Alignment& tree_alignment,
    const parse::Alignment& query_alignment)
{
    if (tree_alignment.names.empty() ||
        tree_alignment.sequences.empty() ||
        tree_alignment.sites == 0) {
        throw std::runtime_error(
            "loadAlignment: tree alignment must contain at least one sequence and one site.");
    }
    if (query_alignment.names.empty() ||
        query_alignment.sequences.empty() ||
        query_alignment.sites == 0) {
        throw std::runtime_error(
            "loadAlignment: query alignment must contain at least one sequence and one site.");
    }
    if (query_alignment.sites != tree_alignment.sites) {
        throw std::runtime_error(
            "loadAlignment: query alignment site count does not match tree alignment.");
    }

    input::validate_alignment_names(tree_alignment, "loadAlignment reference");
    input::validate_alignment_symbols(
        tree_alignment, 4, "loadAlignment reference");
    input::validate_alignment_names(query_alignment, "loadAlignment query");
    input::validate_alignment_symbols(
        query_alignment, 4, "loadAlignment query");
    input::validate_query_reference_name_overlap(
        tree_alignment, query_alignment, "loadAlignment");

    tree_alignment_ = tree_alignment;
    query_alignment_ = query_alignment;
    has_tree_alignment_ = true;
    has_query_alignment_ = true;
    has_custom_pattern_weights_ = false;
    loaded_pattern_weights_.clear();
    clearInitializedCpuState();
    return tree_alignment_;
}

void MlipperSession::setPatternWeights(
    const std::vector<unsigned>& pattern_weights)
{
    if (!has_tree_alignment_) {
        throw std::runtime_error(
            "setPatternWeights: alignment must be loaded before setting pattern weights.");
    }
    if (!pattern_weights.empty() &&
        pattern_weights.size() != tree_alignment_.sites) {
        throw std::runtime_error(
            "setPatternWeights: pattern weight count does not match the loaded alignment site count.");
    }

    loaded_pattern_weights_ = pattern_weights;
    has_custom_pattern_weights_ = true;
    clearInitializedCpuState();
}

const parse::ModelConfig& MlipperSession::loadModel(
    const parse::ModelConfig& model_config,
    bool model_uses_empirical_freqs)
{
    input::validate_model_inputs(model_config);
    model_config_ = model_config;
    model_uses_empirical_freqs_ = model_uses_empirical_freqs;
    has_model_config_ = true;
    clearInitializedCpuState();
    return model_config_;
}

void MlipperSession::ensureGpuResources()
{
    int expected_device = -1;
    if (owned_device_tree_.d_tipchars != nullptr) {
        require_matching_device(
            owned_device_tree_.device_id,
            expected_device,
            "ensureGpuResources",
            "owned DeviceTree device is inconsistent.");
    }
    require_matching_device(
        active_gpu_device_,
        expected_device,
        "ensureGpuResources",
        "session GPU device does not match current session device.");

    if (owned_device_tree_.d_tipchars != nullptr) {
        validateLoadedInputsAgainstDeviceTree(
            owned_device_tree_,
            "ensureGpuResources");
        int tipchars_device = -1;
        if (!query_device_pointer_device(
                owned_device_tree_.d_tipchars,
                &tipchars_device)) {
            throw std::runtime_error(
                "ensureGpuResources: owned DeviceTree tipchars pointer is not a valid CUDA device buffer.");
        }
        require_matching_device(
            tipchars_device,
            expected_device,
            "ensureGpuResources",
            "owned DeviceTree tipchars device is inconsistent.");
        owned_device_tree_.device_id = expected_device;
    }

    if (owned_device_tree_.d_tipchars != nullptr && expected_device >= 0) {
        owned_device_tree_.device_id = expected_device;
    }
}

void MlipperSession::initializeCPU()
{
    if (!has_backbone_tree_) {
        throw std::runtime_error(
            "initializeCPU: backbone tree must be loaded before CPU initialization.");
    }
    if (!has_tree_alignment_) {
        throw std::runtime_error(
            "initializeCPU: alignment must be loaded before CPU initialization.");
    }
    if (!has_model_config_) {
        throw std::runtime_error(
            "initializeCPU: model must be loaded before CPU initialization.");
    }

    input::validate_alignment_names(
        tree_alignment_,
        "initializeCPU");
    input::validate_alignment_symbols(
        tree_alignment_,
        model_config_.states,
        "initializeCPU");

    // Resolve empirical parameters before clearing derived state: empirical
    // frequencies depend on the currently loaded alignment and pattern weights.
    parse::ModelConfig effective_model_config = resolved_model_config(
        model_config_,
        model_uses_empirical_freqs_,
        &tree_alignment_,
        has_custom_pattern_weights_ ? &loaded_pattern_weights_ : nullptr);
    model_config_ = effective_model_config;

    clearInitializedCpuState();

    const int states = effective_model_config.states;
    const int rate_cats = effective_model_config.ncat;
    const size_t sites = tree_alignment_.sites;

    // Build the model and topology independently, then join them in HostPacking.
    // HostPacking is the canonical upload layout consumed by DeviceTree setup.
    cpu_pi_ = model::ensure_normalized_pi(effective_model_config.freqs, states);
    cpu_rate_weights_ = model::build_discrete_gamma_weights(rate_cats);
    cpu_rate_multipliers_ = model::build_gamma_rate_categories(
        effective_model_config.alpha, rate_cats);

    const std::vector<double> q_rowmajor =
        build_gtr_q_matrix(states, effective_model_config, cpu_pi_);

    cpu_tree_ = build_tree_from_newick_with_pll(
        tree_alignment_.names,
        backbone_tree_newick_);

    if (cpu_tree_.preorder.size() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(
            "initializeCPU: tree has too many nodes for stable integer labels");
    }
    int label = 0;
    for (int node_id : cpu_tree_.preorder) {
        cpu_tree_.nodes[node_id].stable_node_label = label++;
    }
    next_stable_node_label_ = static_cast<int>(cpu_tree_.preorder.size());

    cpu_host_packing_ = pack_host_arrays_from_tree_and_msa(
        cpu_tree_,
        tree_alignment_.names,
        tree_alignment_.sequences,
        sites,
        states);
    if (has_custom_pattern_weights_) {
        cpu_host_packing_.pattern_weights = loaded_pattern_weights_;
    } else {
        cpu_host_packing_.pattern_weights.assign(sites, 1u);
    }

    cpu_eig_ = gtr_eigendecomp_cpu(
        q_rowmajor.data(),
        cpu_pi_.data(),
        states);
    fill_pmats_in_host_packing(
        cpu_tree_,
        cpu_host_packing_,
        cpu_eig_,
        cpu_rate_multipliers_,
        states,
        rate_cats);

    has_initialized_cpu_state_ = true;
}

bool MlipperSession::usesModelSiteBatching() const noexcept
{
    return gpu_mode_ == GpuMode::SiteBatchedModel &&
        owned_device_tree_.d_tipchars != nullptr &&
        model_site_batch_size_ > 0;
}

model_optimization::ModelOptimizationContext
MlipperSession::modelOptimizationContext()
{
    return model_optimization::ModelOptimizationContext{
        has_initialized_cpu_state_,
        owned_device_tree_,
        cpu_tree_,
        cpu_host_packing_,
        placement_ops_,
        model_config_,
        model_uses_empirical_freqs_,
        cpu_eig_,
        cpu_rate_multipliers_,
        cpu_pi_,
        model_site_batch_workspace_,
        model_site_batch_size_,
        backbone_tree_newick_,
        usesModelSiteBatching(),
        currentStream(),
    };
}

void MlipperSession::enterGlobalParameterOptimizationMode(
    bool require_single_site_batch,
    GlobalOptimizationBackend backend)
{
    if (!has_initialized_cpu_state_ ||
        owned_device_tree_.d_tipchars == nullptr) {
        throw std::runtime_error(
            "enterGlobalParameterOptimizationMode: CPU and GPU state must be initialized");
    }
    if (owned_device_tree_.sites == 0) {
        throw std::runtime_error(
            "enterGlobalParameterOptimizationMode: alignment has no sites");
    }

    if (backend ==
        GlobalOptimizationBackend::ResidentSequential) {
        model_site_batch_workspace_.release();
        model_site_batch_size_ = 0;
        gpu_mode_ = GpuMode::ResidentTree;
        UpdateTreeClvs(
            owned_device_tree_, cpu_tree_, cpu_host_packing_,
            placement_ops_, currentStream());
        return;
    }

    releasePlacementWorkspace();
    model_site_batch_workspace_.release();
    size_t free_gpu_bytes = 0;
    size_t total_gpu_bytes = 0;
    CUDA_CHECK(cudaMemGetInfo(&free_gpu_bytes, &total_gpu_bytes));
    auto selection = choose_global_optimization_site_batch(
        static_cast<size_t>(owned_device_tree_.N),
        static_cast<size_t>(owned_device_tree_.tips),
        owned_device_tree_.sites,
        static_cast<size_t>(owned_device_tree_.rate_cats),
        static_cast<size_t>(owned_device_tree_.states),
        owned_device_tree_.per_rate_scaling,
        free_gpu_bytes,
        require_single_site_batch
            ? model_optimization::OptimizationWorkspaceKind::AllBranchGradient
            : model_optimization::OptimizationWorkspaceKind::GlobalParameters);
    if (require_single_site_batch &&
        selection.batch_sites < owned_device_tree_.sites) {
        const size_t full_workspace_bytes = mlipper::util::checked_product(
            "global optimization full workspace",
            selection.estimated_bytes_per_site,
            owned_device_tree_.sites);
        // Leave headroom for the resident tree/model allocations and CUDA
        // runtime bookkeeping which are not part of the workspace estimate.
        const size_t usable_free_bytes = free_gpu_bytes - free_gpu_bytes / 5;
        if (full_workspace_bytes > usable_free_bytes) {
            throw std::runtime_error(
                "global branch optimization requires one site batch: "
                "estimated full workspace=" +
                std::to_string(full_workspace_bytes) +
                " bytes, usable free GPU memory=" +
                std::to_string(usable_free_bytes) + " bytes");
        }
        selection.batch_sites = owned_device_tree_.sites;
        selection.batch_count = 1;
        selection.estimated_workspace_bytes = full_workspace_bytes;
    }
    model_site_batch_size_ = selection.batch_sites;
    gpu_mode_ = GpuMode::SiteBatchedModel;
}

FinalModelOptimizationResult MlipperSession::runFinalModelOptimization(
    const FinalModelOptimizationOptions& options)
{
    if (owned_device_tree_.d_tipchars == nullptr) {
        throw std::runtime_error("runFinalModelOptimization: GPU state is not initialized");
    }
    if (model_config_.pinv > 0.0) {
        throw std::runtime_error(
            "runFinalModelOptimization does not yet support models with pinv > 0");
    }
    enterGlobalParameterOptimizationMode(
        options.optimize_branch_lengths, options.backend);
    model_optimization::ModelOptimizationBackend backend(
        modelOptimizationContext());
    (void)backend.evaluate();
    FinalModelOptimizationResult result = global_model_optimizer_.optimize(
        cpu_pi_, model_config_.rates, model_config_.alpha,
        model_config_.states, owned_device_tree_.rate_cats,
        options, backend.callbacks());
    backend.refreshBranchLengthTransitionMatrices(false, false);
    return result;
}

size_t MlipperSession::loadedQueryCount() const noexcept
{
    return query_alignment_.names.size();
}

PlacementBatchResult
MlipperSession::findBestLoadedPlacements()
{
    if (!has_query_alignment_ || query_alignment_.names.empty()) {
        throw std::runtime_error(
            "findBestLoadedPlacements: a non-empty query alignment is required.");
    }
    std::vector<SequenceRecord> queries;
    queries.reserve(query_alignment_.names.size());
    for (size_t index = 0; index < query_alignment_.names.size(); ++index) {
        queries.push_back(buildLoadedQueryRecord(static_cast<int>(index)));
    }
    return findBestPlacements(queries);
}

PlacementBatchResult MlipperSession::findBestPlacements(
    const std::vector<SequenceRecord>& queries)
{
    ensurePlacementReady("findBestPlacements");
    if (queries.empty()) {
        throw std::runtime_error(
            "findBestPlacements: queries must not be empty.");
    }
    if (queries.size() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(
            "findBestPlacements: query count exceeds the DeviceTree integer limit");
    }
    for (const SequenceRecord& query : queries) {
        if (query.sequence.empty()) {
            throw std::runtime_error(
                "findBestPlacements: query sequence must not be empty.");
        }
    }
    DeviceTree& device_tree = owned_device_tree_;
    ensure_device_tree_current_device(
        device_tree,
        "MlipperSession::findBestPlacements");
    ensure_device_tree_query_capacity(
        device_tree,
        static_cast<int>(queries.size()),
        "MlipperSession::findBestPlacements");

    pending_query_names_.clear();
    pending_query_names_.reserve(queries.size());
    std::vector<std::string> query_rows;
    query_rows.reserve(queries.size());
    for (const SequenceRecord& query : queries) {
        pending_query_names_.push_back(query.name);
        query_rows.push_back(query.sequence);
    }
    pending_queries_ =
        divide_and_conquer::build_query_batch_from_placement_queries(
            build_placement_query(
                pending_query_names_,
                query_rows),
            tree_alignment_.sites,
            model_config_.states);
    const size_t query_chars_bytes =
        mlipper::util::checked_allocation_bytes<uint8_t>(
            pending_queries_.query_chars.size(),
            "findBestPlacements query characters");
    CUDA_CHECK(cudaMemcpyAsync(
        device_tree.d_query_chars,
        pending_queries_.query_chars.data(),
        query_chars_bytes,
        cudaMemcpyHostToDevice,
        currentStream()));

    PlacementCommitContext commit_ctx;
    commit_ctx.placement_ops = &placement_ops_;

    std::vector<::RawPlacementResult> placement_results;
    EvaluatePlacementQueries(
        device_tree,
        cpu_eig_,
        cpu_rate_multipliers_,
        commit_ctx,
        &placement_results,
        1,
        false,
        currentStream());
    if (placement_results.size() != queries.size()) {
        throw std::runtime_error(
            "findBestPlacements: placement result count does not match query count.");
    }

    PlacementBatchResult result;
    result.placements.reserve(queries.size());
    result.query_names.reserve(queries.size());
    for (size_t idx = 0; idx < queries.size(); ++idx) {
        const ::RawPlacementResult& placement = placement_results[idx];
        result.placements.push_back(makePlacementResultForTarget(
            placement.target_id,
            placement.loglikelihood,
            placement.distal_length,
            placement.pendant_length,
            queries[idx].name));
        result.query_names.push_back(queries[idx].name);
    }
    result.ranked_placements = std::move(placement_results);
    return result;
}

std::vector<PlacementResult> MlipperSession::placeAndCommitQueries(
    const std::vector<SequenceRecord>& queries)
{
    ensurePlacementReady("placeAndCommitQueries");
    if (queries.empty()) {
        throw std::runtime_error(
            "placeAndCommitQueries: queries must not be empty.");
    }
    if (queries.size() >
        static_cast<size_t>(std::numeric_limits<int>::max() / 2)) {
        throw std::overflow_error(
            "placeAndCommitQueries: query count exceeds insertion ID capacity");
    }
    for (const SequenceRecord& query : queries) {
        if (query.sequence.empty()) {
            throw std::runtime_error(
                "placeAndCommitQueries: query sequence must not be empty.");
        }
    }

    DeviceTree& device_tree = owned_device_tree_;
    ensure_device_tree_current_device(
        device_tree,
        "MlipperSession::placeAndCommitQueries");
    const int query_count = static_cast<int>(queries.size());
    const long long required_nodes =
        static_cast<long long>(device_tree.N) + 2LL * query_count;
    const long long required_tips =
        static_cast<long long>(device_tree.tips) + query_count;
    if (required_nodes > device_tree.capacity_N ||
        required_tips > device_tree.capacity_tips) {
        throw std::runtime_error(
            "placeAndCommitQueries: DeviceTree does not have reserved insertion capacity. Reinitialize GPU with commit_to_tree enabled.");
    }

    ensure_device_tree_query_capacity(
        device_tree,
        query_count,
        "MlipperSession::placeAndCommitQueries");
    pending_query_names_.clear();
    pending_query_names_.reserve(queries.size());
    std::vector<std::string> query_rows;
    query_rows.reserve(queries.size());
    for (const SequenceRecord& query : queries) {
        pending_query_names_.push_back(query.name);
        query_rows.push_back(query.sequence);
    }
    pending_queries_ =
        divide_and_conquer::build_query_batch_from_placement_queries(
            build_placement_query(pending_query_names_, query_rows),
            tree_alignment_.sites,
            model_config_.states);
    const size_t query_chars_bytes =
        mlipper::util::checked_allocation_bytes<uint8_t>(
            pending_queries_.query_chars.size(),
            "placeAndCommitQueries query characters");
    CUDA_CHECK(cudaMemcpyAsync(
        device_tree.d_query_chars,
        pending_queries_.query_chars.data(),
        query_chars_bytes,
        cudaMemcpyHostToDevice,
        currentStream()));

    PlacementCommitContext commit_ctx;
    commit_ctx.tree = &cpu_tree_;
    commit_ctx.host = &cpu_host_packing_;
    commit_ctx.queries = &pending_queries_;
    commit_ctx.placement_ops = &placement_ops_;
    commit_ctx.query_names = &pending_query_names_;
    std::vector<std::string> inserted_query_names(queries.size());
    commit_ctx.inserted_query_names = &inserted_query_names;

    // Commit mode coordinates CPU topology, HostPacking, and DeviceTree updates
    // inside EvaluatePlacementQueries. Each query appends one attachment node
    // followed by one tip, which determines the IDs below. The multi-stage
    // commit is not yet transactional if a later allocation fails.
    const int initial_node_count = device_tree.N;
    std::vector<::RawPlacementResult> kernel_results;
    EvaluatePlacementQueries(
        device_tree,
        cpu_eig_,
        cpu_rate_multipliers_,
        commit_ctx,
        &kernel_results,
        1,
        true,
        currentStream());
    if (kernel_results.size() != queries.size()) {
        throw std::runtime_error(
            "placeAndCommitQueries: placement result count does not match query count.");
    }

    std::vector<PlacementResult> committed_results;
    committed_results.reserve(queries.size());
    for (size_t index = 0; index < queries.size(); ++index) {
        const int internal_id =
            initial_node_count + 2 * static_cast<int>(index);
        const int tip_id = internal_id + 1;
        // Node IDs can change when a topology is rebuilt. Stable labels are the
        // workflow-facing identity used to carry placement anchors across it.
        const int internal_label = next_stable_node_label_++;
        const int tip_label = next_stable_node_label_++;
        cpu_tree_.nodes[static_cast<size_t>(internal_id)].stable_node_label =
            internal_label;
        cpu_tree_.nodes[static_cast<size_t>(tip_id)].stable_node_label = tip_label;

        const ::RawPlacementResult& kernel = kernel_results[index];
        PlacementResult committed;
        committed.child_node_label =
            cpu_tree_.nodes[static_cast<size_t>(kernel.target_id)].stable_node_label;
        const TreeNode& attachment_node =
            cpu_tree_.nodes[static_cast<size_t>(internal_id)];
        const int original_parent_id = attachment_node.parent;
        committed.parent_node_label = original_parent_id >= 0
            ? cpu_tree_.nodes[static_cast<size_t>(original_parent_id)].stable_node_label
            : -1;
        committed.attachment_node_label = internal_label;
        committed.query_tip_label = tip_label;
        committed.parent_side_branch_length =
            static_cast<double>(attachment_node.branch_length_to_parent);
        committed.child_side_branch_length = static_cast<double>(
            cpu_tree_.nodes[static_cast<size_t>(kernel.target_id)]
                .branch_length_to_parent);
        committed.pendant_branch_length = static_cast<double>(
            cpu_tree_.nodes[static_cast<size_t>(tip_id)]
                .branch_length_to_parent);
        committed.log_likelihood = kernel.loglikelihood;
        committed.query_name = inserted_query_names[index].empty()
            ? queries[index].name
            : inserted_query_names[index];

        tree_alignment_.names.push_back(committed.query_name);
        tree_alignment_.sequences.push_back(queries[index].sequence);
        committed_results.push_back(std::move(committed));
    }
    CUDA_CHECK(cudaStreamSynchronize(currentStream()));
    clearPendingPlacement();
    return committed_results;
}

void MlipperSession::runSmallTipBatches(
    const MlipperPlacementParams& params,
    const MlipperLocalSPRParams& local_spr_params,
    const MlipperGpuConfig& gpu_config)
{
    if (!has_initialized_cpu_state_) {
        throw std::runtime_error(
            "runSmallTipBatches: initializeCPU must be called first.");
    }
    if (!has_query_alignment_ || query_alignment_.names.empty()) {
        throw std::runtime_error(
            "runSmallTipBatches: a non-empty query alignment is required.");
    }
    if (!params.commit_to_tree) {
        throw std::runtime_error(
            "runSmallTipBatches: commit_to_tree must be enabled.");
    }

    MlipperPlacementParams effective_params = params;
    if (loadedQueryCount() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(
            "runSmallTipBatches: query count exceeds the workflow integer limit");
    }
    const int total_queries = static_cast<int>(loadedQueryCount());
    const int batch_size = effective_params.insertion_capacity > 0
        ? std::min(effective_params.insertion_capacity, total_queries)
        : total_queries;
    effective_params.insertion_capacity = batch_size;
    initializeGPU(effective_params, gpu_config);

    const MlipperGpuConfig recycle_gpu_config{};
    for (int batch_start = 0;
         batch_start < total_queries;
         batch_start += batch_size) {
        const int batch_end =
            std::min(batch_start + batch_size, total_queries);
        std::vector<SequenceRecord> batch_queries;
        batch_queries.reserve(batch_end - batch_start);
        for (int index = batch_start; index < batch_end; ++index) {
            batch_queries.push_back(buildLoadedQueryRecord(index));
        }

        const std::vector<PlacementResult> committed =
            placeAndCommitQueries(batch_queries);
        if (effective_params.local_spr) {
            runLocalSPR(
                local_spr_params,
                committed);
        }
        if (batch_end < total_queries) {
            // The current resident allocation reserves one batch of insertion
            // slots. Rebuild from the newly committed CPU tree before reusing
            // those slots for the next batch.
            initializeGPU(effective_params, recycle_gpu_config);
        }
    }
}

HostPacking MlipperSession::buildCurrentFullTreeHostPacking() const
{
    HostPacking host_packing = pack_host_arrays_from_tree_and_msa(
        cpu_tree_,
        tree_alignment_.names,
        tree_alignment_.sequences,
        tree_alignment_.sites,
        model_config_.states);
    host_packing.pattern_weights = cpu_host_packing_.pattern_weights;
    if (host_packing.pattern_weights.empty() &&
        !has_custom_pattern_weights_) {
        host_packing.pattern_weights.assign(tree_alignment_.sites, 1u);
    }
    fill_pmats_in_host_packing(
        cpu_tree_,
        host_packing,
        cpu_eig_,
        cpu_rate_multipliers_,
        model_config_.states,
        model_config_.ncat);
    return host_packing;
}

void MlipperSession::rebuildCpuTreeFromCurrentTopology()
{
    // Topology is authoritative in cpu_tree_. Re-serialize it for future
    // rebuilds, regenerate host arrays, then refresh the resident GPU view when
    // this workflow owns one.
    backbone_tree_newick_ = treeio::write_tree_to_newick_string(cpu_tree_);
    cpu_host_packing_ = buildCurrentFullTreeHostPacking();

    if (owned_device_tree_.d_tipchars == nullptr) {
        return;
    }

    DeviceTree& device_tree = owned_device_tree_;
    ensure_device_tree_current_device(
        device_tree,
        "MlipperSession::rebuildCpuTreeFromCurrentTopology");
    reload_device_tree_live_data(
        device_tree,
        cpu_tree_,
        cpu_host_packing_,
        nullptr,
        currentStream());
    UpdateTreeClvs(
        device_tree,
        cpu_tree_,
        cpu_host_packing_,
        placement_ops_,
        currentStream());
}

void MlipperSession::runTopologyRefinement(
    const TopologyRefinementParams& params,
    const std::vector<divide_and_conquer::TreeEdgeEndpoints>& anchors,
    TopologyMoveType move_type,
    const std::vector<NNIOwnedSplit>& owned_nni_splits)
{
    ensurePlacementReady("runTopologyRefinement");
    if (params.radius < 0) {
        throw std::runtime_error(
            "runTopologyRefinement: radius must be >= 0.");
    }
    if (params.cluster_threshold < 0) {
        throw std::runtime_error(
            "runTopologyRefinement: cluster threshold must be >= 0.");
    }
    if (params.topk_per_unit <= 0) {
        throw std::runtime_error(
            "runTopologyRefinement: top-k per unit must be > 0.");
    }
    if (params.rounds <= 0) {
        throw std::runtime_error(
            "runTopologyRefinement: rounds must be >= 1.");
    }
    if (anchors.empty() && move_type != TopologyMoveType::NNI) {
        throw std::runtime_error(
            "runTopologyRefinement: anchors must not be empty.");
    }

    // Small-tip SPR edits the resident tree directly. D&C NNI instead operates
    // on a bounded local tree whose boundary CLVs preserve the likelihood
    // contribution of the omitted full-tree components.
    const bool use_directional_subtree =
        move_type == TopologyMoveType::NNI;
    SplitSet owned_nni;
    std::vector<divide_and_conquer::TreeEdgeEndpoints> nni_central_edges;
    if (move_type == TopologyMoveType::NNI) {
        for (const NNIOwnedSplit& split : owned_nni_splits) {
            owned_nni.insert(SplitFingerprint{
                split.side_size, split.hash_a, split.hash_b});
        }
        const WholeTreeSplitIndex nni_index =
            build_whole_tree_split_index(cpu_tree_);
        nni_central_edges.reserve(owned_nni.size());
        for (const SplitFingerprint& fingerprint : owned_nni) {
            const auto endpoint_it = nni_index.endpoints.find(fingerprint);
            if (endpoint_it == nni_index.endpoints.end()) {
                continue;
            }
            nni_central_edges.push_back({
                endpoint_it->second.first,
                endpoint_it->second.second,
            });
        }
        if (nni_central_edges.size() != owned_nni.size()) {
            throw std::runtime_error(
                "runTopologyRefinement: not every owned NNI split exists in the current tree.");
        }
    }
    TreeBuildResult* active_tree = &cpu_tree_;
    HostPacking* active_host_packing = &cpu_host_packing_;
    DeviceTree* active_device_tree =
        gpu_mode_ == GpuMode::DivideAndConquer
            ? nullptr
            : &owned_device_tree_;
    PlacementOpBuffer* active_placement_ops = &placement_ops_;
    std::vector<std::string>* active_alignment_names = &tree_alignment_.names;
    std::vector<std::string>* active_alignment_rows = &tree_alignment_.sequences;
    divide_and_conquer::PreparedDirectionalSubtree directional_subtree;
    std::vector<int> local_to_global;

    if (use_directional_subtree) {
        directional_subtree =
            divide_and_conquer::prepare_nni_sector_with_directional_boundaries(
                cpu_tree_,
                nni_central_edges,
                tree_alignment_.names,
                tree_alignment_.sequences,
                tree_alignment_.sites,
                model_config_.states,
                cpu_eig_,
                cpu_rate_multipliers_,
                model_config_.ncat);
        local_to_global = directional_subtree.local_to_global;
        directional_subtree.host_packing.pattern_weights =
            cpu_host_packing_.pattern_weights;

        const PlacementQueryBatch empty_queries;
        loadSubtreeWorkspace(
            directional_subtree.subtree,
            directional_subtree.host_packing,
            empty_queries,
            "MlipperSession::runTopologyRefinement",
            true,
            0);
        installStreamingDirectionalBoundaryMessages(
            directional_subtree.boundary_ports);
        UpdateTreeClvsPreservingTipClvs(
            subtree_workspace_.dev,
            directional_subtree.subtree,
            directional_subtree.host_packing,
            subtree_placement_ops_,
            currentStream());
        // From this point the shared refinement engine is agnostic to workflow:
        // these aliases select either resident full-tree state or the prepared
        // directional subtree state.
        active_tree = &directional_subtree.subtree;
        active_host_packing = &directional_subtree.host_packing;
        active_device_tree = &subtree_workspace_.dev;
        active_placement_ops = &subtree_placement_ops_;
        active_alignment_names =
            &directional_subtree.alignment.names;
        active_alignment_rows =
            &directional_subtree.alignment.sequences;
    }
    ensure_device_tree_current_device(
        *active_device_tree,
        "MlipperSession::runTopologyRefinement");

    TopologyRefinementState refinement_state{};
    refinement_state.device = *active_device_tree;
    refinement_state.tree = *active_tree;
    refinement_state.host_packing = *active_host_packing;
    refinement_state.eig = cpu_eig_;
    refinement_state.queries = PlacementQueryBatch{};

    std::string current_tree_newick;
    if (use_directional_subtree) {
        TreeBuildResult output_tree = refinement_state.tree;
        for (const divide_and_conquer::DirectionalBoundaryPort& port :
             directional_subtree.boundary_ports) {
            output_tree.nodes[static_cast<size_t>(port.local_node_id)].is_tip = true;
        }
        current_tree_newick =
            treeio::write_tree_to_newick_string(output_tree);
    } else {
        current_tree_newick =
            treeio::write_tree_to_newick_string(refinement_state.tree);
    }
    LocalSPRSessionWorkspaceSet& refinement_workspaces =
        ensureTopologyRefinementWorkspaces();
    TopologyRefinementRunContext refinement_ctx(
        refinement_state,
        *active_placement_ops,
        currentStream(),
        *active_alignment_names,
        *active_alignment_rows,
        current_tree_newick,
        refinement_state.host_packing.pattern_weights,
        cpu_rate_weights_,
        cpu_rate_multipliers_,
        cpu_pi_,
        tree_alignment_.sites,
        model_config_.states,
        model_config_.ncat,
        model_config_.per_rate_scaling,
        anchors,
        refinement_workspaces);
    refinement_ctx.radius = params.radius;
    refinement_ctx.cluster_threshold = params.cluster_threshold;
    refinement_ctx.topk_per_unit = params.topk_per_unit;
    refinement_ctx.rounds = params.rounds;
    refinement_ctx.move_type = move_type;
    refinement_ctx.symmetric_branch_sweeps =
        gpu_mode_ == GpuMode::DivideAndConquer ? 1 : 0;
    if (move_type == TopologyMoveType::NNI) {
        std::cout << "MLIPPER local topology mode: exact-internal-edge-NNI\n";
    }
    if (use_directional_subtree) {
        refinement_ctx.directional_boundary_node_ids.reserve(
            directional_subtree.boundary_ports.size());
        for (const divide_and_conquer::DirectionalBoundaryPort& port :
             directional_subtree.boundary_ports) {
            refinement_ctx.directional_boundary_node_ids.push_back(
                port.local_node_id);
            if (port.outside_global_node_id >= 0 &&
                port.outside_global_node_id <
                    static_cast<int>(cpu_tree_.nodes.size()) &&
                cpu_tree_.nodes[static_cast<size_t>(
                    port.outside_global_node_id)].is_tip) {
                refinement_ctx.directional_boundary_tip_node_ids.insert(
                    port.local_node_id);
            }
        }
    }
    if (move_type == TopologyMoveType::NNI) {
        if (!owned_nni.empty()) {
            const auto fingerprints =
                build_full_tree_split_fingerprints(cpu_tree_);
            std::unordered_map<SplitFingerprint, int, SplitFingerprintHash>
                local_child_by_owned_split;
            for (const TreeNode& local_node : refinement_state.tree.nodes) {
                if (local_node.parent < 0 || local_node.is_tip ||
                    local_node.id >= static_cast<int>(local_to_global.size()) ||
                    local_node.parent >= static_cast<int>(local_to_global.size())) {
                    continue;
                }
                const int global_id =
                    local_to_global[static_cast<size_t>(local_node.id)];
                const int global_parent = local_to_global[
                    static_cast<size_t>(local_node.parent)];
                if (global_id < 0 || global_parent < 0) continue;
                const SplitFingerprint* fingerprint = nullptr;
                const auto direct = fingerprints.find(global_id);
                if (direct != fingerprints.end() &&
                    cpu_tree_.nodes.at(static_cast<size_t>(global_id)).parent ==
                        global_parent) {
                    fingerprint = &direct->second;
                } else {
                    const auto reverse = fingerprints.find(global_parent);
                    if (reverse != fingerprints.end() &&
                        cpu_tree_.nodes.at(static_cast<size_t>(global_parent)).parent ==
                            global_id) {
                        fingerprint = &reverse->second;
                    }
                }
                if (fingerprint != nullptr && owned_nni.count(*fingerprint) != 0) {
                    auto [it, inserted] = local_child_by_owned_split.emplace(
                        *fingerprint, local_node.id);
                    if (!inserted && local_node.id < it->second) {
                        it->second = local_node.id;
                    }
                }
            }
            if (local_child_by_owned_split.size() != owned_nni.size()) {
                throw std::runtime_error(
                    "compact D&C NNI sector did not contain every owned edge: owned=" +
                    std::to_string(owned_nni.size()) + " local=" +
                    std::to_string(local_child_by_owned_split.size()));
            }
            for (const auto& [fingerprint, local_child_id] :
                 local_child_by_owned_split) {
                (void)fingerprint;
                refinement_ctx.allowed_nni_central_edge_child_ids.insert(
                    local_child_id);
            }
        }
    }
    if (use_directional_subtree) {
        refinement_ctx.on_accept_move =
            [&](const TreeBuildResult& local_tree_before_move,
                int prune_root_id,
                int regraft_child_id,
                double pendant_length,
                double proximal_length) {
                // Accepted local NNI moves must be replayed immediately because
                // cpu_tree_ remains the authoritative topology across sectors.
                replay_local_spr_move_on_full_tree(
                    cpu_tree_,
                    local_tree_before_move,
                    local_to_global,
                    prune_root_id,
                    regraft_child_id,
                    pendant_length,
                    proximal_length);
            };
    }
    run_topology_refinement(refinement_ctx);
    last_topology_refinement_accepted_move_count_ =
        refinement_ctx.accepted_move_count;
    if (use_directional_subtree &&
        refinement_ctx.accepted_move_count > 0 &&
        refinement_ctx.symmetric_branch_sweeps > 0) {
        copy_local_branch_lengths_to_full_tree(
            refinement_ctx.state.tree,
            local_to_global,
            cpu_tree_);
    }
    if (use_directional_subtree) {
        releaseTopologyRefinementWorkspaces();
    }

    if (use_directional_subtree) {
        if (refinement_ctx.accepted_move_count > 0) {
            rebuildCpuTreeFromCurrentTopology();
        }
        return;
    }
    cpu_tree_ = std::move(refinement_state.tree);
    cpu_host_packing_ = std::move(refinement_state.host_packing);
    backbone_tree_newick_ = std::move(current_tree_newick);
}

void MlipperSession::runLocalSPR(
    const MlipperLocalSPRParams& params,
    const std::vector<PlacementResult>& recent_committed_placements)
{
    const TopologyRefinementParams refinement_params{
        params.local_spr_radius,
        params.local_spr_cluster_threshold,
        kSmallTipSPRTopKPerUnit,
        params.local_spr_rounds,
    };
    const std::vector<divide_and_conquer::TreeEdgeEndpoints> anchors =
        build_local_spr_insertion_anchors(
            cpu_tree_,
            recent_committed_placements);
    runTopologyRefinement(
        refinement_params,
        anchors,
        TopologyMoveType::SPR,
        {});
}

void MlipperSession::runSectorNNI(
    const std::vector<NNIOwnedSplit>& owned_splits)
{
    if (gpu_mode_ != GpuMode::DivideAndConquer) {
        throw std::runtime_error(
            "runSectorNNI requires divide-and-conquer mode");
    }
    if (owned_splits.empty()) {
        throw std::invalid_argument("runSectorNNI requires owned edges");
    }
    const TopologyRefinementParams refinement_params{
        kSectorNNIRadius,
        kSectorNNIClusterThreshold,
        kSectorNNITopKPerUnit,
        kSectorNNIRounds,
    };
    runTopologyRefinement(
        refinement_params,
        {},
        TopologyMoveType::NNI,
        owned_splits);
}

DivideAndConquerNNIResult MlipperSession::runDivideAndConquerNNI(
    const DivideAndConquerNNIOptions& options)
{
    if (options.core_edges <= 0 || options.max_sweeps <= 0) {
        throw std::invalid_argument("invalid whole-tree NNI options");
    }
    DivideAndConquerNNIResult result;
    for (int sweep = 1; sweep <= options.max_sweeps; ++sweep) {
        using SweepClock = std::chrono::steady_clock;
        const auto sweep_wall_start = SweepClock::now();
        auto elapsed_ms = [](const SweepClock::time_point& start) {
            return std::chrono::duration<double, std::milli>(
                SweepClock::now() - start).count();
        };
        double split_index_ms = 0.0;
        double partition_ms = 0.0;
        double sector_refinement_ms = 0.0;
        auto phase_start = SweepClock::now();
        const WholeTreeSplitIndex starting_index =
            build_whole_tree_split_index(cpu_tree_);
        split_index_ms += elapsed_ms(phase_start);
        SplitSet starting_splits;
        for (const auto& [fingerprint, endpoints] : starting_index.endpoints) {
            (void)endpoints;
            starting_splits.insert(fingerprint);
        }
        phase_start = SweepClock::now();
        std::vector<SplitSet> cores = partition_owned_edge_cores(
            starting_index, starting_splits, options.core_edges);
        partition_ms += elapsed_ms(phase_start);
        SplitSet covered;
        int sweep_moves = 0;
        std::cout << "MLIPPER whole-tree NNI sweep start: sweep=" << sweep
                  << " sectors=" << cores.size()
                  << " internal_edges=" << starting_splits.size()
                  << "\n";
        for (size_t sector = 0; sector < cores.size();) {
            phase_start = SweepClock::now();
            const WholeTreeSplitIndex current_index =
                build_whole_tree_split_index(cpu_tree_);
            split_index_ms += elapsed_ms(phase_start);
            SplitSet owned;
            for (const SplitFingerprint& fingerprint : cores[sector]) {
                if (!covered.count(fingerprint) &&
                    current_index.endpoints.count(fingerprint)) {
                    owned.insert(fingerprint);
                }
            }
            if (owned.empty()) {
                ++sector;
                continue;
            }
            std::vector<NNIOwnedSplit> owned_split_values;
            owned_split_values.reserve(owned.size());
            for (const SplitFingerprint& fingerprint : owned) {
                owned_split_values.push_back(NNIOwnedSplit{
                    fingerprint.side_size,
                    fingerprint.hash_a,
                    fingerprint.hash_b,
                });
            }
            phase_start = SweepClock::now();
            try {
                runSectorNNI(owned_split_values);
            } catch (const CudaRuntimeError& error) {
                if (error.code() != cudaErrorMemoryAllocation ||
                    owned.size() <= 1) {
                    throw;
                }

                // Reduce only the sector that exceeded available GPU memory.
                releaseDivideAndConquerSubtreeCache();
                model_site_batch_workspace_.release();
                const int primary_device = active_gpu_device_;
                releaseOwnedDeviceTreeMemory();
                MlipperDivideAndConquerParams reset_params;
                reset_params.tip_budget = divide_and_conquer_tip_budget_;
                MlipperGpuConfig reset_gpu;
                reset_gpu.acquire_mode =
                    MlipperGpuAcquireMode::UseCurrentDevice;
                reset_gpu.gpu_id = primary_device;
                initializeDivideAndConquerGPU(reset_params, reset_gpu);
                const int retry_core_edges = std::max(
                    1, static_cast<int>((owned.size() + 1) / 2));
                std::vector<SplitSet> retry_cores =
                    partition_owned_edge_cores(
                        current_index, owned, retry_core_edges);
                if (retry_cores.size() <= 1) {
                    throw;
                }
                std::cout << "MLIPPER whole-tree NNI sector split: sweep="
                          << sweep << " sector=" << (sector + 1)
                          << " owned=" << owned.size()
                          << " retry_core_edges=" << retry_core_edges
                          << " replacement_sectors=" << retry_cores.size()
                          << " reason=device-out-of-memory"
                          << "\n";
                cores.erase(cores.begin() + static_cast<std::ptrdiff_t>(sector));
                cores.insert(
                    cores.begin() + static_cast<std::ptrdiff_t>(sector),
                    std::make_move_iterator(retry_cores.begin()),
                    std::make_move_iterator(retry_cores.end()));
                sector_refinement_ms += elapsed_ms(phase_start);
                continue;
            }
            sector_refinement_ms += elapsed_ms(phase_start);
            const int accepted = last_topology_refinement_accepted_move_count_;
            sweep_moves += accepted;
            result.accepted_moves += accepted;
            ++result.sectors;
            covered.insert(owned.begin(), owned.end());
            std::cout << "MLIPPER whole-tree NNI sector: sweep=" << sweep
                      << " sector=" << (sector + 1)
                      << " owned=" << owned.size()
                      << " covered=" << covered.size()
                      << " moves=" << accepted << "\n";
            // The extracted CLVs and scoring buffers are sector-local and
            // invalid after accepted topology moves have been installed.
            releaseDivideAndConquerSubtreeCache();
            model_site_batch_workspace_.release();
            ++sector;
        }
        phase_start = SweepClock::now();
        const WholeTreeSplitIndex final_index =
            build_whole_tree_split_index(cpu_tree_);
        split_index_ms += elapsed_ms(phase_start);
        size_t missing = 0;
        for (const SplitFingerprint& fingerprint : starting_splits) {
            if (final_index.endpoints.count(fingerprint) &&
                !covered.count(fingerprint)) ++missing;
        }
        ++result.sweeps;
        result.covered_edges = covered.size();
        result.final_internal_edges = final_index.endpoints.size();
        std::cout << "MLIPPER whole-tree NNI sweep done: sweep=" << sweep
                  << " covered=" << covered.size()
                  << " missing=" << missing
                  << " moves=" << sweep_moves << "\n";
        std::cout << "MLIPPER whole-tree NNI sweep timing: sweep=" << sweep
                  << " split_index_ms=" << split_index_ms
                  << " partition_ms=" << partition_ms
                  << " sector_refinement_ms=" << sector_refinement_ms
                  << " wall_ms=" << elapsed_ms(sweep_wall_start) << "\n";
        if (missing != 0) {
            throw std::runtime_error(
                "whole-tree NNI sweep ended with uncovered surviving edges");
        }
        if (sweep_moves == 0) {
            result.converged = true;
            break;
        }
    }
    return result;
}

mlipper::BranchLengthOptimizationResult
MlipperSession::optimizeSectorPartitionedBranchLengths(
    const DivideAndConquerFinalOptimizationOptions& options,
    size_t& sector_count,
    size_t& covered_branch_count)
{
    if (gpu_mode_ != GpuMode::DivideAndConquer || active_gpu_device_ < 0) {
        throw std::runtime_error(
            "sector branch optimization requires initialized D&C GPU state");
    }
    mlipper::BranchLengthOptimizationResult aggregate;
    aggregate.sweeps = options.branch_length_sweeps;
    aggregate.newton_iterations = options.branch_length_newton_iterations;

    const int root_id = cpu_tree_.root_id;
    if (root_id < 0 ||
        root_id >= static_cast<int>(cpu_tree_.nodes.size())) {
        throw std::runtime_error(
            "sector branch optimization requires a valid tree root");
    }
    const int artificial_root_child =
        cpu_tree_.nodes[static_cast<size_t>(root_id)].left;
    const int optimized_root_child =
        cpu_tree_.nodes[static_cast<size_t>(root_id)].right;
    if (artificial_root_child < 0 || optimized_root_child < 0) {
        throw std::runtime_error(
            "sector branch optimization requires a bifurcating root");
    }
    const fp_t combined_root_length = static_cast<fp_t>(
        static_cast<double>(cpu_tree_.nodes[
            static_cast<size_t>(artificial_root_child)]
                .branch_length_to_parent) +
        static_cast<double>(cpu_tree_.nodes[
            static_cast<size_t>(optimized_root_child)]
                .branch_length_to_parent));
    cpu_tree_.nodes[static_cast<size_t>(artificial_root_child)]
        .branch_length_to_parent = fp_t(0);
    cpu_tree_.nodes[static_cast<size_t>(optimized_root_child)]
        .branch_length_to_parent = combined_root_length;

    BranchEdgeSet all_branches;
    for (const TreeNode& node : cpu_tree_.nodes) {
        if (node.parent >= 0 && node.id != artificial_root_child) {
            all_branches.insert(node.id);
        }
    }
    std::vector<BranchEdgeSet> cores = partition_branch_edge_cores(
        cpu_tree_, all_branches, options.core_edges);
    BranchEdgeSet covered;
    sector_count = 0;
    covered_branch_count = 0;
    aggregate.attempted = !cores.empty();
    size_t accepted_sectors = 0;

    for (size_t sector = 0; sector < cores.size(); ++sector) {
        const BranchEdgeSet& owned = cores[sector];
        std::vector<divide_and_conquer::TreeEdgeEndpoints> anchors;
        anchors.reserve(owned.size());
        std::vector<int> ordered_owned(owned.begin(), owned.end());
        std::sort(ordered_owned.begin(), ordered_owned.end());
        for (int global_child : ordered_owned) {
            anchors.push_back(divide_and_conquer::TreeEdgeEndpoints{
                global_child,
                cpu_tree_.nodes[static_cast<size_t>(global_child)].parent,
            });
        }

        divide_and_conquer::PreparedDirectionalSubtree directional =
            divide_and_conquer::
                prepare_anchor_subtree_with_directional_boundaries(
                    cpu_tree_,
                    anchors,
                    std::min(
                        std::max(1, divide_and_conquer_tip_budget_),
                        count_live_tip_nodes(cpu_tree_)),
                    tree_alignment_.names,
                    tree_alignment_.sequences,
                    tree_alignment_.sites,
                    model_config_.states,
                    cpu_eig_,
                    cpu_rate_multipliers_,
                    model_config_.ncat);
        directional.host_packing.pattern_weights =
            cpu_host_packing_.pattern_weights;
        const PlacementQueryBatch empty_queries;
        loadSubtreeWorkspace(
            directional.subtree,
            directional.host_packing,
            empty_queries,
            "MlipperSession::optimizeSectorPartitionedBranchLengths",
            true,
            0);
        installStreamingDirectionalBoundaryMessages(
            directional.boundary_ports);
        UpdateTreeClvsPreservingTipClvs(
            subtree_workspace_.dev,
            directional.subtree,
            directional.host_packing,
            subtree_placement_ops_,
            currentStream());

        std::unordered_map<int, int> local_by_global;
        for (size_t local_id = 0;
             local_id < directional.local_to_global.size();
             ++local_id) {
            const int global_id =
                directional.local_to_global[local_id];
            if (global_id >= 0) {
                local_by_global.emplace(
                    global_id, static_cast<int>(local_id));
            }
        }
        std::vector<int> local_owned_edges;
        local_owned_edges.reserve(ordered_owned.size());
        std::unordered_map<int, int> global_by_local_edge;
        const TreeBuildResult& local_tree =
            directional.subtree;
        for (int global_child : ordered_owned) {
            const int global_parent =
                cpu_tree_.nodes[static_cast<size_t>(global_child)].parent;
            const auto child_it = local_by_global.find(global_child);
            const auto parent_it = local_by_global.find(global_parent);
            if (child_it == local_by_global.end() ||
                parent_it == local_by_global.end()) {
                throw std::runtime_error(
                    "sector branch extraction omitted an owned endpoint");
            }
            int local_edge_child = -1;
            if (local_tree.nodes[static_cast<size_t>(child_it->second)].parent ==
                parent_it->second) {
                local_edge_child = child_it->second;
            } else if (local_tree.nodes[
                           static_cast<size_t>(parent_it->second)].parent ==
                       child_it->second) {
                local_edge_child = parent_it->second;
            }
            if (local_edge_child < 0) {
                throw std::runtime_error(
                    "sector branch extraction did not preserve an owned edge");
            }
            local_owned_edges.push_back(local_edge_child);
            global_by_local_edge.emplace(local_edge_child, global_child);
        }

        mlipper::BranchOptimizationOptions branch_options;
        branch_options.sweeps = options.branch_length_sweeps;
        branch_options.newton_iterations =
            options.branch_length_newton_iterations;
        branch_options.likelihood_tolerance = options.likelihood_tolerance;
        branch_options.update_scheme = mlipper::EdgeUpdateScheme::Jacobi;
        branch_options.clv_retention = mlipper::ClvRetention::PreserveTips;
        branch_options.acceptance_scope =
            mlipper::AcceptanceScope::LocalSubtree;
        const mlipper::BranchLengthOptimizationResult local_result =
            RunSelectedTreeEdgeJacobiBranchLengthOptimization(
                subtree_workspace_.dev,
                directional.subtree,
                directional.host_packing,
                subtree_placement_ops_,
                cpu_eig_,
                cpu_rate_multipliers_,
                local_owned_edges,
                currentStream(),
                branch_options);
        aggregate.derivative_seconds += local_result.derivative_seconds;
        if (local_result.accepted) {
            ++accepted_sectors;
            for (const auto& [local_child, global_child] :
                 global_by_local_edge) {
                cpu_tree_.nodes[static_cast<size_t>(global_child)]
                    .branch_length_to_parent =
                    directional.subtree.nodes[
                        static_cast<size_t>(local_child)]
                            .branch_length_to_parent;
            }
        }
        covered.insert(owned.begin(), owned.end());
        ++sector_count;
        covered_branch_count = covered.size();
        releaseDivideAndConquerSubtreeCache();
    }

    if (covered.size() != all_branches.size()) {
        throw std::runtime_error(
            "sector branch optimization did not cover every tree branch");
    }
    cpu_host_packing_ = buildCurrentFullTreeHostPacking();
    backbone_tree_newick_ =
        treeio::write_tree_to_newick_string(cpu_tree_);
    std::cout << "MLIPPER D&C branch partition: sectors="
              << sector_count
              << " covered_branches=" << covered_branch_count
              << " accepted_sectors=" << accepted_sectors
              << "\n";
    return aggregate;
}

DivideAndConquerFinalOptimizationResult
MlipperSession::runDivideAndConquerFinalOptimization(
    const DivideAndConquerFinalOptimizationOptions& options)
{
    if (gpu_mode_ != GpuMode::DivideAndConquer || active_gpu_device_ < 0) {
        throw std::runtime_error(
            "D&C final optimization requires initialized D&C GPU state");
    }
    if (options.max_rounds <= 0 || options.model_parameter_rounds <= 0 ||
        options.core_edges <= 0 || options.branch_length_sweeps <= 0 ||
        options.branch_length_newton_iterations <= 0 ||
        options.likelihood_tolerance < 0.0) {
        throw std::invalid_argument(
            "invalid D&C final optimization options");
    }

    const int primary_device = active_gpu_device_;
    MlipperGpuConfig reuse_gpu;
    reuse_gpu.acquire_mode = MlipperGpuAcquireMode::UseCurrentDevice;
    reuse_gpu.gpu_id = primary_device;
    MlipperDivideAndConquerParams dnc_params;
    dnc_params.tip_budget = divide_and_conquer_tip_budget_;

    DivideAndConquerFinalOptimizationResult result;
    initializeSiteBatchedModelOptimization(reuse_gpu);
    model_optimization::ModelOptimizationBackend backend(
        modelOptimizationContext());
    result.optimization.initial_log_likelihood =
        backend.evaluate();
    result.optimization.final_log_likelihood =
        result.optimization.initial_log_likelihood;

    for (int round = 0; round < options.max_rounds; ++round) {
        const double round_before =
            result.optimization.final_log_likelihood;
        FinalModelOptimizationOptions model_options;
        model_options.max_rounds = options.model_parameter_rounds;
        model_options.likelihood_tolerance = options.likelihood_tolerance;
        model_options.optimize_model_parameters = true;
        model_options.optimize_branch_lengths = false;
        model_options.backend =
            GlobalOptimizationBackend::SiteBatchedSequential;
        FinalModelOptimizationResult model_result =
            runFinalModelOptimization(model_options);
        result.optimization.frequency_updates.insert(
            result.optimization.frequency_updates.end(),
            std::make_move_iterator(model_result.frequency_updates.begin()),
            std::make_move_iterator(model_result.frequency_updates.end()));
        result.optimization.gtr_rate_updates.insert(
            result.optimization.gtr_rate_updates.end(),
            std::make_move_iterator(model_result.gtr_rate_updates.begin()),
            std::make_move_iterator(model_result.gtr_rate_updates.end()));
        result.optimization.alpha_updates.insert(
            result.optimization.alpha_updates.end(),
            std::make_move_iterator(model_result.alpha_updates.begin()),
            std::make_move_iterator(model_result.alpha_updates.end()));

        const double branch_before = model_result.final_log_likelihood;
        const TreeBuildResult tree_before_branches = cpu_tree_;
        initializeDivideAndConquerGPU(dnc_params, reuse_gpu);
        size_t round_sectors = 0;
        size_t round_covered = 0;
        mlipper::BranchLengthOptimizationResult branch_result =
            optimizeSectorPartitionedBranchLengths(
                options, round_sectors, round_covered);
        result.branch_sectors += round_sectors;
        result.covered_branches = round_covered;

        initializeSiteBatchedModelOptimization(reuse_gpu);
        const double branch_after = backend.evaluate();
        branch_result.log_likelihood_before = branch_before;
        branch_result.log_likelihood_after = branch_after;
        branch_result.accepted = std::isfinite(branch_after) &&
            branch_after >= branch_before -
                optimization::branch_lengths::
                    kFullTreeAuditDecreaseTolerance;
        if (!branch_result.accepted) {
            cpu_tree_ = tree_before_branches;
            cpu_host_packing_ = buildCurrentFullTreeHostPacking();
            backbone_tree_newick_ =
                treeio::write_tree_to_newick_string(cpu_tree_);
            initializeSiteBatchedModelOptimization(reuse_gpu);
            branch_result.log_likelihood_after =
                backend.evaluate();
        }
        result.optimization.branch_updates.push_back(branch_result);
        result.optimization.final_log_likelihood =
            branch_result.log_likelihood_after;
        result.optimization.rounds = round + 1;

        const double improvement =
            result.optimization.final_log_likelihood - round_before;
        std::cout << "MLIPPER D&C final optimization round: round="
                  << (round + 1)
                  << " before=" << round_before
                  << " after="
                  << result.optimization.final_log_likelihood
                  << " improvement=" << improvement
                  << " sectors=" << round_sectors
                  << " covered_branches=" << round_covered
                  << "\n";
        if (improvement <
            -optimization::branch_lengths::
                kFullTreeAuditDecreaseTolerance) {
            throw std::runtime_error(
                "D&C final optimization decreased full-tree likelihood");
        }
        if (improvement <= options.likelihood_tolerance) {
            result.optimization.converged = true;
            break;
        }
    }
    backbone_tree_newick_ =
        treeio::write_tree_to_newick_string(cpu_tree_);
    return result;
}

void MlipperSession::writeTree(
    const std::string& output_path,
    double collapse_internal_epsilon) {
    if (has_initialized_cpu_state_ && !cpu_tree_.nodes.empty()) {
        treeio::write_tree_to_newick_file(
            cpu_tree_,
            output_path,
            collapse_internal_epsilon);
        return;
    }
    if (!has_backbone_tree_) {
        throw std::runtime_error(
            "writeTree: no loaded tree is available for export.");
    }

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error(
            "writeTree: failed to open output path '" + output_path + "'.");
    }
    out << backbone_tree_newick_;
    if (backbone_tree_newick_.empty() || backbone_tree_newick_.back() != ';') {
        out << ';';
    }
}

void MlipperSession::writeJplace(
    const std::string& output_path,
    const std::string& invocation,
    const PlacementBatchResult& placements) const
{
    if (placements.ranked_placements.empty() ||
        placements.query_names.empty()) {
        throw std::runtime_error(
            "writeJplace: placement results must not be empty.");
    }
    if (placements.ranked_placements.size() !=
        placements.query_names.size()) {
        throw std::runtime_error(
            "writeJplace: placement and query-name counts differ.");
    }
    jplaceio::write_jplace(
        output_path,
        cpu_tree_,
        placements.ranked_placements,
        placements.query_names,
        invocation);
}

void MlipperSession::initializeGPU(
    const MlipperPlacementParams& params,
    const MlipperGpuConfig& gpu_config)
{
    initializeGPUWorkspace(params, gpu_config, nullptr, false);
}

void MlipperSession::initializeGPUWorkspace(
    const MlipperPlacementParams& params,
    const MlipperGpuConfig& gpu_config,
    const PlacementQueryBatch* queries,
    bool prepare_site_batched_model)
{
    if (!has_initialized_cpu_state_) {
        throw std::runtime_error(
            "initializeGPU: initializeCPU must be called before GPU initialization.");
    }
    if (!std::isfinite(params.accumulated_lwr_threshold) ||
        params.accumulated_lwr_threshold < 0.0 ||
        params.accumulated_lwr_threshold > 1.0) {
        throw std::runtime_error(
            "accumulated LWR threshold must be zero or in (0, 1]");
    }

    if (cpu_tree_.nodes.size() >
            static_cast<size_t>(std::numeric_limits<int>::max()) ||
        cpu_host_packing_.tip_node_ids.size() >
            static_cast<size_t>(std::numeric_limits<int>::max()) ||
        (queries != nullptr && queries->size() >
            static_cast<size_t>(std::numeric_limits<int>::max())) ||
        (has_query_alignment_ && query_alignment_.names.size() >
            static_cast<size_t>(std::numeric_limits<int>::max()))) {
        throw std::overflow_error(
            "initializeGPU: input dimensions exceed DeviceTree integer limits");
    }

    releaseOwnedDeviceTreeMemory();
    DeviceTree& device_tree = owned_device_tree_;

    try {
        const int query_count =
            (queries == nullptr) ? 0 : static_cast<int>(queries->size());
        const int loaded_query_count =
            has_query_alignment_
                ? static_cast<int>(query_alignment_.names.size())
                : 0;
        const int reserve_query_count =
            params.insertion_capacity > 0
                ? params.insertion_capacity
                : ((query_count > 0) ? query_count : loaded_query_count);
        const int reserve_inserts =
            params.commit_to_tree ? std::max(reserve_query_count, 1) : 0;

        model_site_batch_size_ = prepare_site_batched_model
            ? choose_global_optimization_site_batch(
                cpu_tree_.nodes.size(),
                cpu_host_packing_.tip_node_ids.size(),
                tree_alignment_.sites,
                static_cast<size_t>(model_config_.ncat),
                static_cast<size_t>(model_config_.states),
                model_config_.per_rate_scaling,
                0,
                model_optimization::OptimizationWorkspaceKind::
                    GlobalParameters).batch_sites
            : 0;
        // Site-batched model mode keeps model data resident but intentionally
        // omits full-tree directional CLVs. Placement and topology refinement
        // require the resident-tree mode instead.
        const bool allocate_site_batched_only = model_site_batch_size_ > 0 &&
            model_site_batch_size_ < tree_alignment_.sites &&
            query_count == 0 && !params.commit_to_tree && !params.local_spr;
        gpu_mode_ = allocate_site_batched_only
            ? GpuMode::SiteBatchedModel
            : GpuMode::ResidentTree;
        const size_t estimated_active_sites =
            allocate_site_batched_only
                ? model_site_batch_size_
                : tree_alignment_.sites;
        int estimated_process_memory_mb =
            gpu::estimate_mlipper_gpu_process_memory_mb(
                static_cast<int>(cpu_tree_.nodes.size()),
                static_cast<int>(cpu_host_packing_.tip_node_ids.size()),
                reserve_query_count,
                estimated_active_sites,
                model_config_.states,
                model_config_.ncat,
                model_config_.per_rate_scaling,
                params.commit_to_tree);
        if (params.local_spr) {
            estimated_process_memory_mb =
                gpu::add_local_spr_gpu_process_memory_mb(
                    estimated_process_memory_mb);
        }

        const int selected_device =
            prepareGpuDevice(gpu_config, estimated_process_memory_mb);
        allocate_device_tree_on_current_gpu(
            owned_device_tree_,
            cpu_tree_,
            cpu_host_packing_,
            cpu_eig_,
            cpu_rate_weights_,
            cpu_rate_multipliers_,
            cpu_pi_,
            tree_alignment_.sites,
            model_config_.states,
            model_config_.ncat,
            model_config_.per_rate_scaling,
            queries,
            params.commit_to_tree,
            reserve_inserts,
            !allocate_site_batched_only);
        int expected_device = -1;
        require_matching_device(
            selected_device,
            expected_device,
            "initializeGPU",
            "selected GPU device does not match DeviceTree device.");
        require_matching_device(
            device_tree.device_id,
            expected_device,
            "initializeGPU",
            "allocated DeviceTree device does not match selected GPU.");

        ensureGpuResources();
        placement_ops_.tuning.export_accumulated_lwr_threshold =
            params.accumulated_lwr_threshold;
        if (!allocate_site_batched_only) {
            UpdateTreeClvs(
                device_tree, cpu_tree_, cpu_host_packing_,
                placement_ops_, currentStream());
        }

    } catch (...) {
        try {
            releaseOwnedDeviceTreeMemory();
        } catch (...) {
        }
        throw;
    }
}

void MlipperSession::initializeSiteBatchedModelOptimization(
    const MlipperGpuConfig& gpu_config)
{
    initializeGPUWorkspace(
        MlipperPlacementParams{}, gpu_config, nullptr, true);
    enterGlobalParameterOptimizationMode(
        false, GlobalOptimizationBackend::SiteBatchedSequential);
}

void MlipperSession::initializeDivideAndConquerGPU(
    const MlipperDivideAndConquerParams& params,
    const MlipperGpuConfig& gpu_config)
{
    initializeDivideAndConquerGPUImpl(params, gpu_config, nullptr);
}

void MlipperSession::initializeDivideAndConquerGPU(
    const MlipperDivideAndConquerParams& params,
    gpu::DeviceReservation reservation)
{
    const MlipperGpuConfig use_current_device{};
    initializeDivideAndConquerGPUImpl(
        params,
        use_current_device,
        &reservation);
}

void MlipperSession::initializeDivideAndConquerGPUImpl(
    const MlipperDivideAndConquerParams& params,
    const MlipperGpuConfig& gpu_config,
    gpu::DeviceReservation* transferred_reservation)
{
    if (!has_initialized_cpu_state_) {
        throw std::runtime_error(
            "initializeDivideAndConquerGPU: initializeCPU must be called before GPU initialization.");
    }

    releaseOwnedDeviceTreeMemory();

    try {
        const int subtree_tip_budget =
            std::max(1, params.tip_budget);
        const int estimated_process_memory_mb =
            gpu::estimate_divide_and_conquer_gpu_process_memory_mb(
                subtree_tip_budget,
                tree_alignment_.sites,
                model_config_.states,
                model_config_.ncat,
                model_config_.per_rate_scaling);

        prepareGpuDevice(
            gpu_config,
            estimated_process_memory_mb,
            transferred_reservation);
        ensureGpuResources();

        // D&C materializes only one bounded subtree at a time, so the session
        // owns no full-tree DeviceTree in this mode.
        placement_ops_ = PlacementOpBuffer{};
        subtree_placement_ops_ = PlacementOpBuffer{};
        gpu_mode_ = GpuMode::DivideAndConquer;
        divide_and_conquer_tip_budget_ = subtree_tip_budget;
        model_site_batch_size_ = 0;
        owned_device_tree_.reset();
    } catch (...) {
        try {
            releaseOwnedDeviceTreeMemory();
        } catch (...) {
        }
        throw;
    }
}

void MlipperSession::releaseOwnedDeviceTreeMemory()
{
    releasePlacementOps();
    releaseSubtreeWorkspace();
    releaseTopologyRefinementWorkspaces();
    model_site_batch_workspace_.release();
    owned_device_tree_.reset();
    gpu_mode_ = GpuMode::None;
    model_site_batch_size_ = 0;
}

void MlipperSession::releasePlacementWorkspace()
{
    DeviceTree* device_tree =
        owned_device_tree_.d_tipchars != nullptr
            ? &owned_device_tree_
            : nullptr;
    if (device_tree == nullptr) {
        releasePlacementOps();
        releaseSubtreeWorkspace();
        releaseTopologyRefinementWorkspaces();
        clearPendingPlacement();
        return;
    }

    ensure_device_tree_current_device(
        *device_tree,
        "MlipperSession::releasePlacementWorkspace");
    CUDA_CHECK(cudaStreamSynchronize(currentStream()));

    releasePlacementOps();
    releaseSubtreeWorkspace();
    releaseTopologyRefinementWorkspaces();

    // d_clv_up owns one contiguous allocation containing both the upward and
    // downward resident pools.  The global-parameter optimizer constructs its
    // own site-batched upward CLV, so release the complete resident allocation.
    fp_t* const resident_clv_allocation = device_tree->d_clv_up;
    if (device_tree->d_edge_midpoint_clv == resident_clv_allocation) {
        device_tree->d_edge_midpoint_clv = nullptr;
    }
    if (device_tree->d_edge_outside_clv == resident_clv_allocation) {
        device_tree->d_edge_outside_clv = nullptr;
    }
    cuda_free_and_clear(device_tree->d_clv_up);
    device_tree->d_clv_down = nullptr;
    if (device_tree->d_edge_outside_clv == device_tree->d_edge_midpoint_clv) {
        device_tree->d_edge_outside_clv = nullptr;
    }
    cuda_free_and_clear(device_tree->d_edge_midpoint_clv);
    cuda_free_and_clear(device_tree->d_edge_outside_clv);

    // The four resident scaler views share d_scaler_storage.
    cuda_free_and_clear(device_tree->d_scaler_storage);
    device_tree->d_site_scaler_up = nullptr;
    device_tree->d_site_scaler_down = nullptr;
    device_tree->d_edge_midpoint_scaler = nullptr;
    device_tree->d_edge_outside_scaler = nullptr;

    cuda_free_and_clear(device_tree->d_query_chars);
    cuda_free_and_clear(device_tree->d_query_clv);
    cuda_free_and_clear(device_tree->d_query_pmat);
    device_tree->placement_queries = 0;
    device_tree->query_capacity = 0;
    clearPendingPlacement();
}

} // namespace mlipper
