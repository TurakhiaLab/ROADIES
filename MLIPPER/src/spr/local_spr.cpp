#include "local_spr.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "io/tree_newick.hpp"
#include "placement/placement.cuh"
#include "util/mlipper_util.h"

namespace mltreeio = mlipper::treeio;

namespace {

// Local SPR uses a few small records to move between three stages:
// 1. define the local search region around inserted queries,
// 2. score hypothetical prune/regraft edits,
// 3. validate and commit the best edits back to the main tree.

struct PruneInfo {
    int pruned_id = -1;
    int free_internal_id = -1;
    int sibling_id = -1;
    int grandparent_id = -1;
    fp_t pruned_branch_length = fp_t(0);
};

struct LocalSprInsertionAnchor {
    std::string query_name;
    int tip_id = -1;
    int anchor_id = -1;
};

struct LocalSprRepairUnit {
    int unit_id = -1;
    std::vector<std::string> query_names;
    std::vector<int> anchor_ids;
    std::vector<int> anchor_indices;
    std::vector<char> envelope_mask;
    std::vector<int> envelope_nodes;
};

struct LocalSprCandidateMove {
    int repair_unit_id = -1;
    int prune_root_id = -1;
    int regraft_child_id = -1;
    int regraft_parent_id = -1;
    int old_parent_id = -1;
    double approx_gain = -std::numeric_limits<double>::infinity();
    double pendant_length = 0.0;
    double proximal_length = 0.0;
    std::vector<int> subtree_nodes;
    std::vector<int> regraft_path_nodes;
};

struct LocalSprSearchSummary {
    size_t unit_count = 0;
    size_t enumerated_candidates = 0;
    size_t retained_candidates = 0;
    size_t selected_candidates = 0;
};

struct LocalSprPruneRootWorkItem {
    int prune_root_id = -1;
    int skeleton_distance = std::numeric_limits<int>::max();
    std::vector<int> subtree_nodes;
    std::vector<int> legal_inner_candidate_edges;
};

struct LocalSprPlacementPassResult {
    PlacementResult placement_result;
    std::vector<NodeOpInfo> required_downward_update_ops;
};

struct LocalSprScoringWorkspace {
    DeviceTree scoring_dev{};
    PlacementOpBuffer tree_ops{};
    PlacementOpBuffer candidate_ops{};
    int previous_main_pmat_node = -1;
};

struct LocalSprEvalWorkspace {
    BuildToGpuResult res{};
    PlacementOpBuffer ops{};
    bool initialized = false;
};

struct IntDisjointSet {
    std::vector<int> parent;
    std::vector<int> rank;

    explicit IntDisjointSet(int n) : parent((size_t)n), rank((size_t)n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(int x) {
        if (parent[(size_t)x] != x) {
            parent[(size_t)x] = find(parent[(size_t)x]);
        }
        return parent[(size_t)x];
    }

    void unite(int a, int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) return;
        if (rank[(size_t)ra] < rank[(size_t)rb]) std::swap(ra, rb);
        parent[(size_t)rb] = ra;
        if (rank[(size_t)ra] == rank[(size_t)rb]) {
            ++rank[(size_t)ra];
        }
    }
};

// Basic tree maintenance and assertion helpers.

void rebuild_traversals(TreeBuildResult& T) {
    T.preorder.clear();
    T.postorder.clear();
    if (T.root_id < 0 || T.root_id >= static_cast<int>(T.nodes.size())) return;

    std::vector<int> stack;
    stack.push_back(T.root_id);
    while (!stack.empty()) {
        const int id = stack.back();
        stack.pop_back();
        T.preorder.push_back(id);
        const TreeNode& node = T.nodes[static_cast<size_t>(id)];
        if (!node.is_tip) {
            if (node.right >= 0) stack.push_back(node.right);
            if (node.left >= 0) stack.push_back(node.left);
        }
    }

    std::vector<std::pair<int, bool>> st;
    st.emplace_back(T.root_id, false);
    while (!st.empty()) {
        const auto cur = st.back();
        st.pop_back();
        const int id = cur.first;
        const bool visited = cur.second;
        if (id < 0) continue;
        const TreeNode& node = T.nodes[static_cast<size_t>(id)];
        if (visited) {
            T.postorder.push_back(id);
        } else {
            st.emplace_back(id, true);
            if (!node.is_tip) {
                if (node.right >= 0) st.emplace_back(node.right, false);
                if (node.left >= 0) st.emplace_back(node.left, false);
            }
        }
    }
}

[[noreturn]] void local_spr_fail(const std::string& message) {
    throw std::runtime_error("Local SPR subtree assertion failed: " + message);
}

void local_spr_assert(bool condition, const std::string& message) {
    if (!condition) {
        local_spr_fail(message);
    }
}

// Topology and envelope helpers used to describe local subtrees.

void collect_subtree_node_ids(
    const TreeBuildResult& tree,
    int node_id,
    std::vector<int>& node_ids)
{
    if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) return;
    node_ids.clear();
    std::vector<int> stack;
    stack.push_back(node_id);
    while (!stack.empty()) {
        const int cur = stack.back();
        stack.pop_back();
        if (cur < 0 || cur >= static_cast<int>(tree.nodes.size())) continue;
        node_ids.push_back(cur);
        const TreeNode& node = tree.nodes[(size_t)cur];
        if (node.right >= 0) stack.push_back(node.right);
        if (node.left >= 0) stack.push_back(node.left);
    }
}

bool subtree_fully_inside_mask(
    const TreeBuildResult& tree,
    int node_id,
    const std::vector<char>& mask,
    std::vector<int>* node_ids_out = nullptr)
{
    if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) return false;
    if (mask.empty()) return false;

    std::vector<int> node_ids;
    std::vector<int>* sink = node_ids_out ? node_ids_out : &node_ids;
    collect_subtree_node_ids(tree, node_id, *sink);
    for (int cur : *sink) {
        if (cur < 0 || cur >= static_cast<int>(mask.size()) || !mask[(size_t)cur]) {
            return false;
        }
    }
    return true;
}

// Build the minimal downward-update op lists needed to rescore a local edit on GPU.

void rebuild_host_topology_from_tree_local(
    const TreeBuildResult& tree,
    HostPacking& host)
{
    host.postorder = tree.postorder;
    host.preorder = tree.preorder;
    host.parent.resize(tree.nodes.size(), -1);
    host.left.resize(tree.nodes.size(), -1);
    host.right.resize(tree.nodes.size(), -1);
    host.is_tip.resize(tree.nodes.size(), 0);
    host.blen.resize(tree.nodes.size(), fp_t(0));
    for (size_t node_idx = 0; node_idx < tree.nodes.size(); ++node_idx) {
        const TreeNode& node = tree.nodes[node_idx];
        host.parent[node_idx] = node.parent;
        host.left[node_idx] = node.left;
        host.right[node_idx] = node.right;
        host.is_tip[node_idx] = node.is_tip ? 1 : 0;
        host.blen[node_idx] = node.branch_length_to_parent;
    }
}

static inline void append_selected_downward_op(
    std::vector<NodeOpInfo>& ops,
    int parent_id,
    int left_id,
    int right_id,
    bool left_is_tip,
    bool right_is_tip,
    const std::vector<int>& node_to_tip,
    uint8_t dir_tag)
{
    NodeOpInfo op{};
    op.parent_id = parent_id;
    op.left_id = left_id;
    op.right_id = right_id;
    op.left_tip_index = left_is_tip ? node_to_tip[static_cast<size_t>(left_id)] : -1;
    op.right_tip_index = right_is_tip ? node_to_tip[static_cast<size_t>(right_id)] : -1;
    op.clv_pool = static_cast<uint8_t>(CLV_POOL_DOWN);
    op.dir_tag = dir_tag;

    const bool target_is_left = (dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const bool target_is_tip = target_is_left ? left_is_tip : right_is_tip;
    const bool sibling_is_tip = target_is_left ? right_is_tip : left_is_tip;
    if (target_is_tip && sibling_is_tip) {
        op.op_type = static_cast<int>(OP_DOWN_TIP_TIP);
    } else if (target_is_tip) {
        op.op_type = static_cast<int>(OP_DOWN_TIP_INNER);
    } else if (sibling_is_tip) {
        op.op_type = static_cast<int>(OP_DOWN_INNER_TIP);
    } else {
        op.op_type = static_cast<int>(OP_DOWN_INNER_INNER);
    }
    ops.push_back(op);
}

std::vector<int> build_local_spr_node_to_tip(
    const TreeBuildResult& tree,
    const HostPacking& host)
{
    std::vector<int> node_to_tip(tree.nodes.size(), -1);
    for (int tip_idx = 0; tip_idx < static_cast<int>(host.tip_node_ids.size()); ++tip_idx) {
        const int node_id = host.tip_node_ids[static_cast<size_t>(tip_idx)];
        if (node_id >= 0 && node_id < static_cast<int>(node_to_tip.size())) {
            node_to_tip[static_cast<size_t>(node_id)] = tip_idx;
        }
    }
    return node_to_tip;
}

void build_selected_downward_ops(
    const TreeBuildResult& tree,
    const std::vector<int>& node_to_tip,
    const std::vector<int>& target_child_ids,
    std::vector<NodeOpInfo>& host_ops)
{
    host_ops.clear();
    host_ops.reserve(target_child_ids.size());
    std::unordered_set<int> seen;
    for (int target_child_id : target_child_ids) {
        if (!seen.insert(target_child_id).second) continue;
        if (target_child_id < 0 || target_child_id >= static_cast<int>(tree.nodes.size())) continue;
        const int parent_id = tree.nodes[static_cast<size_t>(target_child_id)].parent;
        if (parent_id < 0 || parent_id >= static_cast<int>(tree.nodes.size())) continue;
        const TreeNode& parent = tree.nodes[static_cast<size_t>(parent_id)];
        const int left_id = parent.left;
        const int right_id = parent.right;
        if (left_id < 0 || right_id < 0) continue;

        const bool left_is_tip = tree.nodes[static_cast<size_t>(left_id)].is_tip;
        const bool right_is_tip = tree.nodes[static_cast<size_t>(right_id)].is_tip;
        if (left_id == target_child_id) {
            append_selected_downward_op(
                host_ops,
                parent_id,
                left_id,
                right_id,
                left_is_tip,
                right_is_tip,
                node_to_tip,
                static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
        } else if (right_id == target_child_id) {
            append_selected_downward_op(
                host_ops,
                parent_id,
                left_id,
                right_id,
                left_is_tip,
                right_is_tip,
                node_to_tip,
                static_cast<uint8_t>(CLV_DIR_DOWN_RIGHT));
        } else {
            local_spr_fail("target child does not hang under its recorded parent");
        }
    }
}

void build_required_downward_update_ops(
    const TreeBuildResult& tree,
    const std::vector<int>& node_to_tip,
    const std::vector<int>& target_child_ids,
    std::vector<NodeOpInfo>& host_ops)
{
    std::vector<char> closure_mask(tree.nodes.size(), 0);
    for (int target_child_id : target_child_ids) {
        if (target_child_id < 0 || target_child_id >= static_cast<int>(tree.nodes.size())) continue;
        for (int node_id = target_child_id;
             node_id >= 0 && node_id < static_cast<int>(tree.nodes.size());
             node_id = tree.nodes[static_cast<size_t>(node_id)].parent) {
            const int parent_id = tree.nodes[static_cast<size_t>(node_id)].parent;
            if (parent_id < 0) break;
            closure_mask[static_cast<size_t>(node_id)] = 1;
        }
    }

    host_ops.clear();
    host_ops.reserve(target_child_ids.size() * 4);
    for (int node_id : tree.preorder) {
        if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) continue;
        if (!closure_mask[static_cast<size_t>(node_id)]) continue;

        const int parent_id = tree.nodes[static_cast<size_t>(node_id)].parent;
        if (parent_id < 0 || parent_id >= static_cast<int>(tree.nodes.size())) continue;

        const TreeNode& parent = tree.nodes[static_cast<size_t>(parent_id)];
        const int left_id = parent.left;
        const int right_id = parent.right;
        if (left_id < 0 || right_id < 0) continue;

        const bool left_is_tip = tree.nodes[static_cast<size_t>(left_id)].is_tip;
        const bool right_is_tip = tree.nodes[static_cast<size_t>(right_id)].is_tip;
        if (left_id == node_id) {
            append_selected_downward_op(
                host_ops,
                parent_id,
                left_id,
                right_id,
                left_is_tip,
                right_is_tip,
                node_to_tip,
                static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
        } else if (right_id == node_id) {
            append_selected_downward_op(
                host_ops,
                parent_id,
                left_id,
                right_id,
                left_is_tip,
                right_is_tip,
                node_to_tip,
                static_cast<uint8_t>(CLV_DIR_DOWN_RIGHT));
        } else {
            local_spr_fail("closure node does not hang under its recorded parent");
        }
    }
}

// Deduplicate node ops so staged outward expansion only applies genuinely new work.

struct LocalSprNodeOpKey {
    int parent_id = -1;
    int left_id = -1;
    int right_id = -1;
    int left_tip_index = -1;
    int right_tip_index = -1;
    int op_type = OP_TIP_TIP;
    uint8_t clv_pool = static_cast<uint8_t>(CLV_POOL_UP);
    uint8_t dir_tag = static_cast<uint8_t>(CLV_DIR_UP);

    bool operator==(const LocalSprNodeOpKey& other) const {
        return parent_id == other.parent_id &&
               left_id == other.left_id &&
               right_id == other.right_id &&
               left_tip_index == other.left_tip_index &&
               right_tip_index == other.right_tip_index &&
               op_type == other.op_type &&
               clv_pool == other.clv_pool &&
               dir_tag == other.dir_tag;
    }
};

struct LocalSprNodeOpKeyHash {
    size_t operator()(const LocalSprNodeOpKey& key) const {
        size_t h = static_cast<size_t>(key.parent_id + 0x9e3779b9);
        auto mix = [&](size_t value) {
            h ^= value + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        };
        mix(static_cast<size_t>(key.left_id));
        mix(static_cast<size_t>(key.right_id));
        mix(static_cast<size_t>(key.left_tip_index));
        mix(static_cast<size_t>(key.right_tip_index));
        mix(static_cast<size_t>(key.op_type));
        mix(static_cast<size_t>(key.clv_pool));
        mix(static_cast<size_t>(key.dir_tag));
        return h;
    }
};

LocalSprNodeOpKey make_local_spr_node_op_key(const NodeOpInfo& op) {
    LocalSprNodeOpKey key;
    key.parent_id = op.parent_id;
    key.left_id = op.left_id;
    key.right_id = op.right_id;
    key.left_tip_index = op.left_tip_index;
    key.right_tip_index = op.right_tip_index;
    key.op_type = op.op_type;
    key.clv_pool = op.clv_pool;
    key.dir_tag = op.dir_tag;
    return key;
}

std::vector<NodeOpInfo> filter_local_spr_new_ops(
    const std::vector<NodeOpInfo>& ops,
    const std::vector<NodeOpInfo>& already_computed_ops)
{
    std::unordered_set<LocalSprNodeOpKey, LocalSprNodeOpKeyHash> seen;
    seen.reserve(already_computed_ops.size() * 2 + ops.size());
    for (const NodeOpInfo& op : already_computed_ops) {
        seen.insert(make_local_spr_node_op_key(op));
    }

    std::vector<NodeOpInfo> filtered_ops;
    filtered_ops.reserve(ops.size());
    for (const NodeOpInfo& op : ops) {
        const LocalSprNodeOpKey key = make_local_spr_node_op_key(op);
        if (!seen.insert(key).second) {
            continue;
        }
        filtered_ops.push_back(op);
    }
    return filtered_ops;
}

// Host/GPU conversion helpers for evaluating a fully materialized candidate tree.

std::vector<char> build_local_spr_subtree_mask(
    const TreeBuildResult& tree,
    const std::vector<int>& subtree_nodes)
{
    std::vector<char> subtree_mask(tree.nodes.size(), 0);
    for (int node_id : subtree_nodes) {
        if (node_id >= 0 && node_id < static_cast<int>(subtree_mask.size())) {
            subtree_mask[static_cast<size_t>(node_id)] = 1;
        }
    }
    return subtree_mask;
}

std::vector<int> filter_local_spr_candidate_edges(
    const TreeBuildResult& pruned_tree,
    const std::vector<char>& envelope_mask,
    const std::vector<char>& subtree_mask,
    const std::vector<int>& edge_candidates)
{
    std::vector<int> legal_candidate_edges;
    legal_candidate_edges.reserve(edge_candidates.size());
    for (int edge_child : edge_candidates) {
        if (edge_child < 0 ||
            edge_child >= static_cast<int>(pruned_tree.nodes.size())) {
            continue;
        }
        const int edge_parent = pruned_tree.nodes[static_cast<size_t>(edge_child)].parent;
        if (edge_parent < 0) continue;
        if (!envelope_mask[static_cast<size_t>(edge_child)] ||
            !envelope_mask[static_cast<size_t>(edge_parent)]) {
            continue;
        }
        if (subtree_mask[static_cast<size_t>(edge_child)] ||
            subtree_mask[static_cast<size_t>(edge_parent)]) {
            continue;
        }
        legal_candidate_edges.push_back(edge_child);
    }
    return legal_candidate_edges;
}

HostPacking build_local_spr_tree_host_packing(
    const LocalSprBatchRunContext& ctx,
    const TreeBuildResult& tree)
{
    HostPacking host_pack = pack_host_arrays_from_tree_and_msa(
        tree,
        ctx.current_names,
        ctx.current_rows,
        ctx.sites,
        ctx.states);
    host_pack.pattern_weights = ctx.pattern_weights_arg;
    fill_pmats_in_host_packing(
        tree,
        host_pack,
        ctx.res.eig,
        ctx.rate_multipliers,
        ctx.states,
        ctx.rate_cats);
    return host_pack;
}

void release_local_spr_eval_workspace(
    LocalSprEvalWorkspace& workspace,
    cudaStream_t stream)
{
    if (!workspace.initialized) {
        return;
    }
    free_placement_op_buffer(workspace.ops, stream);
    cudaStreamSynchronize(stream);
    free_device_tree(workspace.res.dev);
    workspace = LocalSprEvalWorkspace{};
}

double evaluate_local_spr_tree_loglikelihood(
    const LocalSprBatchRunContext& ctx,
    const TreeBuildResult& candidate_tree,
    LocalSprEvalWorkspace& workspace)
{
    HostPacking candidate_host_pack =
        build_local_spr_tree_host_packing(ctx, candidate_tree);
    if (!workspace.initialized) {
        workspace.res.tree = candidate_tree;
        workspace.res.hostPack = std::move(candidate_host_pack);
        workspace.res.eig = ctx.res.eig;
        workspace.res.queries = PlacementQueryBatch{};
        workspace.res.dev = upload_to_gpu(
            workspace.res.tree,
            workspace.res.hostPack,
            workspace.res.eig,
            ctx.rate_weights,
            ctx.rate_multipliers,
            ctx.pi,
            ctx.sites,
            ctx.states,
            ctx.rate_cats,
            ctx.per_rate_scaling,
            nullptr,
            false);
        workspace.initialized = true;
    } else {
        workspace.res.tree = candidate_tree;
        workspace.res.hostPack = std::move(candidate_host_pack);
        reload_device_tree_live_data(
            workspace.res.dev,
            workspace.res.tree,
            workspace.res.hostPack,
            nullptr,
            ctx.stream);
    }
    if (workspace.res.tree.nodes.empty() || workspace.res.dev.N == 0) {
        throw std::runtime_error("Local SPR eval produced empty tree/device structures.");
    }
    if (workspace.res.tree.root_id < 0) {
        throw std::runtime_error("Local SPR eval produced tree with invalid root_id.");
    }
    UpdateTreeClvs(
        workspace.res.dev,
        workspace.res.tree,
        workspace.res.hostPack,
        workspace.ops,
        ctx.stream);
    return root_likelihood::compute_root_loglikelihood_total(
        workspace.res.dev,
        workspace.res.tree.root_id,
        workspace.res.dev.d_pattern_weights_u,
        nullptr,
        0.0,
        0);
}

// Candidate legality checks ensure a proposed regraft stays inside the local envelope
// and still matches the topology snapshot it was derived from.

void local_spr_assert_candidate_legal(
    const TreeBuildResult& tree,
    const LocalSprRepairUnit& unit,
    int prune_root_id,
    const std::vector<int>& subtree_nodes,
    const std::vector<int>& legal_edges,
    int regraft_child_id)
{
    local_spr_assert(
        prune_root_id >= 0 && prune_root_id < static_cast<int>(tree.nodes.size()),
        "prune root out of range");
    local_spr_assert(
        prune_root_id != tree.root_id,
        "prune root cannot be the tree root");
    local_spr_assert(
        tree.nodes[static_cast<size_t>(prune_root_id)].parent >= 0,
        "prune root must have a parent");

    std::vector<int> subtree_nodes_check;
    local_spr_assert(
        subtree_fully_inside_mask(tree, prune_root_id, unit.envelope_mask, &subtree_nodes_check),
        "committed subtree is not fully contained in its repair envelope");
    local_spr_assert(
        subtree_nodes_check == subtree_nodes,
        "subtree node set drifted between enumeration and commit");

    local_spr_assert(
        regraft_child_id >= 0 && regraft_child_id < static_cast<int>(tree.nodes.size()),
        "regraft child is out of range");
    const int regraft_parent_id = tree.nodes[static_cast<size_t>(regraft_child_id)].parent;
    local_spr_assert(regraft_parent_id >= 0, "regraft edge has no parent");
    local_spr_assert(
        unit.envelope_mask[static_cast<size_t>(regraft_child_id)] &&
        unit.envelope_mask[static_cast<size_t>(regraft_parent_id)],
        "regraft edge escaped the repair envelope");

    std::vector<char> subtree_mask(tree.nodes.size(), 0);
    for (int node_id : subtree_nodes) {
        if (node_id >= 0 && node_id < static_cast<int>(subtree_mask.size())) {
            subtree_mask[static_cast<size_t>(node_id)] = 1;
        }
    }
    local_spr_assert(
        !subtree_mask[static_cast<size_t>(regraft_child_id)] &&
        !subtree_mask[static_cast<size_t>(regraft_parent_id)],
        "self-regraft detected: target edge lies inside pruned subtree");
    local_spr_assert(
        std::find(legal_edges.begin(), legal_edges.end(), regraft_child_id) != legal_edges.end(),
        "selected regraft edge is not legal under the bounded-radius search");
}

// Graph-distance helpers drive the bounded-radius local search around each prune site.

std::vector<int> build_tree_path_nodes(const TreeBuildResult& tree, int start, int end) {
    std::vector<int> path_a;
    std::unordered_map<int, int> pos_a;
    int cur = start;
    while (cur >= 0) {
        pos_a[cur] = static_cast<int>(path_a.size());
        path_a.push_back(cur);
        cur = tree.nodes[(size_t)cur].parent;
    }

    std::vector<int> path_b;
    int lca = -1;
    cur = end;
    while (cur >= 0) {
        auto it = pos_a.find(cur);
        if (it != pos_a.end()) {
            lca = cur;
            break;
        }
        path_b.push_back(cur);
        cur = tree.nodes[(size_t)cur].parent;
    }

    std::vector<int> out;
    if (lca < 0) return out;
    const int lca_pos = pos_a[lca];
    out.insert(out.end(), path_a.begin(), path_a.begin() + lca_pos + 1);
    std::reverse(path_b.begin(), path_b.end());
    out.insert(out.end(), path_b.begin(), path_b.end());
    return out;
}

std::vector<int> bfs_distances(const TreeBuildResult& tree, int start) {
    const int node_count = static_cast<int>(tree.nodes.size());
    std::vector<int> dist((size_t)node_count, -1);
    if (start < 0 || start >= node_count) return dist;
    std::vector<int> queue;
    queue.reserve((size_t)node_count);
    dist[(size_t)start] = 0;
    queue.push_back(start);
    for (size_t qi = 0; qi < queue.size(); ++qi) {
        const int cur = queue[qi];
        const int next_dist = dist[(size_t)cur] + 1;
        const TreeNode& node = tree.nodes[cur];
        const int neighbors[3] = {node.parent, node.left, node.right};
        for (int neighbor : neighbors) {
            if (neighbor < 0 || neighbor >= node_count) continue;
            if (dist[(size_t)neighbor] >= 0) continue;
            dist[(size_t)neighbor] = next_dist;
            queue.push_back(neighbor);
        }
    }
    return dist;
}

std::vector<LocalSprInsertionAnchor> build_local_spr_insertion_anchors(
    const TreeBuildResult& tree,
    const std::vector<std::string>& inserted_names)
{
    std::vector<LocalSprInsertionAnchor> anchors;
    anchors.reserve(inserted_names.size());
    for (const std::string& name : inserted_names) {
        auto tip_it = tree.tip_node_by_name.find(name);
        if (tip_it == tree.tip_node_by_name.end()) {
            throw std::runtime_error("Local SPR could not find inserted tip in tree: " + name);
        }
        const int tip_id = tip_it->second;
        if (tip_id < 0 || tip_id >= static_cast<int>(tree.nodes.size())) {
            throw std::runtime_error("Local SPR tip id out of range for: " + name);
        }
        const int anchor_id = tree.nodes[(size_t)tip_id].parent;
        if (anchor_id < 0 || anchor_id >= static_cast<int>(tree.nodes.size())) {
            throw std::runtime_error("Local SPR inserted tip is missing an anchor parent: " + name);
        }
        anchors.push_back(LocalSprInsertionAnchor{name, tip_id, anchor_id});
    }
    return anchors;
}

std::vector<LocalSprRepairUnit> build_local_spr_repair_units(
    const TreeBuildResult& tree,
    const std::vector<LocalSprInsertionAnchor>& anchors,
    int cluster_threshold,
    int envelope_radius)
{
    std::vector<LocalSprRepairUnit> units;
    const int anchor_count = static_cast<int>(anchors.size());
    if (anchor_count == 0) return units;

    cluster_threshold = std::max(0, cluster_threshold);
    IntDisjointSet dsu(anchor_count);
    std::vector<std::vector<int>> anchor_dist((size_t)anchor_count);
    for (int i = 0; i < anchor_count; ++i) {
        anchor_dist[(size_t)i] = bfs_distances(tree, anchors[(size_t)i].anchor_id);
    }

    for (int i = 0; i < anchor_count; ++i) {
        for (int j = i + 1; j < anchor_count; ++j) {
            const int d = anchor_dist[(size_t)i][(size_t)anchors[(size_t)j].anchor_id];
            if (d >= 0 && d <= cluster_threshold) {
                dsu.unite(i, j);
            }
        }
    }

    std::unordered_map<int, int> root_to_unit;
    root_to_unit.reserve((size_t)anchor_count);
    for (int i = 0; i < anchor_count; ++i) {
        const int root = dsu.find(i);
        auto [it, inserted] = root_to_unit.emplace(root, static_cast<int>(units.size()));
        if (inserted) {
            LocalSprRepairUnit unit;
            unit.unit_id = static_cast<int>(units.size());
            unit.envelope_mask.assign(tree.nodes.size(), 0);
            units.push_back(std::move(unit));
        }
        LocalSprRepairUnit& unit = units[(size_t)it->second];
        unit.query_names.push_back(anchors[(size_t)i].query_name);
        unit.anchor_ids.push_back(anchors[(size_t)i].anchor_id);
        unit.anchor_indices.push_back(i);
    }

    for (LocalSprRepairUnit& unit : units) {
        for (int anchor_idx : unit.anchor_indices) {
            const std::vector<int>& dist = anchor_dist[(size_t)anchor_idx];
            for (int node_id = 0; node_id < static_cast<int>(tree.nodes.size()); ++node_id) {
                if (dist[(size_t)node_id] >= 0 && dist[(size_t)node_id] <= envelope_radius) {
                    unit.envelope_mask[(size_t)node_id] = 1;
                }
            }
        }
        for (int node_id = 0; node_id < static_cast<int>(tree.nodes.size()); ++node_id) {
            if (unit.envelope_mask[(size_t)node_id]) {
                unit.envelope_nodes.push_back(node_id);
            }
        }
    }

    return units;
}

void keep_local_spr_topk(
    std::vector<LocalSprCandidateMove>& topk,
    LocalSprCandidateMove candidate,
    int topk_limit)
{
    topk.push_back(std::move(candidate));
    std::sort(
        topk.begin(),
        topk.end(),
        [](const LocalSprCandidateMove& lhs, const LocalSprCandidateMove& rhs) {
            return lhs.approx_gain > rhs.approx_gain;
        });
    if (static_cast<int>(topk.size()) > topk_limit) {
        topk.resize((size_t)topk_limit);
    }
}

std::vector<int> select_local_spr_seed_edges(
    const PlacementResult& placement_result,
    const TreeBuildResult& tree,
    const std::vector<int>& center_dist,
    double baseline_logL,
    int seed_limit,
    int required_edge_radius)
{
    std::vector<std::pair<double, int>> ranked_edges;
    ranked_edges.reserve(placement_result.top_placements.size());
    for (const PlacementResult::RankedPlacement& placement :
         placement_result.top_placements) {
        const int edge_child_id = placement.target_id;
        if (edge_child_id < 0 || edge_child_id >= static_cast<int>(tree.nodes.size())) {
            continue;
        }
        const int edge_parent_id = tree.nodes[(size_t)edge_child_id].parent;
        if (edge_parent_id < 0) {
            continue;
        }
        if (required_edge_radius >= 0) {
            const int edge_child_dist = center_dist[(size_t)edge_child_id];
            const int edge_parent_dist = center_dist[(size_t)edge_parent_id];
            int edge_radius = -1;
            if (edge_child_dist < 0) {
                edge_radius = edge_parent_dist;
            } else if (edge_parent_dist < 0) {
                edge_radius = edge_child_dist;
            } else {
                edge_radius = std::min(edge_child_dist, edge_parent_dist);
            }
            if (edge_radius != required_edge_radius) {
                continue;
            }
        }
        const double approx_gain = placement.loglikelihood - baseline_logL;
        if (!(approx_gain > 0.0)) {
            continue;
        }
        ranked_edges.emplace_back(approx_gain, edge_child_id);
    }

    std::sort(
        ranked_edges.begin(),
        ranked_edges.end(),
        [](const std::pair<double, int>& lhs, const std::pair<double, int>& rhs) {
            if (lhs.first == rhs.first) return lhs.second < rhs.second;
            return lhs.first > rhs.first;
        });

    std::vector<int> seed_edges;
    seed_edges.reserve(ranked_edges.size());
    std::unordered_set<int> seen;
    for (const auto& [gain, edge_child_id] : ranked_edges) {
        (void)gain;
        if (!seen.insert(edge_child_id).second) continue;
        seed_edges.push_back(edge_child_id);
        if (seed_limit > 0 && static_cast<int>(seed_edges.size()) >= seed_limit) {
            break;
        }
    }
    return seed_edges;
}

std::vector<LocalSprCandidateMove> select_local_spr_candidates(
    const std::vector<LocalSprCandidateMove>& ranked_candidates,
    int node_count)
{
    std::vector<LocalSprCandidateMove> selected;
    std::vector<char> selected_subtree_mask((size_t)std::max(0, node_count), 0);
    std::vector<char> selected_path_mask((size_t)std::max(0, node_count), 0);
    std::unordered_set<int> used_units;

    for (const LocalSprCandidateMove& candidate : ranked_candidates) {
        if (!used_units.insert(candidate.repair_unit_id).second) {
            continue;
        }

        bool conflict = false;
        for (int node_id : candidate.subtree_nodes) {
            if (node_id >= 0 &&
                node_id < node_count &&
                selected_subtree_mask[(size_t)node_id]) {
                conflict = true;
                break;
            }
        }
        if (conflict) {
            used_units.erase(candidate.repair_unit_id);
            continue;
        }
        if (candidate.regraft_child_id >= 0 &&
            candidate.regraft_child_id < node_count &&
            selected_subtree_mask[(size_t)candidate.regraft_child_id]) {
            used_units.erase(candidate.repair_unit_id);
            continue;
        }
        if (candidate.regraft_parent_id >= 0 &&
            candidate.regraft_parent_id < node_count &&
            selected_subtree_mask[(size_t)candidate.regraft_parent_id]) {
            used_units.erase(candidate.repair_unit_id);
            continue;
        }
        for (int node_id : candidate.regraft_path_nodes) {
            if (node_id >= 0 &&
                node_id < node_count &&
                selected_path_mask[(size_t)node_id]) {
                conflict = true;
                break;
            }
        }
        if (conflict) {
            used_units.erase(candidate.repair_unit_id);
            continue;
        }

        selected.push_back(candidate);
        for (int node_id : candidate.subtree_nodes) {
            if (node_id >= 0 && node_id < node_count) {
                selected_subtree_mask[(size_t)node_id] = 1;
            }
        }
        for (int node_id : candidate.regraft_path_nodes) {
            if (node_id >= 0 && node_id < node_count) {
                selected_path_mask[(size_t)node_id] = 1;
            }
        }
    }

    return selected;
}

bool local_spr_candidate_still_legal(
    const TreeBuildResult& tree,
    const LocalSprRepairUnit& unit,
    const LocalSprCandidateMove& candidate,
    std::vector<int>* current_subtree_nodes_out = nullptr)
{
    if (candidate.prune_root_id < 0 ||
        candidate.prune_root_id >= static_cast<int>(tree.nodes.size()) ||
        candidate.prune_root_id == tree.root_id) {
        return false;
    }
    const int current_parent =
        tree.nodes[(size_t)candidate.prune_root_id].parent;
    if (current_parent < 0 ||
        (candidate.old_parent_id >= 0 && current_parent != candidate.old_parent_id)) {
        return false;
    }

    std::vector<int> current_subtree_nodes_storage;
    std::vector<int>* current_subtree_nodes =
        current_subtree_nodes_out != nullptr
            ? current_subtree_nodes_out
            : &current_subtree_nodes_storage;
    if (!subtree_fully_inside_mask(
            tree,
            candidate.prune_root_id,
            unit.envelope_mask,
            current_subtree_nodes)) {
        return false;
    }

    if (candidate.regraft_child_id < 0 ||
        candidate.regraft_child_id >= static_cast<int>(tree.nodes.size())) {
        return false;
    }
    const int current_regraft_parent =
        tree.nodes[(size_t)candidate.regraft_child_id].parent;
    if (current_regraft_parent < 0 ||
        (candidate.regraft_parent_id >= 0 &&
         current_regraft_parent != candidate.regraft_parent_id)) {
        return false;
    }
    if (!unit.envelope_mask[(size_t)candidate.regraft_child_id] ||
        !unit.envelope_mask[(size_t)current_regraft_parent]) {
        return false;
    }

    std::vector<char> current_subtree_mask(tree.nodes.size(), 0);
    for (int node_id : *current_subtree_nodes) {
        if (node_id >= 0 &&
            node_id < static_cast<int>(current_subtree_mask.size())) {
            current_subtree_mask[(size_t)node_id] = 1;
        }
    }
    if (current_subtree_mask[(size_t)candidate.regraft_child_id] ||
        current_subtree_mask[(size_t)current_regraft_parent]) {
        return false;
    }

    return true;
}

// These two routines are the actual topology edit primitives:
// prune_subtree_for_spr removes a local subtree from the current tree,
// regraft_subtree_for_spr inserts it onto a new target edge.

bool prune_subtree_for_spr(TreeBuildResult& tree, int pruned_id, PruneInfo& info) {
    if (pruned_id < 0 || pruned_id >= static_cast<int>(tree.nodes.size())) return false;
    if (pruned_id == tree.root_id) return false;
    const int parent_id = tree.nodes[pruned_id].parent;
    if (parent_id < 0 || parent_id >= static_cast<int>(tree.nodes.size())) return false;
    const TreeNode& parent = tree.nodes[parent_id];
    const int sibling_id = (parent.left == pruned_id) ? parent.right :
                           (parent.right == pruned_id ? parent.left : -1);
    if (sibling_id < 0 || sibling_id >= static_cast<int>(tree.nodes.size())) return false;

    const int grandparent_id = parent.parent;
    info.pruned_id = pruned_id;
    info.free_internal_id = parent_id;
    info.sibling_id = sibling_id;
    info.grandparent_id = grandparent_id;
    info.pruned_branch_length = tree.nodes[pruned_id].branch_length_to_parent;

    tree.nodes[pruned_id].parent = -1;

    if (grandparent_id >= 0) {
        TreeNode& grandparent = tree.nodes[grandparent_id];
        if (grandparent.left == parent_id) {
            grandparent.left = sibling_id;
        } else if (grandparent.right == parent_id) {
            grandparent.right = sibling_id;
        } else {
            return false;
        }
        TreeNode& sibling = tree.nodes[sibling_id];
        sibling.parent = grandparent_id;
        sibling.branch_length_to_parent += parent.branch_length_to_parent;
    } else {
        tree.root_id = sibling_id;
        TreeNode& sibling = tree.nodes[sibling_id];
        sibling.parent = -1;
        sibling.branch_length_to_parent = fp_t(0);
    }

    TreeNode& free_internal = tree.nodes[parent_id];
    free_internal.left = -1;
    free_internal.right = -1;
    free_internal.parent = -1;
    free_internal.is_tip = false;
    free_internal.name.clear();
    free_internal.branch_length_to_parent = fp_t(0);
    rebuild_traversals(tree);
    return true;
}

std::vector<char> build_unit_anchor_skeleton_mask(
    const TreeBuildResult& tree,
    const std::vector<int>& anchor_ids)
{
    std::vector<char> mask(tree.nodes.size(), 0);
    if (anchor_ids.empty()) return mask;

    const int root_anchor = anchor_ids.front();
    if (root_anchor >= 0 && root_anchor < static_cast<int>(tree.nodes.size())) {
        mask[(size_t)root_anchor] = 1;
    }

    for (int anchor_id : anchor_ids) {
        if (anchor_id < 0 || anchor_id >= static_cast<int>(tree.nodes.size())) continue;
        const std::vector<int> path_nodes =
            build_tree_path_nodes(tree, root_anchor, anchor_id);
        for (int node_id : path_nodes) {
            if (node_id < 0 || node_id >= static_cast<int>(mask.size())) continue;
            mask[(size_t)node_id] = 1;
        }
    }
    return mask;
}

std::vector<int> multi_source_bfs_distances(
    const TreeBuildResult& tree,
    const std::vector<char>& source_mask)
{
    const int node_count = static_cast<int>(tree.nodes.size());
    std::vector<int> dist((size_t)node_count, -1);
    if (static_cast<int>(source_mask.size()) != node_count) return dist;

    std::vector<int> queue;
    queue.reserve((size_t)node_count);
    for (int node_id = 0; node_id < node_count; ++node_id) {
        if (!source_mask[(size_t)node_id]) continue;
        dist[(size_t)node_id] = 0;
        queue.push_back(node_id);
    }

    for (size_t qi = 0; qi < queue.size(); ++qi) {
        const int cur = queue[qi];
        const int next_dist = dist[(size_t)cur] + 1;
        const TreeNode& node = tree.nodes[(size_t)cur];
        const int neighbors[3] = {node.parent, node.left, node.right};
        for (int neighbor : neighbors) {
            if (neighbor < 0 || neighbor >= node_count) continue;
            if (dist[(size_t)neighbor] >= 0) continue;
            dist[(size_t)neighbor] = next_dist;
            queue.push_back(neighbor);
        }
    }
    return dist;
}

int edge_distance_from_endpoint(const std::vector<int>& dist, int child_id, int parent_id) {
    if (child_id < 0 || parent_id < 0) return -1;
    const int d_child = dist[(size_t)child_id];
    const int d_parent = dist[(size_t)parent_id];
    if (d_child < 0) return d_parent;
    if (d_parent < 0) return d_child;
    return std::min(d_child, d_parent);
}

std::vector<int> collect_candidate_edges(
    const TreeBuildResult& tree,
    int endpoint_a,
    int endpoint_b,
    int radius,
    int exclude_a,
    int exclude_b)
{
    std::vector<int> edges;
    if (radius < 0) return edges;
    const int node_count = static_cast<int>(tree.nodes.size());
    std::vector<int> dist_a = bfs_distances(tree, endpoint_a);
    std::vector<int> dist_b;
    if (endpoint_b >= 0 && endpoint_b != endpoint_a) {
        dist_b = bfs_distances(tree, endpoint_b);
    }
    for (int node_id = 0; node_id < node_count; ++node_id) {
        const int parent_id = tree.nodes[node_id].parent;
        if (parent_id < 0) continue;
        if (node_id == exclude_a || node_id == exclude_b) continue;
        if (parent_id == exclude_a || parent_id == exclude_b) continue;
        int best = edge_distance_from_endpoint(dist_a, node_id, parent_id);
        if (!dist_b.empty()) {
            const int alt = edge_distance_from_endpoint(dist_b, node_id, parent_id);
            if (best < 0 || (alt >= 0 && alt < best)) best = alt;
        }
        if (best >= 0 && best <= radius) {
            edges.push_back(node_id);
        }
    }
    return edges;
}

std::vector<int> compute_center_distances(
    const TreeBuildResult& tree,
    int endpoint_a,
    int endpoint_b)
{
    std::vector<int> best = bfs_distances(tree, endpoint_a);
    if (endpoint_b >= 0 && endpoint_b != endpoint_a) {
        const std::vector<int> alt = bfs_distances(tree, endpoint_b);
        if (best.size() != alt.size()) {
            throw std::runtime_error("compute_center_distances produced mismatched BFS arrays.");
        }
        for (size_t i = 0; i < best.size(); ++i) {
            if (best[i] < 0) {
                best[i] = alt[i];
            } else if (alt[i] >= 0 && alt[i] < best[i]) {
                best[i] = alt[i];
            }
        }
    }
    return best;
}

int edge_distance_from_center(
    const std::vector<int>& center_dist,
    int child_id,
    int parent_id)
{
    if (child_id < 0 || parent_id < 0) return -1;
    const int d_child = center_dist[(size_t)child_id];
    const int d_parent = center_dist[(size_t)parent_id];
    if (d_child < 0) return d_parent;
    if (d_parent < 0) return d_child;
    return std::min(d_child, d_parent);
}

std::vector<int> collect_outward_edges_from_seed(
    const TreeBuildResult& tree,
    const std::vector<int>& center_dist,
    int seed_edge_child_id,
    int min_radius_exclusive,
    int max_radius,
    int exclude_a,
    int exclude_b)
{
    std::vector<int> edges;
    if (max_radius <= min_radius_exclusive) return edges;
    if (seed_edge_child_id < 0 ||
        seed_edge_child_id >= static_cast<int>(tree.nodes.size())) {
        return edges;
    }

    const int seed_parent_id = tree.nodes[(size_t)seed_edge_child_id].parent;
    if (seed_parent_id < 0 ||
        seed_parent_id >= static_cast<int>(tree.nodes.size())) {
        return edges;
    }

    const int child_dist = center_dist[(size_t)seed_edge_child_id];
    const int parent_dist = center_dist[(size_t)seed_parent_id];
    if (child_dist < 0 || parent_dist < 0 || child_dist == parent_dist) {
        return edges;
    }

    struct FrontierState {
        int node_id = -1;
        int prev_id = -1;
    };

    std::vector<FrontierState> queue;
    queue.reserve(tree.nodes.size());
    std::vector<char> visited(tree.nodes.size(), 0);
    std::vector<char> edge_seen(tree.nodes.size(), 0);

    auto push_start = [&](int node_id, int prev_id) {
        if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) return;
        if (visited[(size_t)node_id]) return;
        visited[(size_t)node_id] = 1;
        queue.push_back(FrontierState{node_id, prev_id});
    };

    if (child_dist > parent_dist) {
        push_start(seed_edge_child_id, seed_parent_id);
    } else {
        push_start(seed_parent_id, seed_edge_child_id);
    }

    auto maybe_record_edge = [&](int node_a, int node_b) {
        int edge_child_id = -1;
        int edge_parent_id = -1;
        if (tree.nodes[(size_t)node_a].parent == node_b) {
            edge_child_id = node_a;
            edge_parent_id = node_b;
        } else if (tree.nodes[(size_t)node_b].parent == node_a) {
            edge_child_id = node_b;
            edge_parent_id = node_a;
        } else {
            return;
        }

        if (edge_child_id == exclude_a || edge_child_id == exclude_b) return;
        if (edge_parent_id == exclude_a || edge_parent_id == exclude_b) return;
        if (edge_seen[(size_t)edge_child_id]) return;

        const int edge_radius =
            edge_distance_from_center(center_dist, edge_child_id, edge_parent_id);
        if (edge_radius <= min_radius_exclusive || edge_radius > max_radius) return;

        edge_seen[(size_t)edge_child_id] = 1;
        edges.push_back(edge_child_id);
    };

    for (size_t qi = 0; qi < queue.size(); ++qi) {
        const FrontierState state = queue[qi];
        const int cur = state.node_id;
        const int cur_dist = center_dist[(size_t)cur];
        const TreeNode& node = tree.nodes[(size_t)cur];
        const int neighbors[3] = {node.parent, node.left, node.right};
        for (int neighbor : neighbors) {
            if (neighbor < 0 || neighbor == state.prev_id) continue;
            if (neighbor >= static_cast<int>(tree.nodes.size())) continue;
            const int neighbor_dist = center_dist[(size_t)neighbor];
            if (neighbor_dist < 0 || neighbor_dist <= cur_dist) continue;

            maybe_record_edge(cur, neighbor);
            if (!visited[(size_t)neighbor]) {
                visited[(size_t)neighbor] = 1;
                queue.push_back(FrontierState{neighbor, cur});
            }
        }
    }

    return edges;
}

void regraft_subtree_for_spr(
    TreeBuildResult& tree,
    const PruneInfo& info,
    int target_child_id,
    bool use_lengths,
    double pendant_length,
    double proximal_length)
{
    const int parent_id = tree.nodes[target_child_id].parent;
    if (parent_id < 0) {
        throw std::runtime_error("SPR regraft target has no parent edge.");
    }
    const int free_internal_id = info.free_internal_id;
    const int pruned_id = info.pruned_id;

    TreeNode& parent = tree.nodes[parent_id];
    if (parent.left == target_child_id) {
        parent.left = free_internal_id;
    } else if (parent.right == target_child_id) {
        parent.right = free_internal_id;
    } else {
        throw std::runtime_error("SPR regraft target not found under its parent.");
    }

    TreeNode& internal = tree.nodes[free_internal_id];
    internal.parent = parent_id;
    internal.left = target_child_id;
    internal.right = pruned_id;
    internal.is_tip = false;
    internal.name.clear();

    TreeNode& target_child = tree.nodes[target_child_id];
    const fp_t total_length = target_child.branch_length_to_parent;
    fp_t proximal = total_length * fp_t(0.5);
    fp_t distal = total_length - proximal;
    if (use_lengths) {
        double prox = 0.0;
        double dist = 0.0;
        normalize_split_branch_lengths(
            static_cast<double>(total_length),
            proximal_length,
            OPT_BRANCH_LEN_MIN,
            prox,
            dist);
        proximal = static_cast<fp_t>(prox);
        distal = static_cast<fp_t>(dist);
    }
    internal.branch_length_to_parent = proximal;
    target_child.branch_length_to_parent = distal;
    target_child.parent = free_internal_id;

    TreeNode& pruned = tree.nodes[pruned_id];
    pruned.parent = free_internal_id;
    pruned.branch_length_to_parent =
        use_lengths
            ? static_cast<fp_t>(sanitize_branch_length(pendant_length))
            : info.pruned_branch_length;
    rebuild_traversals(tree);
}

// The search pass temporarily overrides a few placement-tuning environment variables
// so local SPR scoring uses a predictable, lightweight configuration.

struct EnvSnapshot {
    const char* key = nullptr;
    std::string value;
    bool has_value = false;
};

EnvSnapshot snapshot_env(const char* key) {
    const char* value = std::getenv(key);
    if (value) {
        return EnvSnapshot{key, std::string(value), true};
    }
    return EnvSnapshot{key, std::string(), false};
}

void restore_env(const EnvSnapshot& snap) {
    if (snap.has_value) {
        setenv(snap.key, snap.value.c_str(), 1);
    } else {
        unsetenv(snap.key);
    }
}

struct LocalSprEnvGuard {
    EnvSnapshot full_opt = snapshot_env("MLIPPER_FULL_OPT_PASSES");
    EnvSnapshot export_topk = snapshot_env("MLIPPER_EXPORT_PLACEMENT_TOPK");
    EnvSnapshot double_rerank = snapshot_env("MLIPPER_DOUBLE_RERANK");

    LocalSprEnvGuard(bool local_spr_fast, size_t node_count) {
        setenv("MLIPPER_FULL_OPT_PASSES", local_spr_fast ? "1" : "4", 1);
        setenv(
            "MLIPPER_EXPORT_PLACEMENT_TOPK",
            std::to_string(std::max(1, static_cast<int>(node_count) * 2)).c_str(),
            1);
        setenv("MLIPPER_DOUBLE_RERANK", "0", 1);
    }

    ~LocalSprEnvGuard() {
        restore_env(full_opt);
        restore_env(double_rerank);
        restore_env(export_topk);
    }
};

struct LocalSprSearchContext {
    BuildToGpuResult& res;
    cudaStream_t stream = nullptr;
    const std::vector<unsigned>& pattern_weights_arg;
    const std::vector<double>& rate_weights;
    const std::vector<double>& rate_multipliers;
    const std::vector<double>& pi;
    size_t sites = 0;
    int states = 0;
    int rate_cats = 0;
    bool per_rate_scaling = false;
    int local_spr_radius = 4;
    int local_spr_topk_per_unit = 0;
};

// Repair-unit construction groups nearby inserted queries into one local search region.

std::vector<LocalSprRepairUnit> prepare_local_spr_repair_units(
    const TreeBuildResult& base_tree,
    const std::vector<std::string>& inserted_names,
    int cluster_threshold,
    int radius)
{
    const std::vector<LocalSprInsertionAnchor> anchors =
        build_local_spr_insertion_anchors(base_tree, inserted_names);
    return build_local_spr_repair_units(
        base_tree,
        anchors,
        cluster_threshold,
        radius);
}

std::vector<LocalSprPruneRootWorkItem> build_local_spr_prune_root_work_items(
    const TreeBuildResult& base_tree,
    const LocalSprRepairUnit& unit,
    int inner_search_radius,
    const std::vector<int>& skeleton_distance_by_node)
{
    std::vector<LocalSprPruneRootWorkItem> work_items;
    work_items.reserve(unit.envelope_nodes.size());

    for (int candidate_prune_root_id : unit.envelope_nodes) {
        if (candidate_prune_root_id < 0 ||
            candidate_prune_root_id >= static_cast<int>(base_tree.nodes.size()) ||
            candidate_prune_root_id == base_tree.root_id ||
            base_tree.nodes[(size_t)candidate_prune_root_id].parent < 0) {
            continue;
        }

        std::vector<int> subtree_nodes;
        if (!subtree_fully_inside_mask(
                base_tree,
                candidate_prune_root_id,
                unit.envelope_mask,
                &subtree_nodes)) {
            continue;
        }

        const std::vector<char> subtree_mask =
            build_local_spr_subtree_mask(base_tree, subtree_nodes);

        TreeBuildResult pruned_tree = base_tree;
        PruneInfo prune_info;
        if (!prune_subtree_for_spr(
                pruned_tree,
                candidate_prune_root_id,
                prune_info)) {
            continue;
        }

        std::vector<int> inner_candidate_edges = collect_candidate_edges(
            pruned_tree,
            prune_info.sibling_id,
            prune_info.grandparent_id,
            inner_search_radius,
            prune_info.pruned_id,
            prune_info.free_internal_id);
        std::vector<int> legal_inner_candidate_edges =
            filter_local_spr_candidate_edges(
                pruned_tree,
                unit.envelope_mask,
                subtree_mask,
                inner_candidate_edges);
        if (legal_inner_candidate_edges.empty()) {
            continue;
        }

        int skeleton_distance = std::numeric_limits<int>::max();
        if (candidate_prune_root_id <
            static_cast<int>(skeleton_distance_by_node.size())) {
            const int dist =
                skeleton_distance_by_node[(size_t)candidate_prune_root_id];
            if (dist >= 0) {
                skeleton_distance = dist;
            }
        }

        LocalSprPruneRootWorkItem work_item;
        work_item.prune_root_id = candidate_prune_root_id;
        work_item.skeleton_distance = skeleton_distance;
        work_item.subtree_nodes = std::move(subtree_nodes);
        work_item.legal_inner_candidate_edges =
            std::move(legal_inner_candidate_edges);
        work_items.push_back(std::move(work_item));
    }

    std::sort(
        work_items.begin(),
        work_items.end(),
        [](const LocalSprPruneRootWorkItem& lhs,
           const LocalSprPruneRootWorkItem& rhs) {
            if (lhs.skeleton_distance != rhs.skeleton_distance) {
                return lhs.skeleton_distance < rhs.skeleton_distance;
            }
            if (lhs.legal_inner_candidate_edges.size() !=
                rhs.legal_inner_candidate_edges.size()) {
                return lhs.legal_inner_candidate_edges.size() >
                       rhs.legal_inner_candidate_edges.size();
            }
            if (lhs.subtree_nodes.size() != rhs.subtree_nodes.size()) {
                return lhs.subtree_nodes.size() <
                       rhs.subtree_nodes.size();
            }
            return lhs.prune_root_id < rhs.prune_root_id;
        });

    return work_items;
}

// Workspace lifecycle and placement-pass helpers for the candidate ranking phase.

void release_local_spr_scoring_workspace(
    LocalSprScoringWorkspace& workspace,
    cudaStream_t stream)
{
    free_placement_op_buffer(workspace.candidate_ops, stream);
    free_placement_op_buffer(workspace.tree_ops, stream);
    cudaStreamSynchronize(stream);
    free_device_tree(workspace.scoring_dev);
    workspace = LocalSprScoringWorkspace{};
}

LocalSprPlacementPassResult run_local_spr_placement_pass(
    const LocalSprSearchContext& ctx,
    TreeBuildResult& pruned_tree,
    HostPacking& pruned_host,
    const std::vector<int>& node_to_tip,
    const std::vector<int>& candidate_edges,
    int upward_start_node,
    int prune_root_id,
    LocalSprScoringWorkspace& workspace,
    const std::vector<NodeOpInfo>* already_updated_ops)
{
    LocalSprPlacementPassResult pass_result;
    std::vector<NodeOpInfo> candidate_ops;
    build_selected_downward_ops(
        pruned_tree,
        node_to_tip,
        candidate_edges,
        candidate_ops);
    build_required_downward_update_ops(
        pruned_tree,
        node_to_tip,
        candidate_edges,
        pass_result.required_downward_update_ops);
    local_spr_assert(
        !candidate_ops.empty(),
        "local SPR scoring produced zero candidate ops");
    local_spr_assert(
        !pass_result.required_downward_update_ops.empty(),
        "local SPR scoring produced zero required downward update ops");

    if (already_updated_ops == nullptr) {
        UpdateTreeClvsAfterPrune(
            workspace.scoring_dev,
            pruned_tree,
            pruned_host,
            workspace.tree_ops,
            upward_start_node,
            pass_result.required_downward_update_ops,
            ctx.stream);
    } else {
        std::vector<NodeOpInfo> new_update_ops =
            filter_local_spr_new_ops(
                pass_result.required_downward_update_ops,
                *already_updated_ops);
        if (!new_update_ops.empty()) {
            UpdateTreeClvsAfterPrune(
                workspace.scoring_dev,
                pruned_tree,
                pruned_host,
                workspace.tree_ops,
                -1,
                new_update_ops,
                ctx.stream);
        }
    }

    copy_unscaled_up_clv_to_query_slot(
        ctx.res.dev,
        prune_root_id,
        workspace.scoring_dev,
        0,
        ctx.stream);

    UploadPlacementOps(
        workspace.candidate_ops,
        candidate_ops,
        ctx.stream);

    DeviceTree query_view =
        make_query_view(workspace.scoring_dev, 0);
    pass_result.placement_result = PlacementEvaluationKernel(
        query_view,
        workspace.candidate_ops.d_ops,
        workspace.candidate_ops.num_ops,
        1,
        ctx.stream,
        false);
    local_spr_assert(
        pass_result.placement_result.top_placements.size() ==
            candidate_ops.size(),
        "local SPR scoring did not return a score for every edge");
    return pass_result;
}

double find_local_spr_baseline_loglikelihood(
    const PlacementResult& placement_result,
    int baseline_target_id)
{
    for (const PlacementResult::RankedPlacement& placement :
         placement_result.top_placements) {
        if (placement.target_id == baseline_target_id) {
            return placement.loglikelihood;
        }
    }
    return -std::numeric_limits<double>::infinity();
}

// Convert placement scores into concrete SPR proposals that can later be validated
// against the full tree likelihood.

void append_positive_gain_local_spr_candidates(
    const PlacementResult& placement_result,
    const TreeBuildResult& base_tree,
    const TreeBuildResult& pruned_tree,
    const LocalSprRepairUnit& unit,
    int prune_root_id,
    const std::vector<int>& subtree_nodes,
    double baseline_logL,
    int candidate_pool_limit,
    std::vector<LocalSprCandidateMove>& unit_topk)
{
    for (const PlacementResult::RankedPlacement& placement :
         placement_result.top_placements) {
        if (placement.target_id < 0 ||
            placement.target_id >= static_cast<int>(pruned_tree.nodes.size())) {
            continue;
        }
        const int edge_child = placement.target_id;
        const int edge_parent =
            pruned_tree.nodes[(size_t)edge_child].parent;
        if (edge_parent < 0) continue;

        const double approx_gain =
            placement.loglikelihood - baseline_logL;
        if (!(approx_gain > 0.0)) {
            continue;
        }

        LocalSprCandidateMove candidate;
        candidate.repair_unit_id = unit.unit_id;
        candidate.prune_root_id = prune_root_id;
        candidate.regraft_child_id = edge_child;
        candidate.regraft_parent_id = edge_parent;
        candidate.old_parent_id =
            base_tree.nodes[(size_t)prune_root_id].parent;
        candidate.approx_gain = approx_gain;
        candidate.pendant_length = placement.pendant_length;
        candidate.proximal_length =
            static_cast<double>(
                pruned_tree.nodes[(size_t)edge_child].branch_length_to_parent) -
            placement.proximal_length;
        candidate.subtree_nodes = subtree_nodes;
        const int path_start =
            candidate.old_parent_id >= 0
                ? candidate.old_parent_id
                : prune_root_id;
        candidate.regraft_path_nodes =
            build_tree_path_nodes(base_tree, path_start, edge_child);
        keep_local_spr_topk(
            unit_topk,
            std::move(candidate),
            candidate_pool_limit);
    }
}

// Score every plausible regraft edge for one chosen prune root, optionally expand
// outward from the best inner-ring edges, and accumulate the top local candidates.

void evaluate_local_spr_prune_root(
    const LocalSprSearchContext& ctx,
    const TreeBuildResult& base_tree,
    const LocalSprRepairUnit& unit,
    const LocalSprPruneRootWorkItem& work_item,
    int inner_search_radius,
    bool staged_radius_expansion,
    int outer_seed_limit,
    int candidate_pool_limit,
    LocalSprScoringWorkspace& workspace,
    std::vector<LocalSprCandidateMove>& unit_topk,
    LocalSprSearchSummary& search_summary)
{
    const int prune_root_id = work_item.prune_root_id;
    const std::vector<int>& subtree_nodes = work_item.subtree_nodes;
    const std::vector<char> subtree_mask =
        build_local_spr_subtree_mask(base_tree, subtree_nodes);

    TreeBuildResult pruned_tree = base_tree;
    PruneInfo prune_info;
    if (!prune_subtree_for_spr(pruned_tree, prune_root_id, prune_info)) {
        return;
    }

    const std::vector<int>& legal_candidate_edges =
        work_item.legal_inner_candidate_edges;
    search_summary.enumerated_candidates += legal_candidate_edges.size();
    if (legal_candidate_edges.empty()) {
        return;
    }

    HostPacking pruned_host = ctx.res.hostPack;
    rebuild_host_topology_from_tree_local(pruned_tree, pruned_host);
    pruned_host.pattern_weights = ctx.pattern_weights_arg;
    const std::vector<int> node_to_tip =
        build_local_spr_node_to_tip(pruned_tree, pruned_host);

    int changed_nodes[1] = { prune_info.sibling_id };
    fill_pmats_in_host_packing(
        pruned_tree,
        pruned_host,
        ctx.res.eig,
        ctx.rate_multipliers,
        ctx.states,
        ctx.rate_cats,
        changed_nodes,
        1);
    reload_device_tree_live_data_local_spr(
        workspace.scoring_dev,
        pruned_tree,
        pruned_host,
        ctx.res.hostPack,
        prune_info.sibling_id,
        workspace.previous_main_pmat_node,
        nullptr,
        ctx.stream,
        nullptr);

    // Stage 1: score the inner regraft edges around the prune site.
    copy_upward_state(
        ctx.res.dev,
        workspace.scoring_dev,
        ctx.stream);
    LocalSprPlacementPassResult inner_pass =
        run_local_spr_placement_pass(
            ctx,
            pruned_tree,
            pruned_host,
            node_to_tip,
            legal_candidate_edges,
            prune_info.grandparent_id,
            prune_root_id,
            workspace,
            nullptr);
    PlacementResult placement_result = std::move(inner_pass.placement_result);

    const double baseline_logL =
        (prune_info.grandparent_id >= 0)
            ? find_local_spr_baseline_loglikelihood(
                  placement_result,
                  prune_info.sibling_id)
            : -std::numeric_limits<double>::infinity();
    if (!std::isfinite(baseline_logL)) {
        return;
    }

    // Stage 2: if the inner ring looks promising, expand outward once.
    if (staged_radius_expansion) {
        const std::vector<int> center_dist = compute_center_distances(
            pruned_tree,
            prune_info.sibling_id,
            prune_info.grandparent_id);
        const std::vector<int> outer_seed_edges =
            select_local_spr_seed_edges(
                placement_result,
                pruned_tree,
                center_dist,
                baseline_logL,
                outer_seed_limit,
                inner_search_radius);
        if (!outer_seed_edges.empty()) {
            std::vector<int> outer_candidate_edges;
            outer_candidate_edges.reserve(
                outer_seed_edges.size() *
                std::max(1, ctx.local_spr_radius - inner_search_radius));
            std::unordered_set<int> seen_outer_edges;
            for (int seed_edge_child_id : outer_seed_edges) {
                const std::vector<int> expanded_edges =
                    collect_outward_edges_from_seed(
                        pruned_tree,
                        center_dist,
                        seed_edge_child_id,
                        inner_search_radius,
                        ctx.local_spr_radius,
                        prune_info.pruned_id,
                        prune_info.free_internal_id);
                for (int edge_child_id : expanded_edges) {
                    if (!seen_outer_edges.insert(edge_child_id).second) {
                        continue;
                    }
                    outer_candidate_edges.push_back(edge_child_id);
                }
            }

            const std::vector<int> outer_legal_candidate_edges =
                filter_local_spr_candidate_edges(
                    pruned_tree,
                    unit.envelope_mask,
                    subtree_mask,
                    outer_candidate_edges);
            search_summary.enumerated_candidates +=
                outer_legal_candidate_edges.size();
            if (!outer_legal_candidate_edges.empty()) {
                LocalSprPlacementPassResult outer_pass =
                    run_local_spr_placement_pass(
                        ctx,
                        pruned_tree,
                        pruned_host,
                        node_to_tip,
                        outer_legal_candidate_edges,
                        -1,
                        prune_root_id,
                        workspace,
                        &inner_pass.required_downward_update_ops);
                placement_result.top_placements.insert(
                    placement_result.top_placements.end(),
                    std::make_move_iterator(
                        outer_pass.placement_result.top_placements.begin()),
                    std::make_move_iterator(
                        outer_pass.placement_result.top_placements.end()));
                std::sort(
                    placement_result.top_placements.begin(),
                    placement_result.top_placements.end(),
                    [](const PlacementResult::RankedPlacement& lhs,
                       const PlacementResult::RankedPlacement& rhs) {
                        return lhs.loglikelihood > rhs.loglikelihood;
                    });
            }
        }
    }

    // Stage 3: convert positive-gain placements into commit candidates.
    append_positive_gain_local_spr_candidates(
        placement_result,
        base_tree,
        pruned_tree,
        unit,
        prune_root_id,
        subtree_nodes,
        baseline_logL,
        candidate_pool_limit,
        unit_topk);
}

// Search phase: enumerate candidate prune roots inside each repair unit and keep
// only the strongest local SPR moves for later validation.

std::vector<LocalSprCandidateMove> rank_local_spr_candidates(
    const LocalSprSearchContext& ctx,
    const TreeBuildResult& base_tree,
    const std::vector<LocalSprRepairUnit>& repair_units,
    LocalSprSearchSummary& search_summary)
{
    std::vector<LocalSprCandidateMove> ranked_candidates;
    PlacementQueryBatch subtree_query_batch;
    subtree_query_batch.count = 1;
    subtree_query_batch.branch_lengths.assign(1, fp_t(0.5));
    subtree_query_batch.query_chars.assign(
        ctx.sites,
        static_cast<uint8_t>(ctx.states == 4 ? 15 : 4));
    LocalSprScoringWorkspace scoring_workspace;

    try {
        scoring_workspace.scoring_dev = upload_to_gpu(
            base_tree,
            ctx.res.hostPack,
            ctx.res.eig,
            ctx.rate_weights,
            ctx.rate_multipliers,
            ctx.pi,
            ctx.sites,
            ctx.states,
            ctx.rate_cats,
            ctx.per_rate_scaling,
            &subtree_query_batch,
            false);

        const int inner_search_radius = std::min(ctx.local_spr_radius, 2);
        const bool staged_radius_expansion =
            ctx.local_spr_radius > inner_search_radius;
        const int local_spr_candidate_pool_limit = ctx.local_spr_topk_per_unit;
        const int local_spr_outer_seed_limit = ctx.local_spr_topk_per_unit;
        const bool local_spr_use_tier_early_stop = true;
        ranked_candidates.reserve(
            repair_units.size() * static_cast<size_t>(
                std::max(1, local_spr_candidate_pool_limit)));
        for (const LocalSprRepairUnit& unit : repair_units) {
            std::vector<LocalSprCandidateMove> unit_topk;
            const std::vector<char> anchor_skeleton_mask =
                build_unit_anchor_skeleton_mask(base_tree, unit.anchor_ids);
            const std::vector<int> skeleton_distance_by_node =
                multi_source_bfs_distances(base_tree, anchor_skeleton_mask);
            std::vector<LocalSprPruneRootWorkItem> prune_root_work_items =
                build_local_spr_prune_root_work_items(
                    base_tree,
                    unit,
                    inner_search_radius,
                    skeleton_distance_by_node);
            if (prune_root_work_items.empty()) {
                continue;
            }

            bool stop_after_next_tier = false;
            for (size_t work_idx = 0; work_idx < prune_root_work_items.size();) {
                const int tier_distance =
                    prune_root_work_items[work_idx].skeleton_distance;
                size_t tier_end = work_idx;
                while (tier_end < prune_root_work_items.size() &&
                       prune_root_work_items[tier_end].skeleton_distance ==
                           tier_distance) {
                    ++tier_end;
                }

                const size_t prev_topk_size = unit_topk.size();
                for (; work_idx < tier_end; ++work_idx) {
                    evaluate_local_spr_prune_root(
                        ctx,
                        base_tree,
                        unit,
                        prune_root_work_items[work_idx],
                        inner_search_radius,
                        staged_radius_expansion,
                        local_spr_outer_seed_limit,
                        local_spr_candidate_pool_limit,
                        scoring_workspace,
                        unit_topk,
                        search_summary);
                }

                const bool tier_found_positive =
                    unit_topk.size() > prev_topk_size;
                if (local_spr_use_tier_early_stop && stop_after_next_tier) {
                    break;
                }
                if (local_spr_use_tier_early_stop && tier_found_positive) {
                    stop_after_next_tier = true;
                }
            }
            search_summary.retained_candidates += unit_topk.size();
            ranked_candidates.insert(
                ranked_candidates.end(),
                std::make_move_iterator(unit_topk.begin()),
                std::make_move_iterator(unit_topk.end()));
        }
    } catch (...) {
        release_local_spr_scoring_workspace(scoring_workspace, ctx.stream);
        throw;
    }
    release_local_spr_scoring_workspace(scoring_workspace, ctx.stream);
    return ranked_candidates;
}

// Validation phase: replay the proposed topology edit on a fresh tree copy and
// accept it only if the full-tree likelihood really improves.

int validate_local_spr_candidates(
    std::vector<LocalSprCandidateMove> validation_candidates,
    const std::vector<LocalSprRepairUnit>& repair_units,
    bool dynamic_validation_conflicts,
    int local_spr_radius,
    const std::function<double(const TreeBuildResult&)>& eval_logL,
    TreeBuildResult& base_tree,
    double& current_logL)
{
    const double logL_accept_eps = 0.0;
    int accepted_candidates = 0;
    for (size_t candidate_idx = 0;
         candidate_idx < validation_candidates.size();
         ++candidate_idx) {
        const LocalSprCandidateMove candidate =
            validation_candidates[candidate_idx];
        if (candidate.repair_unit_id < 0 ||
            candidate.repair_unit_id >= static_cast<int>(repair_units.size())) {
            continue;
        }

        const LocalSprRepairUnit& unit =
            repair_units[(size_t)candidate.repair_unit_id];
        if (dynamic_validation_conflicts &&
            !local_spr_candidate_still_legal(
                base_tree,
                unit,
                candidate,
                nullptr)) {
            continue;
        }

        TreeBuildResult candidate_tree = base_tree;
        PruneInfo prune_info;
        if (!prune_subtree_for_spr(candidate_tree, candidate.prune_root_id, prune_info)) {
            continue;
        }

        std::vector<int> current_subtree_nodes;
        collect_subtree_node_ids(candidate_tree, prune_info.pruned_id, current_subtree_nodes);
        const std::vector<char> current_subtree_mask =
            build_local_spr_subtree_mask(candidate_tree, current_subtree_nodes);
        if (candidate.regraft_child_id < 0 ||
            candidate.regraft_child_id >= static_cast<int>(candidate_tree.nodes.size())) {
            continue;
        }
        const int regraft_parent =
            candidate_tree.nodes[(size_t)candidate.regraft_child_id].parent;
        if (regraft_parent < 0) continue;
        if (current_subtree_mask[(size_t)candidate.regraft_child_id] ||
            current_subtree_mask[(size_t)regraft_parent]) {
            continue;
        }
        if (!unit.envelope_mask[(size_t)candidate.regraft_child_id] ||
            !unit.envelope_mask[(size_t)regraft_parent]) {
            continue;
        }

        if (!dynamic_validation_conflicts) {
            std::vector<int> legal_edges = collect_candidate_edges(
                candidate_tree,
                prune_info.sibling_id,
                prune_info.grandparent_id,
                local_spr_radius,
                prune_info.pruned_id,
                prune_info.free_internal_id);
            local_spr_assert_candidate_legal(
                base_tree,
                unit,
                candidate.prune_root_id,
                candidate.subtree_nodes,
                legal_edges,
                candidate.regraft_child_id);
        }

        regraft_subtree_for_spr(
            candidate_tree,
            prune_info,
            candidate.regraft_child_id,
            true,
            candidate.pendant_length,
            candidate.proximal_length);
        local_spr_assert(
            subtree_fully_inside_mask(
                candidate_tree,
                prune_info.pruned_id,
                unit.envelope_mask),
            "committed subtree escaped repair envelope after regraft");
        const double candidate_logL = eval_logL(candidate_tree);
        if (candidate_logL > current_logL + logL_accept_eps) {
            current_logL = candidate_logL;
            base_tree = std::move(candidate_tree);
            ++accepted_candidates;
            if (dynamic_validation_conflicts &&
                candidate_idx + 1 < validation_candidates.size()) {
                auto keep_begin =
                    validation_candidates.begin() +
                    static_cast<std::ptrdiff_t>(candidate_idx + 1);
                keep_begin = std::remove_if(
                    keep_begin,
                    validation_candidates.end(),
                    [&](const LocalSprCandidateMove& pending_candidate) {
                        if (pending_candidate.repair_unit_id < 0 ||
                            pending_candidate.repair_unit_id >=
                                static_cast<int>(repair_units.size())) {
                            return true;
                        }
                        return !local_spr_candidate_still_legal(
                            base_tree,
                            repair_units[(size_t)pending_candidate.repair_unit_id],
                            pending_candidate,
                            nullptr);
                    });
                validation_candidates.erase(
                    keep_begin,
                    validation_candidates.end());
            }
        }
    }
    return accepted_candidates;
}

// Once a validation round accepts one or more moves, rebuild the authoritative
// host/device tree state so the next round starts from a consistent baseline.

void rebuild_after_local_spr_round(
    LocalSprBatchRunContext& ctx,
    const TreeBuildResult& base_tree)
{
    free_placement_op_buffer(ctx.placement_ops, ctx.stream);
    cudaStreamSynchronize(ctx.stream);

    ctx.placement_ops = PlacementOpBuffer{};
    ctx.placement_ops.profile_commit_timing = ctx.profile_batch_timing;
    const PlacementQueryBatch empty_queries;
    HostPacking rebuilt_host_pack =
        build_local_spr_tree_host_packing(ctx, base_tree);
    ctx.res.tree = base_tree;
    ctx.res.hostPack = std::move(rebuilt_host_pack);
    ctx.res.queries = PlacementQueryBatch{};
    reload_device_tree_live_data(
        ctx.res.dev,
        ctx.res.tree,
        ctx.res.hostPack,
        &empty_queries,
        ctx.stream);
    UpdateTreeClvs(
        ctx.res.dev,
        ctx.res.tree,
        ctx.res.hostPack,
        ctx.placement_ops,
        ctx.stream);
}

} // namespace

// Public entrypoint: run several rounds of local SPR search, validation, and rebuild
// on top of the current committed tree for this batch of inserted queries.

void run_local_spr_batch_refinement(LocalSprBatchRunContext& ctx) {
    if (ctx.inserted_names.empty()) {
        std::cout << "Local SPR skipped: no inserted placements in this batch.\n";
        return;
    }

    LocalSprEnvGuard env_guard(ctx.local_spr_fast, ctx.res.tree.nodes.size());
    int local_spr_rounds_executed = 0;
    LocalSprEvalWorkspace eval_workspace;

    const LocalSprSearchContext search_ctx{
        ctx.res,
        ctx.stream,
        ctx.pattern_weights_arg,
        ctx.rate_weights,
        ctx.rate_multipliers,
        ctx.pi,
        ctx.sites,
        ctx.states,
        ctx.rate_cats,
        ctx.per_rate_scaling,
        ctx.local_spr_radius,
        ctx.local_spr_topk_per_unit,
    };

    try {
        for (int local_spr_round = 1;
             local_spr_round <= ctx.local_spr_rounds;
             ++local_spr_round) {
            ++local_spr_rounds_executed;
            double current_logL = root_likelihood::compute_root_loglikelihood_total(
                ctx.res.dev,
                ctx.res.tree.root_id,
                ctx.res.dev.d_pattern_weights_u,
                nullptr,
                0.0,
                0);
            TreeBuildResult base_tree = ctx.res.tree;

            const std::vector<LocalSprRepairUnit> repair_units =
                prepare_local_spr_repair_units(
                    base_tree,
                    ctx.inserted_names,
                    ctx.local_spr_cluster_threshold,
                    ctx.local_spr_radius);

            LocalSprSearchSummary search_summary;
            search_summary.unit_count = repair_units.size();
            std::vector<LocalSprCandidateMove> ranked_candidates =
                rank_local_spr_candidates(
                    search_ctx,
                    base_tree,
                    repair_units,
                    search_summary);

            std::sort(
                ranked_candidates.begin(),
                ranked_candidates.end(),
                [](const LocalSprCandidateMove& lhs, const LocalSprCandidateMove& rhs) {
                    return lhs.approx_gain > rhs.approx_gain;
                });

            std::vector<LocalSprCandidateMove> validation_candidates;
            if (ctx.local_spr_dynamic_validation_conflicts) {
                validation_candidates = ranked_candidates;
            } else {
                validation_candidates = select_local_spr_candidates(
                    ranked_candidates,
                    static_cast<int>(base_tree.nodes.size()));
            }
            search_summary.selected_candidates = validation_candidates.size();

            std::cout << "Local SPR round " << local_spr_round
                      << "/" << ctx.local_spr_rounds
                      << " subtree repair units: " << search_summary.unit_count
                      << ", enumerated candidates: " << search_summary.enumerated_candidates
                      << ", retained candidates: " << search_summary.retained_candidates
                      << ", selected: " << search_summary.selected_candidates
                      << " (radius=" << ctx.local_spr_radius
                      << ", cluster=" << ctx.local_spr_cluster_threshold
                      << ", topk=" << ctx.local_spr_topk_per_unit
                      << ", selection="
                      << (ctx.local_spr_dynamic_validation_conflicts
                              ? "dynamic-validation"
                              : "one-per-unit")
                      << ", pool=per-unit-topk"
                      << ")\n";

            const int accepted_this_round = validate_local_spr_candidates(
                std::move(validation_candidates),
                repair_units,
                ctx.local_spr_dynamic_validation_conflicts,
                ctx.local_spr_radius,
                [&](const TreeBuildResult& candidate_tree) -> double {
                    return evaluate_local_spr_tree_loglikelihood(
                        ctx,
                        candidate_tree,
                        eval_workspace);
                },
                base_tree,
                current_logL);
            if (accepted_this_round <= 0) {
                std::cout << "Local SPR round " << local_spr_round
                          << " accepted no candidates; stopping.\n";
                break;
            }
            std::cout << "Local SPR round " << local_spr_round
                      << " accepted " << accepted_this_round
                      << " candidate(s).\n";
            ctx.current_tree_newick = mltreeio::write_tree_to_newick_string(base_tree);
            rebuild_after_local_spr_round(
                ctx,
                base_tree);
        }
    } catch (...) {
        release_local_spr_eval_workspace(eval_workspace, ctx.stream);
        throw;
    }

    release_local_spr_eval_workspace(eval_workspace, ctx.stream);
    std::cout << "Local SPR joint refinement done (radius=" << ctx.local_spr_radius
              << ", rounds=" << local_spr_rounds_executed
              << "/" << ctx.local_spr_rounds << ").\n";
}
