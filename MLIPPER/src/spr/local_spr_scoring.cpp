#include "local_spr_internal.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

#include <tbb/parallel_for.h>

#include "placement/placement.cuh"
#include "nni.hpp"
#include "util/checked_size.hpp"

namespace {

bool local_spr_edge_touches_forbidden_node(
    const TopologyRefinementSearchContext& ctx,
    const TreeBuildResult& tree,
    int edge_child)
{
    if (edge_child < 0 ||
        edge_child >= static_cast<int>(tree.nodes.size())) {
        return true;
    }
    const int edge_parent = tree.nodes[static_cast<size_t>(edge_child)].parent;
    if (edge_parent < 0) {
        return true;
    }
    const auto& forbidden = ctx.forbidden_regraft_node_ids;
    return std::find(forbidden.begin(), forbidden.end(), edge_child) !=
               forbidden.end() ||
           std::find(forbidden.begin(), forbidden.end(), edge_parent) !=
               forbidden.end();
}

// Placement targets are identified by the child endpoint of an edge. After a
// root-adjacent prune, the sibling is promoted to the artificial root and no
// longer identifies an edge. Use its left-child edge as the stable baseline
// representation in that one case.
int canonical_unchanged_topology_edge_child(
    const TreeBuildResult& pruned_tree,
    int sibling_id)
{
    if (sibling_id < 0 ||
        sibling_id >= static_cast<int>(pruned_tree.nodes.size())) {
        return -1;
    }

    const TreeNode& sibling =
        pruned_tree.nodes[static_cast<size_t>(sibling_id)];
    if (sibling.parent >= 0) {
        return sibling_id;
    }
    return sibling.is_tip ? -1 : sibling.left;
}

} // namespace

PlacementQueryBatch build_local_spr_scoring_query_batch(
    const TopologyRefinementSearchContext& ctx,
    int query_capacity)
{
    PlacementQueryBatch query_batch;
    query_batch.count = static_cast<size_t>(std::max(1, query_capacity));
    const size_t query_char_count = mlipper::util::checked_product(
        "local SPR scoring query characters",
        query_batch.count,
        ctx.sites);
    query_batch.query_chars.assign(
        query_char_count,
        static_cast<uint8_t>(ctx.states == 4 ? 15 : 4));
    return query_batch;
}

std::vector<LocalSPRPruneRootWorkItem> build_local_spr_prune_root_work_items(
    const TreeBuildResult& base_tree,
    const LocalSPRRepairUnit& unit,
    int inner_search_radius,
    const std::vector<int>& skeleton_distance_by_node,
    mlipper::TopologyMoveType move_type,
    const std::unordered_set<int>& allowed_nni_central_edge_child_ids)
{
    std::vector<LocalSPRPruneRootWorkItem> work_items;
    work_items.reserve(unit.envelope_nodes.size());

    for (int candidate_prune_root_id : unit.envelope_nodes) {
        if (candidate_prune_root_id < 0 ||
            candidate_prune_root_id >= static_cast<int>(base_tree.nodes.size()) ||
            candidate_prune_root_id == base_tree.root_id ||
            base_tree.nodes[(size_t)candidate_prune_root_id].parent < 0) {
            continue;
        }
        if (move_type == mlipper::TopologyMoveType::NNI) {
            const int central_child_id =
                base_tree.nodes[(size_t)candidate_prune_root_id].parent;
            if (allowed_nni_central_edge_child_ids.count(central_child_id) == 0) {
                continue;
            }
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

        std::vector<int> inner_candidate_edges;
        if (move_type == mlipper::TopologyMoveType::NNI) {
            inner_candidate_edges = mlipper::nni::candidate_regraft_edges(
                pruned_tree, prune_info);
        } else {
            inner_candidate_edges = collect_candidate_edges(
                pruned_tree,
                prune_info.sibling_id,
                prune_info.grandparent_id,
                inner_search_radius,
                prune_info.pruned_id,
                prune_info.free_internal_id);
        }
        std::vector<int> legal_inner_candidate_edges =
            filter_local_spr_candidate_edges(
                pruned_tree,
                unit.envelope_mask,
                subtree_mask,
                inner_candidate_edges);
        // Keep the unchanged-topology score even when the envelope excludes
        // one endpoint. Every alternate move needs this zero-gain reference.
        const int baseline_edge_child_id =
            canonical_unchanged_topology_edge_child(
                pruned_tree,
                prune_info.sibling_id);
        if (baseline_edge_child_id >= 0 &&
            std::find(
                legal_inner_candidate_edges.begin(),
                legal_inner_candidate_edges.end(),
                baseline_edge_child_id) == legal_inner_candidate_edges.end()) {
            legal_inner_candidate_edges.push_back(baseline_edge_child_id);
        }
        if (legal_inner_candidate_edges.empty()) {
            continue;
        }

        LocalSPRPruneRootWorkItem work_item;
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
        [](const LocalSPRPruneRootWorkItem& lhs,
           const LocalSPRPruneRootWorkItem& rhs) {
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

std::vector<LocalSPRPruneRootWorkItem> prepare_local_spr_unit_work_items(
    const TreeBuildResult& base_tree,
    const LocalSPRRepairUnit& unit,
    int inner_search_radius,
    const TopologyRefinementSearchContext& ctx)
{
    std::vector<char> anchor_skeleton_mask(base_tree.nodes.size(), 0);
    if (ctx.move_type == mlipper::TopologyMoveType::NNI) {
        for (int central_child_id : ctx.allowed_nni_central_edge_child_ids) {
            if (central_child_id < 0 ||
                central_child_id >= static_cast<int>(base_tree.nodes.size())) {
                continue;
            }
            anchor_skeleton_mask[static_cast<size_t>(central_child_id)] = 1;
            const int parent_id = base_tree.nodes[
                static_cast<size_t>(central_child_id)].parent;
            if (parent_id >= 0) {
                anchor_skeleton_mask[static_cast<size_t>(parent_id)] = 1;
            }
        }
    } else {
        anchor_skeleton_mask =
            build_unit_anchor_skeleton_mask(base_tree, unit.anchors);
    }
    const std::vector<int> skeleton_distance_by_node =
        mlipper::divide_and_conquer::multi_source_bfs_distances(
            base_tree,
            anchor_skeleton_mask);
    return build_local_spr_prune_root_work_items(
        base_tree,
        unit,
        inner_search_radius,
        skeleton_distance_by_node,
        ctx.move_type,
        ctx.allowed_nni_central_edge_child_ids);
}

size_t next_local_spr_scratch_capacity(
    size_t current_capacity,
    size_t required_ops)
{
    size_t capacity = std::max(required_ops, static_cast<size_t>(256));
    if (current_capacity > capacity) {
        capacity = current_capacity;
    }
    while (capacity < required_ops) {
        capacity *= 2;
    }
    return capacity;
}

template <typename T>
void cuda_free_if_set(T*& ptr, cudaStream_t stream)
{
    if (ptr == nullptr) {
        return;
    }
    CUDA_CHECK(cudaStreamSynchronize(stream));
    CUDA_CHECK(cudaFree(ptr));
    ptr = nullptr;
}

template <typename T>
void cuda_malloc_bytes(T*& ptr, size_t bytes)
{
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&ptr), bytes));
}

template <typename T>
void cuda_copy_device_to_device_async(
    T* dst,
    const T* src,
    size_t elems,
    cudaStream_t stream)
{
    CUDA_CHECK(cudaMemcpyAsync(
        dst,
        src,
        mlipper::util::checked_allocation_bytes<T>(
            elems, "local SPR device-to-device copy"),
        cudaMemcpyDeviceToDevice,
        stream));
}

void release_local_spr_placement_scratch(
    PlacementScratchOverride& scratch,
    cudaStream_t stream)
{
    cuda_free_if_set(scratch.d_sumtable, stream);
    cuda_free_if_set(scratch.d_likelihoods, stream);
    cuda_free_if_set(scratch.d_query_pmat, stream);
    cuda_free_if_set(scratch.d_prev_loglk, stream);
    cuda_free_if_set(scratch.d_active_ops, stream);
    cuda_free_if_set(scratch.d_pmat_mid_prox, stream);
    cuda_free_if_set(scratch.d_pmat_mid_dist, stream);
    cuda_free_if_set(scratch.d_edge_midpoint_clv, stream);
    cuda_free_if_set(scratch.d_edge_midpoint_scaler, stream);
    scratch = PlacementScratchOverride{};
}

void release_local_spr_embedded_placement_scratch(
    DeviceTree& dev,
    cudaStream_t stream)
{
    if (dev.d_query_pmat != nullptr ||
        dev.d_pmat_mid_prox != nullptr ||
        dev.d_pmat_mid_dist != nullptr) {
        ensure_device_tree_current_device(
            dev,
            "release_local_spr_embedded_placement_scratch");
    }
    cuda_free_if_set(dev.d_query_pmat, stream);
    cuda_free_if_set(dev.d_pmat_mid_prox, stream);
    cuda_free_if_set(dev.d_pmat_mid_dist, stream);
}

void ensure_local_spr_placement_scratch_capacity(
    const DeviceTree& dev,
    size_t required_ops,
    PlacementScratchOverride& scratch,
    cudaStream_t stream)
{
    if (required_ops == 0) {
        return;
    }
    if (scratch.d_sumtable != nullptr &&
        scratch.d_likelihoods != nullptr &&
        scratch.d_query_pmat != nullptr &&
        scratch.d_prev_loglk != nullptr &&
        scratch.d_active_ops != nullptr &&
        scratch.d_pmat_mid_prox != nullptr &&
        scratch.d_pmat_mid_dist != nullptr &&
        scratch.sumtable_capacity_ops >= required_ops &&
        scratch.likelihood_capacity_ops >= required_ops &&
        scratch.query_pmat_capacity_ops >= required_ops &&
        scratch.ranking_state_capacity_ops >= required_ops &&
        scratch.midpoint_pmat_capacity_nodes >=
            static_cast<size_t>(std::max(dev.N, 0))) {
        return;
    }

    ensure_device_tree_current_device(
        dev,
        "ensure_local_spr_placement_scratch_capacity");

    const size_t sumtable_stride =
        dev.sites *
        static_cast<size_t>(dev.rate_cats) *
        static_cast<size_t>(dev.states);
    const size_t query_pmat_stride = dev.pmat_per_node_elems();
    const size_t midpoint_pmat_stride = dev.pmat_per_node_elems();
    const size_t required_nodes = static_cast<size_t>(std::max(dev.N, 0));
    if (sumtable_stride == 0) {
        throw std::runtime_error(
            "local SPR cannot allocate placement scratch for an empty device shape.");
    }
    if (query_pmat_stride == 0) {
        throw std::runtime_error(
            "local SPR cannot allocate query PMAT scratch for an empty device shape.");
    }
    if (midpoint_pmat_stride == 0 || required_nodes == 0) {
        throw std::runtime_error(
            "local SPR cannot allocate midpoint PMAT scratch for an empty device shape.");
    }

    const size_t target_capacity = next_local_spr_scratch_capacity(
        std::max(
            std::max(
                scratch.sumtable_capacity_ops,
                scratch.likelihood_capacity_ops),
            std::max(
                scratch.query_pmat_capacity_ops,
                scratch.ranking_state_capacity_ops)),
        required_ops);
    size_t target_node_capacity = std::max(
        required_nodes,
        scratch.midpoint_pmat_capacity_nodes);
    if (target_node_capacity == 0) {
        target_node_capacity = required_nodes;
    }

    release_local_spr_placement_scratch(scratch, stream);
    cuda_malloc_bytes(scratch.d_sumtable, mlipper::util::checked_product(
        "local SPR sumtable", sizeof(fp_t), sumtable_stride, target_capacity));
    cuda_malloc_bytes(scratch.d_likelihoods,
        mlipper::util::checked_allocation_bytes<fp_t>(target_capacity,
            "local SPR likelihoods"));
    cuda_malloc_bytes(scratch.d_query_pmat, mlipper::util::checked_product(
        "local SPR query matrices", sizeof(fp_t), query_pmat_stride, target_capacity));
    cuda_malloc_bytes(scratch.d_prev_loglk,
        mlipper::util::checked_allocation_bytes<fp_t>(target_capacity,
            "local SPR previous likelihoods"));
    cuda_malloc_bytes(scratch.d_active_ops,
        mlipper::util::checked_allocation_bytes<int>(target_capacity,
            "local SPR active operations"));
    cuda_malloc_bytes(scratch.d_pmat_mid_prox, mlipper::util::checked_product(
        "local SPR proximal midpoint matrices", sizeof(fp_t), midpoint_pmat_stride,
        target_node_capacity));
    cuda_malloc_bytes(scratch.d_pmat_mid_dist, mlipper::util::checked_product(
        "local SPR distal midpoint matrices", sizeof(fp_t), midpoint_pmat_stride,
        target_node_capacity));
    scratch.sumtable_capacity_ops = target_capacity;
    scratch.likelihood_capacity_ops = target_capacity;
    scratch.query_pmat_capacity_ops = target_capacity;
    scratch.ranking_state_capacity_ops = target_capacity;
    scratch.midpoint_pmat_capacity_nodes = target_node_capacity;
}

void release_local_spr_scoring_workspace(
    LocalSPRScoringWorkspace& workspace,
    cudaStream_t stream)
{
    cuda_free_if_set(workspace.d_upward_restore_nodes, stream);
    free_placement_op_buffer(workspace.candidate_ops, stream);
    free_placement_op_buffer(workspace.tree_ops, stream);
    release_local_spr_placement_scratch(
        workspace.placement_scratch,
        stream);
    CUDA_CHECK(cudaStreamSynchronize(stream));
    release_subtree_workspace(workspace.subtree_workspace);
    workspace = LocalSPRScoringWorkspace{};
}

void restore_local_spr_upward_state(
    const DeviceTree& source,
    LocalSPRScoringWorkspace& workspace,
    cudaStream_t stream)
{
    if (!workspace.upward_state_initialized) {
        copy_upward_state(source, workspace.subtree_workspace.dev, stream);
        workspace.upward_state_initialized = true;
        workspace.dirty_upward_nodes.clear();
        return;
    }
    auto& nodes = workspace.dirty_upward_nodes;
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    if (nodes.empty()) return;
    if (nodes.size() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error(
            "Local SPR dirty-node list exceeds kernel indexing.");
    }
    if (workspace.upward_restore_capacity < static_cast<int>(nodes.size())) {
        cuda_free_if_set(workspace.d_upward_restore_nodes, stream);
        workspace.upward_restore_capacity = static_cast<int>(nodes.size());
        cuda_malloc_bytes(
            workspace.d_upward_restore_nodes,
            mlipper::util::checked_allocation_bytes<int>(
                static_cast<size_t>(workspace.upward_restore_capacity),
                "local SPR upward restore nodes"));
    }
    CUDA_CHECK(cudaMemcpyAsync(
        workspace.d_upward_restore_nodes,
        nodes.data(),
        mlipper::util::checked_allocation_bytes<int>(
            nodes.size(), "local SPR upward restore node upload"),
        cudaMemcpyHostToDevice,
        stream));
    copy_selected_upward_state(
        source,
        workspace.subtree_workspace.dev,
        workspace.d_upward_restore_nodes,
        static_cast<int>(nodes.size()),
        stream);
    nodes.clear();
}

void remember_local_spr_dirty_upward_nodes(LocalSPRScoringWorkspace& workspace)
{
    for (const NodeOpInfo& op : workspace.tree_ops.upward_ops_host) {
        if (op.parent_id >= 0) {
            workspace.dirty_upward_nodes.push_back(op.parent_id);
        }
    }
}

size_t estimate_local_spr_scoring_lane_bytes(const DeviceTree& dev)
{
    const size_t nodes = mlipper::util::checked_add_size(
        static_cast<size_t>(std::max(dev.N, 0)), 2,
        "local SPR lane node capacity");
    const size_t tips = mlipper::util::checked_add_size(
        static_cast<size_t>(std::max(dev.tips, 0)), 1,
        "local SPR lane tip capacity");
    const size_t sites = dev.sites;
    const size_t rates = static_cast<size_t>(std::max(dev.rate_cats, 1));
    const size_t states = static_cast<size_t>(std::max(dev.states, 1));
    const size_t clv_elems = mlipper::util::checked_product(
        "local SPR lane CLVs", nodes, sites, rates, states);
    const size_t pmat_elems = mlipper::util::checked_product(
        "local SPR lane PMATs", nodes, rates, states, states);
    const size_t scaler_per_node = dev.per_rate_scaling
        ? mlipper::util::checked_mul_size(
              sites, rates, "local SPR lane scaler stride")
        : sites;

    size_t bytes = 0;
    const auto add = [&](size_t value, const char* label) {
        bytes = mlipper::util::checked_add_size(bytes, value, label);
    };
    add(mlipper::util::checked_product(
            "local SPR four CLV pools", sizeof(fp_t), clv_elems, 4),
        "local SPR lane bytes");
    add(mlipper::util::checked_product(
            "local SPR four PMAT pools", sizeof(fp_t), pmat_elems, 4),
        "local SPR lane bytes");
    add(mlipper::util::checked_product(
            "local SPR four scaler pools", sizeof(unsigned), nodes,
            scaler_per_node, 4),
        "local SPR lane bytes");
    add(mlipper::util::checked_product(
            "local SPR tip characters", sizeof(uint8_t), tips, sites),
        "local SPR lane bytes");
    add(mlipper::util::checked_allocation_bytes<int>(
            tips, "local SPR tip node IDs"),
        "local SPR lane bytes");
    add(mlipper::util::checked_product(
            "local SPR branch arrays", sizeof(fp_t), nodes, 5),
        "local SPR lane bytes");

    // Cover query buffers, traversal operations, placement scratch, allocator
    // rounding, and small model arrays that are not part of the dominant pools.
    add(bytes / 4, "local SPR lane safety overhead");
    return std::max<size_t>(bytes, 1);
}

int choose_local_spr_scoring_lane_count(const DeviceTree& device)
{
    ensure_device_tree_current_device(
        device, "choose_local_spr_scoring_lane_count");
    size_t free_bytes = 0;
    size_t total_bytes = 0;
    CUDA_CHECK(cudaMemGetInfo(&free_bytes, &total_bytes));
    (void)total_bytes;

    int cpu_limit = 8;
    if (const unsigned hw_threads = std::thread::hardware_concurrency()) {
        cpu_limit = std::min(cpu_limit, static_cast<int>(hw_threads));
    }
    cpu_limit = std::max(1, cpu_limit);

    const size_t lane_bytes = estimate_local_spr_scoring_lane_bytes(device);
    constexpr size_t kMinimumHeadroomBytes = size_t{256} * 1024 * 1024;
    // Ranking also owns a compact batch workspace. Reserve at least one
    // lane-sized allocation for it, in addition to fixed CUDA headroom.
    const size_t reserved_bytes = mlipper::util::checked_add_size(
        kMinimumHeadroomBytes,
        lane_bytes,
        "local SPR GPU memory reserve");
    if (free_bytes <= reserved_bytes || free_bytes - reserved_bytes < lane_bytes) {
        throw std::runtime_error(
            "Local SPR does not have enough free GPU memory for one scoring "
            "lane plus the compact batch workspace and 256 MiB headroom: free=" +
            std::to_string(free_bytes) + " estimated_lane=" +
            std::to_string(lane_bytes) + " reserved=" +
            std::to_string(reserved_bytes));
    }

    const size_t memory_limit = (free_bytes - reserved_bytes) / lane_bytes;
    const int lanes = std::max(
        1,
        std::min(
            cpu_limit,
            memory_limit > static_cast<size_t>(std::numeric_limits<int>::max())
                ? std::numeric_limits<int>::max()
                : static_cast<int>(memory_limit)));
    if (lanes < cpu_limit) {
        std::cout << "Local SPR scoring lanes reduced from " << cpu_limit
                  << " to " << lanes
                  << " by the GPU memory budget (free=" << free_bytes
                  << ", estimated_lane=" << lane_bytes
                  << ", reserved=" << reserved_bytes << ").\n";
    }
    return lanes;
}

void release_local_spr_batch_scoring_workspace(
    LocalSPRBatchScoringWorkspace& workspace,
    cudaStream_t stream)
{
    const bool has_workspace_state =
        workspace.initialized ||
        workspace.candidate_ops.d_ops != nullptr ||
        workspace.subtree_workspace.dev.device_id >= 0 ||
        workspace.subtree_workspace.dev.d_tipchars != nullptr;
    if (!has_workspace_state) {
        return;
    }
    if (workspace.candidate_ops.d_ops != nullptr) {
        free_placement_op_buffer(workspace.candidate_ops, stream);
    }
    release_local_spr_placement_scratch(
        workspace.placement_scratch,
        stream);
    CUDA_CHECK(cudaStreamSynchronize(stream));
    release_subtree_workspace(workspace.subtree_workspace);
    workspace = LocalSPRBatchScoringWorkspace{};
}

void release_local_spr_scoring_lanes(
    std::vector<LocalSPRScoringLane>& lanes)
{
    for (LocalSPRScoringLane& lane : lanes) {
        if (lane.stream != nullptr) {
            release_local_spr_scoring_workspace(lane.workspace, lane.stream);
            if (lane.ready_event != nullptr) {
                (void)cudaEventDestroy(lane.ready_event);
                lane.ready_event = nullptr;
            }
            (void)cudaStreamDestroy(lane.stream);
            lane.stream = nullptr;
        }
    }
    lanes.clear();
}

void release_local_spr_ranking_workspace(
    LocalSPRRankingWorkspace& workspace)
{
    int original_device = -1;
    (void)cudaGetDevice(&original_device);
    if (workspace.owner_device_id >= 0) {
        CUDA_CHECK(cudaSetDevice(workspace.owner_device_id));
    }
    if (workspace.owner_stream != nullptr) {
        release_local_spr_batch_scoring_workspace(
            workspace.batch_scoring_workspace,
            workspace.owner_stream);
    } else {
        workspace.batch_scoring_workspace = LocalSPRBatchScoringWorkspace{};
    }
    release_local_spr_scoring_lanes(workspace.scoring_lanes);
    workspace.owner_stream = nullptr;
    workspace.owner_device_id = -1;
    workspace.scoring_lane_count = 0;
    if (original_device >= 0) (void)cudaSetDevice(original_device);
}

bool local_spr_scoring_lane_matches_context(
    const LocalSPRScoringLane& lane,
    const TopologyRefinementSearchContext& ctx,
    const TreeBuildResult& base_tree,
    int query_capacity)
{
    if (lane.stream == nullptr || lane.ready_event == nullptr) {
        return false;
    }
    const DeviceTree& dev = lane.workspace.subtree_workspace.dev;
    const int tip_count = static_cast<int>(ctx.state.host_packing.tip_node_ids.size());
    return dev.device_id == ctx.state.device.device_id &&
           dev.sites == ctx.sites &&
           dev.states == ctx.states &&
           dev.rate_cats == ctx.rate_cats &&
           dev.per_rate_scaling == ctx.per_rate_scaling &&
           dev.capacity_N >= static_cast<int>(base_tree.nodes.size()) &&
           dev.capacity_tips >= tip_count &&
           dev.query_capacity >= query_capacity;
}

void ensure_local_spr_ranking_workspace(
    const TopologyRefinementSearchContext& ctx,
    const TreeBuildResult& base_tree,
    int scoring_lane_count,
    LocalSPRRankingWorkspace& workspace)
{
    const bool direct_quartet_enabled =
        ctx.move_type == mlipper::TopologyMoveType::NNI &&
        ctx.per_rate_scaling;
    const int desired_lane_count = direct_quartet_enabled
        ? 0
        : scoring_lane_count;
    constexpr int query_capacity = 1;
    const PlacementQueryBatch subtree_query_batch =
        build_local_spr_scoring_query_batch(ctx, query_capacity);
    const bool compatible =
        workspace.scoring_lane_count == desired_lane_count &&
        workspace.scoring_lanes.size() == static_cast<size_t>(desired_lane_count) &&
        std::all_of(
            workspace.scoring_lanes.begin(),
            workspace.scoring_lanes.end(),
            [&](const LocalSPRScoringLane& lane) {
                return local_spr_scoring_lane_matches_context(
                    lane,
                    ctx,
                    base_tree,
                    query_capacity);
            });

    if (!compatible) {
        release_local_spr_ranking_workspace(workspace);
        workspace.owner_stream = ctx.stream;
        workspace.owner_device_id = ctx.state.device.device_id;
        workspace.scoring_lane_count = desired_lane_count;
        ensure_device_tree_current_device(
            ctx.state.device,
            "ensure_local_spr_ranking_workspace");
        CUDA_CHECK(cudaStreamSynchronize(ctx.stream));
        workspace.scoring_lanes.reserve(static_cast<size_t>(desired_lane_count));
        try {
            for (int lane_idx = 0; lane_idx < desired_lane_count; ++lane_idx) {
                workspace.scoring_lanes.emplace_back();
                LocalSPRScoringLane& lane = workspace.scoring_lanes.back();
                CUDA_CHECK(cudaStreamCreateWithFlags(
                    &lane.stream,
                    cudaStreamNonBlocking));
                CUDA_CHECK(cudaEventCreateWithFlags(
                    &lane.ready_event,
                    cudaEventDisableTiming));
                load_subtree_workspace(
                    lane.workspace.subtree_workspace,
                    SubtreeWorkspaceLoadConfig{
                        base_tree,
                        ctx.state.host_packing,
                        ctx.state.eig,
                        ctx.rate_weights,
                        ctx.rate_multipliers,
                        ctx.pi,
                        ctx.sites,
                        ctx.states,
                        ctx.rate_cats,
                        ctx.per_rate_scaling,
                        &subtree_query_batch,
                        true,
                        query_capacity,
                        query_capacity,
                    },
                    lane.stream,
                    "ensure_local_spr_ranking_workspace");
            }
        } catch (...) {
            release_local_spr_ranking_workspace(workspace);
            throw;
        }
        return;
    }

    workspace.owner_stream = ctx.stream;
    workspace.owner_device_id = ctx.state.device.device_id;
    ensure_device_tree_current_device(
        ctx.state.device,
        "ensure_local_spr_ranking_workspace");
    CUDA_CHECK(cudaStreamSynchronize(ctx.stream));
    for (LocalSPRScoringLane& lane : workspace.scoring_lanes) {
        ensure_local_spr_scoring_lane_current_device(
            lane,
            "ensure_local_spr_ranking_workspace");
        CUDA_CHECK(cudaStreamSynchronize(lane.stream));
        load_subtree_workspace(
            lane.workspace.subtree_workspace,
            SubtreeWorkspaceLoadConfig{
                base_tree,
                ctx.state.host_packing,
                ctx.state.eig,
                ctx.rate_weights,
                ctx.rate_multipliers,
                ctx.pi,
                ctx.sites,
                ctx.states,
                ctx.rate_cats,
                ctx.per_rate_scaling,
                &subtree_query_batch,
                true,
                query_capacity,
                query_capacity,
            },
            lane.stream,
            "ensure_local_spr_ranking_workspace");
        lane.workspace.previous_main_pmat_node = -1;
    }
}

void build_compact_local_spr_batch_tree(
    const TopologyRefinementSearchContext& ctx,
    int node_capacity,
    TreeBuildResult& batch_tree,
    HostPacking& batch_host)
{
    const int compact_nodes = std::max(1, node_capacity);
    // The batch workspace needs DeviceTree-shaped storage, not a biological
    // topology. A dummy root plus independent node slots supplies that storage;
    // each candidate op later rewrites its compact target IDs explicitly.
    batch_tree = TreeBuildResult{};
    batch_host = HostPacking{};
    batch_tree.root_id = 0;
    batch_tree.nodes.resize(static_cast<size_t>(compact_nodes));
    batch_tree.postorder.reserve(static_cast<size_t>(compact_nodes));
    batch_tree.preorder.reserve(static_cast<size_t>(compact_nodes));

    for (int node_id = 0; node_id < compact_nodes; ++node_id) {
        TreeNode node;
        node.id = node_id;
        node.parent = (node_id == 0) ? -1 : 0;
        node.left = -1;
        node.right = -1;
        node.is_tip = (node_id == 0);
        node.branch_length_to_parent =
            (node_id == 0) ? fp_t(0) : fp_t(DEFAULT_BRANCH_LENGTH);
        if (node_id == 0) {
            node.name = "__local_spr_batch_dummy_root";
        }
        batch_tree.nodes[static_cast<size_t>(node_id)] = std::move(node);
        batch_tree.postorder.push_back(node_id);
        batch_tree.preorder.push_back(node_id);
    }
    batch_tree.tip_node_by_name["__local_spr_batch_dummy_root"] = 0;

    const size_t node_count = batch_tree.nodes.size();
    batch_host.postorder = batch_tree.postorder;
    batch_host.preorder = batch_tree.preorder;
    batch_host.parent.resize(node_count, -1);
    batch_host.left.resize(node_count, -1);
    batch_host.right.resize(node_count, -1);
    batch_host.is_tip.resize(node_count, 0);
    batch_host.blen.resize(node_count, fp_t(0));
    for (size_t node_idx = 0; node_idx < node_count; ++node_idx) {
        const TreeNode& node = batch_tree.nodes[node_idx];
        batch_host.parent[node_idx] = node.parent;
        batch_host.left[node_idx] = node.left;
        batch_host.right[node_idx] = node.right;
        batch_host.is_tip[node_idx] = node.is_tip ? 1 : 0;
        batch_host.blen[node_idx] = node.branch_length_to_parent;
    }

    batch_host.tip_node_ids = {0};
    batch_host.tipchars.assign(
        ctx.sites,
        static_cast<uint8_t>(ctx.states == 4 ? 15 : 4));

    const size_t pmat_elems = mlipper::util::checked_product(
        "local SPR compact PMAT elements",
        node_count,
        ctx.rate_cats,
        ctx.states,
        ctx.states);
    batch_host.pmats.assign(pmat_elems, fp_t(0));
    batch_host.pmats_mid.assign(pmat_elems, fp_t(0));
    batch_host.pmats_mid_prox.assign(pmat_elems, fp_t(0));
    batch_host.pmats_mid_dist.assign(pmat_elems, fp_t(0));
    batch_host.pattern_weights = ctx.pattern_weights_arg;
    if (batch_host.pattern_weights.empty()) {
        batch_host.pattern_weights.assign(ctx.sites, 1u);
    }
}

void ensure_local_spr_batch_scoring_workspace(
    const TopologyRefinementSearchContext& ctx,
    int node_capacity,
    int query_capacity,
    LocalSPRBatchScoringWorkspace& workspace)
{
    const int compact_node_capacity = std::max(1, node_capacity);
    const int compact_query_capacity = std::max(1, query_capacity);
    const int shrink_threshold = compact_node_capacity >
            std::numeric_limits<int>::max() / 4
        ? std::numeric_limits<int>::max()
        : std::max(compact_node_capacity * 4, 256);
    if (workspace.initialized &&
        workspace.node_capacity >= compact_node_capacity &&
        workspace.batch_capacity >= compact_query_capacity &&
        workspace.node_capacity <= shrink_threshold) {
        return;
    }
    release_local_spr_batch_scoring_workspace(workspace, ctx.stream);

    workspace.node_capacity = compact_node_capacity;
    workspace.batch_capacity = compact_query_capacity;
    build_compact_local_spr_batch_tree(
        ctx,
        compact_node_capacity,
        workspace.tree,
        workspace.host);

    PlacementQueryBatch query_batch;
    query_batch.count = static_cast<size_t>(compact_query_capacity);
    const size_t query_char_count = mlipper::util::checked_product(
        "local SPR compact query characters",
        query_batch.count,
        ctx.sites);
    query_batch.query_chars.assign(
        query_char_count,
        static_cast<uint8_t>(ctx.states == 4 ? 15 : 4));

    load_subtree_workspace(
        workspace.subtree_workspace,
        SubtreeWorkspaceLoadConfig{
            workspace.tree,
            workspace.host,
            ctx.state.eig,
            ctx.rate_weights,
            ctx.rate_multipliers,
            ctx.pi,
            ctx.sites,
            ctx.states,
            ctx.rate_cats,
            ctx.per_rate_scaling,
            &query_batch,
            false,
            compact_query_capacity,
            0,
        },
        ctx.stream,
        "ensure_local_spr_batch_scoring_workspace");
    release_local_spr_embedded_placement_scratch(
        workspace.subtree_workspace.dev,
        ctx.stream);
    workspace.initialized = true;
}

NodeOpInfo make_compact_local_spr_candidate_op(
    NodeOpInfo op,
    int compact_target_id,
    int query_idx)
{
    op.parent_id = 0;
    op.left_id = compact_target_id;
    op.right_id = compact_target_id;
    if (op.dir_tag != static_cast<uint8_t>(CLV_DIR_DOWN_LEFT)) {
        op.dir_tag = static_cast<uint8_t>(CLV_DIR_DOWN_RIGHT);
    }
    op.left_tip_index = -1;
    op.right_tip_index = query_idx;
    op.op_type = static_cast<int>(OP_DOWN_INNER_INNER);
    op.clv_pool = static_cast<uint8_t>(CLV_POOL_DOWN);
    return op;
}

void copy_local_spr_target_context_to_batch_device(
    const DeviceTree& src,
    int src_target_id,
    DeviceTree& dst,
    int dst_target_id,
    cudaStream_t stream)
{
    if (src_target_id < 0 || src_target_id >= src.N ||
        dst_target_id < 0 || dst_target_id >= dst.N) {
        throw std::runtime_error("local SPR compact batch copy received invalid target ids.");
    }
    if (src.sites != dst.sites ||
        src.states != dst.states ||
        src.rate_cats != dst.rate_cats ||
        src.per_rate_scaling != dst.per_rate_scaling) {
        throw std::runtime_error("local SPR compact batch copy received incompatible device shapes.");
    }

    const size_t src_id = static_cast<size_t>(src_target_id);
    const size_t dst_id = static_cast<size_t>(dst_target_id);
    const size_t per_node = src.per_node_elems();
    const size_t scaler_span = src.scaler_elems();

    cuda_copy_device_to_device_async(
        dst.d_clv_up + dst_id * per_node,
        src.d_clv_up + src_id * per_node,
        per_node,
        stream);
    cuda_copy_device_to_device_async(
        dst.d_edge_outside_clv + dst_id * per_node,
        src.d_edge_outside_clv + src_id * per_node,
        per_node,
        stream);
    cuda_copy_device_to_device_async(
        dst.d_site_scaler_up + dst_id * scaler_span,
        src.d_site_scaler_up + src_id * scaler_span,
        scaler_span,
        stream);
    cuda_copy_device_to_device_async(
        dst.d_edge_outside_scaler + dst_id * scaler_span,
        src.d_edge_outside_scaler + src_id * scaler_span,
        scaler_span,
        stream);
}

LocalSPRPlacementPassResult run_local_spr_placement_pass(
    const TopologyRefinementSearchContext& ctx,
    TreeBuildResult& pruned_tree,
    HostPacking& pruned_host,
    const std::vector<int>& node_to_tip,
    const std::vector<int>& candidate_edges,
    int upward_start_node,
    int prune_root_id,
    LocalSPRScoringWorkspace& workspace,
    const std::vector<NodeOpInfo>* already_updated_ops)
{
    LocalSPRPlacementPassResult pass_result;
    workspace.candidate_ops.tuning = ctx.placement_tuning;
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
            workspace.subtree_workspace.dev,
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
                workspace.subtree_workspace.dev,
                pruned_tree,
                pruned_host,
                workspace.tree_ops,
                -1,
                new_update_ops,
                ctx.stream);
        }
    }
    remember_local_spr_dirty_upward_nodes(workspace);

    copy_unscaled_up_clv_to_query_slot(
        ctx.state.device,
        prune_root_id,
        workspace.subtree_workspace.dev,
        0,
        ctx.stream);

    UploadPlacementOps(
        workspace.candidate_ops,
        candidate_ops,
        ctx.stream);
    ensure_local_spr_placement_scratch_capacity(
        workspace.subtree_workspace.dev,
        candidate_ops.size(),
        workspace.placement_scratch,
        ctx.stream);

    DeviceTree query_view =
        make_query_view(workspace.subtree_workspace.dev, 0);
    pass_result.placement_result = EvaluatePlacementCandidates(
        query_view,
        workspace.candidate_ops.d_ops,
        workspace.candidate_ops.num_ops,
        1,
        ctx.stream,
        false,
        workspace.candidate_ops.tuning,
        &workspace.placement_scratch);
    local_spr_assert(
        pass_result.placement_result.top_placements.size() ==
            candidate_ops.size(),
        "local SPR scoring did not return a score for every edge");
    return pass_result;
}

// Return the score obtained by attaching the pruned subtree back where it
// started. Callers use this as the zero-gain reference for every candidate;
// negative infinity means that the unchanged topology was not scored.
double find_unchanged_topology_loglikelihood(
    const RawPlacementResult& placement_result,
    const TreeBuildResult& pruned_tree,
    int baseline_target_id)
{
    const int baseline_edge_child_id =
        canonical_unchanged_topology_edge_child(
            pruned_tree,
            baseline_target_id);
    for (const RawPlacementResult::RankedPlacement& placement :
         placement_result.top_placements) {
        if (placement.target_id == baseline_edge_child_id) {
            return placement.loglikelihood;
        }
    }
    return -std::numeric_limits<double>::infinity();
}

void append_staged_radius_expansion(
    const TopologyRefinementSearchContext& ctx,
    const LocalSPRRepairUnit& unit,
    TreeBuildResult& pruned_tree,
    HostPacking& pruned_host,
    const std::vector<int>& node_to_tip,
    const std::vector<char>& subtree_mask,
    const PruneInfo& prune_info,
    int prune_root_id,
    int inner_search_radius,
    int outer_seed_limit,
    double baseline_logL,
    LocalSPRScoringWorkspace& workspace,
    const std::vector<NodeOpInfo>& already_updated_ops,
    RawPlacementResult& placement_result,
    size_t& enumerated_candidates)
{
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
    if (outer_seed_edges.empty()) {
        return;
    }

    std::vector<int> outer_candidate_edges;
    outer_candidate_edges.reserve(mlipper::util::checked_product(
        "local SPR staged candidate reserve",
        outer_seed_edges.size(),
        static_cast<size_t>(
            std::max(1, ctx.radius - inner_search_radius))));
    std::unordered_set<int> seen_outer_edges;
    for (int seed_edge_child_id : outer_seed_edges) {
        const std::vector<int> expanded_edges =
            collect_outward_edges_from_seed(
                pruned_tree,
                center_dist,
                seed_edge_child_id,
                inner_search_radius,
                ctx.radius,
                prune_info.pruned_id,
                prune_info.free_internal_id);
        for (int edge_child_id : expanded_edges) {
            if (seen_outer_edges.insert(edge_child_id).second) {
                outer_candidate_edges.push_back(edge_child_id);
            }
        }
    }

    std::vector<int> outer_legal_candidate_edges =
        filter_local_spr_candidate_edges(
            pruned_tree,
            unit.envelope_mask,
            subtree_mask,
            outer_candidate_edges);
    outer_legal_candidate_edges.erase(
        std::remove_if(
            outer_legal_candidate_edges.begin(),
            outer_legal_candidate_edges.end(),
            [&](int edge_child) {
                return local_spr_edge_touches_forbidden_node(
                    ctx, pruned_tree, edge_child);
            }),
        outer_legal_candidate_edges.end());
    enumerated_candidates = mlipper::util::checked_add_size(
        enumerated_candidates,
        outer_legal_candidate_edges.size(),
        "local SPR enumerated candidate count");
    if (outer_legal_candidate_edges.empty()) {
        return;
    }

    LocalSPRPlacementPassResult outer_pass =
        run_local_spr_placement_pass(
            ctx,
            pruned_tree,
            pruned_host,
            node_to_tip,
            outer_legal_candidate_edges,
            -1,
            prune_root_id,
            workspace,
            &already_updated_ops);
    placement_result.top_placements.insert(
        placement_result.top_placements.end(),
        std::make_move_iterator(
            outer_pass.placement_result.top_placements.begin()),
        std::make_move_iterator(
            outer_pass.placement_result.top_placements.end()));
    std::sort(
        placement_result.top_placements.begin(),
        placement_result.top_placements.end(),
        [](const RawPlacementResult::RankedPlacement& lhs,
           const RawPlacementResult::RankedPlacement& rhs) {
            return lhs.loglikelihood > rhs.loglikelihood;
        });
}

void append_positive_gain_local_spr_candidates(
    const RawPlacementResult& placement_result,
    const TreeBuildResult& base_tree,
    const TreeBuildResult& pruned_tree,
    const LocalSPRRepairUnit& unit,
    int prune_root_id,
    const std::vector<int>& subtree_nodes,
    double baseline_logL,
    int baseline_edge_child_id,
    int candidate_pool_limit,
    std::vector<LocalSPRCandidateMove>& unit_topk)
{
    for (const RawPlacementResult::RankedPlacement& placement :
         placement_result.top_placements) {
        if (placement.target_id < 0 ||
            placement.target_id >= static_cast<int>(pruned_tree.nodes.size())) {
            continue;
        }
        const int edge_child = placement.target_id;
        if (edge_child == baseline_edge_child_id) {
            continue;
        }
        const int edge_parent =
            pruned_tree.nodes[(size_t)edge_child].parent;
        if (edge_parent < 0) continue;

        const double approx_gain =
            placement.loglikelihood - baseline_logL;
        if (!(approx_gain > 0.0)) {
            continue;
        }

        LocalSPRCandidateMove candidate;
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
            placement.distal_length;
        candidate.subtree_nodes = subtree_nodes;
        const int path_start =
            candidate.old_parent_id >= 0
                ? candidate.old_parent_id
                : prune_root_id;
        candidate.regraft_path_nodes =
            mlipper::divide_and_conquer::build_tree_path_nodes(
                base_tree,
                path_start,
                edge_child);
        keep_local_spr_topk(
            unit_topk,
            std::move(candidate),
            candidate_pool_limit);
    }
}

LocalSPRPruneRootEvaluationResult evaluate_local_spr_prune_root(
    const TopologyRefinementSearchContext& ctx,
    const TreeBuildResult& base_tree,
    const LocalSPRRepairUnit& unit,
    const LocalSPRPruneRootWorkItem& work_item,
    int inner_search_radius,
    bool staged_radius_expansion,
    int outer_seed_limit,
    int candidate_pool_limit,
    LocalSPRScoringWorkspace& workspace)
{
    ensure_device_tree_current_device(
        workspace.subtree_workspace.dev,
        "evaluate_local_spr_prune_root");
    LocalSPRPruneRootEvaluationResult result;
    const int prune_root_id = work_item.prune_root_id;
    const std::vector<int>& subtree_nodes = work_item.subtree_nodes;
    const std::vector<char> subtree_mask =
        build_local_spr_subtree_mask(base_tree, subtree_nodes);

    TreeBuildResult pruned_tree = base_tree;
    PruneInfo prune_info;
    if (!prune_subtree_for_spr(pruned_tree, prune_root_id, prune_info)) {
        return result;
    }

    const std::vector<int>& legal_candidate_edges =
        work_item.legal_inner_candidate_edges;
    result.enumerated_candidates += legal_candidate_edges.size();
    if (legal_candidate_edges.empty()) {
        return result;
    }

    HostPacking pruned_host = ctx.state.host_packing;
    rebuild_host_topology_from_tree_local(pruned_tree, pruned_host);
    pruned_host.pattern_weights = ctx.pattern_weights_arg;
    const std::vector<int> node_to_tip =
        build_local_spr_node_to_tip(pruned_tree, pruned_host);

    int changed_nodes[1] = { prune_info.sibling_id };
    fill_pmats_in_host_packing(
        pruned_tree,
        pruned_host,
        ctx.state.eig,
        ctx.rate_multipliers,
        ctx.states,
        ctx.rate_cats,
        changed_nodes,
        1);
    reload_device_tree_live_data_local_spr(
        workspace.subtree_workspace.dev,
        pruned_tree,
        pruned_host,
        ctx.state.host_packing,
        prune_info.sibling_id,
        workspace.previous_main_pmat_node,
        nullptr,
        ctx.stream);

    restore_local_spr_upward_state(ctx.state.device, workspace, ctx.stream);
    LocalSPRPlacementPassResult inner_pass =
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
    RawPlacementResult placement_result = std::move(inner_pass.placement_result);

    const double baseline_logL = find_unchanged_topology_loglikelihood(
        placement_result,
        pruned_tree,
        prune_info.sibling_id);
    if (!std::isfinite(baseline_logL)) {
        return result;
    }

    if (staged_radius_expansion) {
        append_staged_radius_expansion(
            ctx,
            unit,
            pruned_tree,
            pruned_host,
            node_to_tip,
            subtree_mask,
            prune_info,
            prune_root_id,
            inner_search_radius,
            outer_seed_limit,
            baseline_logL,
            workspace,
            inner_pass.required_downward_update_ops,
            placement_result,
            result.enumerated_candidates);
    }

    append_positive_gain_local_spr_candidates(
        placement_result,
        base_tree,
        pruned_tree,
        unit,
        prune_root_id,
        subtree_nodes,
        baseline_logL,
        prune_info.sibling_id,
        candidate_pool_limit,
        result.candidates);
    return result;
}

void merge_local_spr_prune_root_evaluation(
    LocalSPRPruneRootEvaluationResult&& evaluation,
    int candidate_pool_limit,
    std::vector<LocalSPRCandidateMove>& unit_topk,
    LocalSPRSearchSummary& search_summary)
{
    search_summary.enumerated_candidates += evaluation.enumerated_candidates;
    for (LocalSPRCandidateMove& candidate : evaluation.candidates) {
        keep_local_spr_topk(
            unit_topk,
            std::move(candidate),
            candidate_pool_limit);
    }
}

LocalSPRPreparedPruneRoot prepare_local_spr_prune_root_inner_context(
    const TopologyRefinementSearchContext& ctx,
    const TreeBuildResult& base_tree,
    const LocalSPRPruneRootWorkItem& work_item,
    int lane_index,
    LocalSPRScoringWorkspace& workspace)
{
    ensure_device_tree_current_device(
        workspace.subtree_workspace.dev,
        "prepare_local_spr_prune_root_inner_context");
    LocalSPRPreparedPruneRoot prepared;
    prepared.lane_index = lane_index;
    prepared.prune_root_id = work_item.prune_root_id;
    prepared.subtree_nodes = work_item.subtree_nodes;
    prepared.subtree_mask =
        build_local_spr_subtree_mask(base_tree, prepared.subtree_nodes);

    prepared.pruned_tree = base_tree;
    if (!prune_subtree_for_spr(
            prepared.pruned_tree,
            prepared.prune_root_id,
            prepared.prune_info)) {
        return prepared;
    }

    const std::vector<int>& legal_candidate_edges =
        work_item.legal_inner_candidate_edges;
    prepared.enumerated_candidates += legal_candidate_edges.size();
    if (legal_candidate_edges.empty()) {
        return prepared;
    }

    prepared.pruned_host = ctx.state.host_packing;
    rebuild_host_topology_from_tree_local(prepared.pruned_tree, prepared.pruned_host);
    prepared.pruned_host.pattern_weights = ctx.pattern_weights_arg;
    prepared.node_to_tip =
        build_local_spr_node_to_tip(prepared.pruned_tree, prepared.pruned_host);

    int changed_nodes[1] = { prepared.prune_info.sibling_id };
    fill_pmats_in_host_packing(
        prepared.pruned_tree,
        prepared.pruned_host,
        ctx.state.eig,
        ctx.rate_multipliers,
        ctx.states,
        ctx.rate_cats,
        changed_nodes,
        1);
    reload_device_tree_live_data_local_spr(
        workspace.subtree_workspace.dev,
        prepared.pruned_tree,
        prepared.pruned_host,
        ctx.state.host_packing,
        prepared.prune_info.sibling_id,
        workspace.previous_main_pmat_node,
        nullptr,
        ctx.stream);

    restore_local_spr_upward_state(ctx.state.device, workspace, ctx.stream);
    build_selected_downward_ops(
        prepared.pruned_tree,
        prepared.node_to_tip,
        legal_candidate_edges,
        prepared.inner_candidate_ops);
    build_required_downward_update_ops(
        prepared.pruned_tree,
        prepared.node_to_tip,
        legal_candidate_edges,
        prepared.inner_required_update_ops);
    local_spr_assert(
        !prepared.inner_candidate_ops.empty(),
        "local SPR batched inner scoring produced zero candidate ops");
    local_spr_assert(
        !prepared.inner_required_update_ops.empty(),
        "local SPR batched inner scoring produced zero required downward update ops");

    UpdateTreeClvsAfterPrune(
        workspace.subtree_workspace.dev,
        prepared.pruned_tree,
        prepared.pruned_host,
        workspace.tree_ops,
        prepared.prune_info.grandparent_id,
        prepared.inner_required_update_ops,
        ctx.stream);
    remember_local_spr_dirty_upward_nodes(workspace);
    prepared.valid = true;
    return prepared;
}

// Pack candidate edge contexts from independent pruned trees into one compact
// DeviceTree-shaped workspace. The compact IDs exist only for this launch;
// results are split and remapped to each source topology before returning.
std::vector<RawPlacementResult> score_local_spr_prepared_inner_batch(
    const TopologyRefinementSearchContext& ctx,
    const std::vector<LocalSPRPreparedPruneRoot>& prepared_roots,
    std::vector<LocalSPRScoringLane>& scoring_lanes,
    LocalSPRBatchScoringWorkspace& batch_workspace)
{
    std::vector<RawPlacementResult> split_results(prepared_roots.size());
    if (prepared_roots.empty()) {
        return split_results;
    }

    size_t total_ops = 0;
    for (const LocalSPRPreparedPruneRoot& prepared : prepared_roots) {
        if (prepared.valid) {
            total_ops = mlipper::util::checked_add_size(
                total_ops,
                prepared.inner_candidate_ops.size(),
                "local SPR compact batch operation count");
        }
    }
    if (total_ops == 0) {
        return split_results;
    }
    if (total_ops + 1 > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("local SPR compact batch has too many candidate ops.");
    }
    if (prepared_roots.size() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error(
            "Local SPR compact batch has too many prune roots.");
    }

    const int compact_node_capacity = static_cast<int>(total_ops) + 1;
    ensure_local_spr_batch_scoring_workspace(
        ctx,
        compact_node_capacity,
        static_cast<int>(prepared_roots.size()),
        batch_workspace);

    // Compact target IDs index copied edge contexts; the two maps recover the
    // originating prune-root batch and real target ID after one GPU scoring call.
    std::vector<NodeOpInfo> batched_ops;
    batched_ops.reserve(total_ops);
    std::vector<fp_t> compact_branch_lengths(
        static_cast<size_t>(compact_node_capacity),
        fp_t(0));
    std::vector<int> compact_to_batch_idx(
        static_cast<size_t>(batch_workspace.node_capacity),
        -1);
    std::vector<int> compact_to_local_target_id(
        static_cast<size_t>(batch_workspace.node_capacity),
        -1);

    for (size_t batch_idx = 0; batch_idx < prepared_roots.size(); ++batch_idx) {
        const LocalSPRPreparedPruneRoot& prepared = prepared_roots[batch_idx];
        if (!prepared.valid) continue;
        local_spr_assert(
            prepared.lane_index >= 0 &&
                prepared.lane_index < static_cast<int>(scoring_lanes.size()),
            "local SPR batched scoring has invalid lane index");
        LocalSPRScoringLane& lane =
            scoring_lanes[static_cast<size_t>(prepared.lane_index)];
        CUDA_CHECK(cudaStreamWaitEvent(ctx.stream, lane.ready_event, 0));

        copy_unscaled_up_clv_to_query_slot(
            ctx.state.device,
            prepared.prune_root_id,
            batch_workspace.subtree_workspace.dev,
            static_cast<int>(batch_idx),
            ctx.stream);

        for (NodeOpInfo op : prepared.inner_candidate_ops) {
            const int local_target_id = local_spr_target_id_from_op(op);
            local_spr_assert(
                local_target_id >= 0 &&
                    local_target_id <
                        lane.workspace.subtree_workspace.dev.N,
                "local SPR compact batch received invalid local target id");
            const int compact_target_id = static_cast<int>(batched_ops.size()) + 1;
            compact_branch_lengths[static_cast<size_t>(compact_target_id)] =
                prepared.pruned_tree.nodes[static_cast<size_t>(local_target_id)]
                    .branch_length_to_parent;
            copy_local_spr_target_context_to_batch_device(
                lane.workspace.subtree_workspace.dev,
                local_target_id,
                batch_workspace.subtree_workspace.dev,
                compact_target_id,
                ctx.stream);
            compact_to_batch_idx[static_cast<size_t>(compact_target_id)] =
                static_cast<int>(batch_idx);
            compact_to_local_target_id[static_cast<size_t>(compact_target_id)] =
                local_target_id;
            batched_ops.push_back(
                make_compact_local_spr_candidate_op(
                    op,
                    compact_target_id,
                    static_cast<int>(batch_idx)));
        }
    }

    if (batched_ops.empty()) {
        return split_results;
    }

    CUDA_CHECK(cudaMemcpyAsync(
        batch_workspace.subtree_workspace.dev.d_blen,
        compact_branch_lengths.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            compact_branch_lengths.size(),
            "local SPR compact branch-length upload"),
        cudaMemcpyHostToDevice,
        ctx.stream));

    PlacementTuningConfig batch_tuning = ctx.placement_tuning;
    batch_tuning.export_placement_topk =
        std::max(static_cast<int>(batched_ops.size()), batch_tuning.export_placement_topk);
    batch_workspace.candidate_ops.tuning = batch_tuning;
    UploadPlacementOps(
        batch_workspace.candidate_ops,
        batched_ops,
        ctx.stream);
    ensure_local_spr_placement_scratch_capacity(
        batch_workspace.subtree_workspace.dev,
        batched_ops.size(),
        batch_workspace.placement_scratch,
        ctx.stream);

    RawPlacementResult batch_result = EvaluatePlacementCandidates(
        batch_workspace.subtree_workspace.dev,
        batch_workspace.candidate_ops.d_ops,
        batch_workspace.candidate_ops.num_ops,
        1,
        ctx.stream,
        false,
        batch_tuning,
        &batch_workspace.placement_scratch);
    local_spr_assert(
        batch_result.top_placements.size() == batched_ops.size(),
        "local SPR batched inner scoring did not return every candidate");

    for (RawPlacementResult::RankedPlacement placement :
         batch_result.top_placements) {
        if (placement.target_id <= 0 ||
            placement.target_id >= static_cast<int>(compact_to_batch_idx.size())) {
            continue;
        }
        const int batch_idx =
            compact_to_batch_idx[static_cast<size_t>(placement.target_id)];
        if (batch_idx < 0 ||
            batch_idx >= static_cast<int>(prepared_roots.size())) {
            continue;
        }
        const int local_target_id =
            compact_to_local_target_id[static_cast<size_t>(placement.target_id)];
        if (local_target_id < 0) {
            continue;
        }
        placement.target_id = local_target_id;
        split_results[static_cast<size_t>(batch_idx)].top_placements.push_back(
            placement);
    }

    for (size_t batch_idx = 0; batch_idx < split_results.size(); ++batch_idx) {
        RawPlacementResult& result = split_results[batch_idx];
        if (result.top_placements.empty()) continue;
        result.target_id = result.top_placements.front().target_id;
        result.loglikelihood = result.top_placements.front().loglikelihood;
        result.distal_length = result.top_placements.front().distal_length;
        result.pendant_length = result.top_placements.front().pendant_length;
    }

    return split_results;
}

// Score NNI alternatives without materializing one pruned DeviceTree per root.
// The resident source tree supplies upward/downward messages; compact target
// contexts reconstruct the two affected edge environments and one batched
// placement call returns scores that are remapped to each prune-root result.
std::vector<LocalSPRPruneRootEvaluationResult> score_direct_nni_batch(
    const TopologyRefinementSearchContext& ctx,
    const TreeBuildResult& base_tree,
    const LocalSPRRepairUnit& unit,
    const std::vector<LocalSPRPruneRootWorkItem>& work_items,
    size_t range_begin,
    size_t range_end,
    int candidate_pool_limit,
    LocalSPRBatchScoringWorkspace& batch_workspace)
{
    if (range_begin > range_end || range_end > work_items.size()) {
        throw std::invalid_argument("Invalid direct NNI scoring range.");
    }
    const DeviceTree& source_dev = ctx.state.device;
    const size_t root_count = range_end - range_begin;
    std::vector<LocalSPRPruneRootEvaluationResult> results(root_count);
    if (root_count == 0) return results;
    if (!ctx.per_rate_scaling) {
        throw std::runtime_error(
            "direct NNI scorer currently requires per-rate scaling");
    }

    const size_t compact_node_count = mlipper::util::checked_add_size(
        mlipper::util::checked_product(
            "direct NNI compact node count", root_count, size_t{2}),
        size_t{1},
        "direct NNI compact node count");
    if (compact_node_count >
            static_cast<size_t>(std::numeric_limits<int>::max()) ||
        root_count > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("Direct NNI batch exceeds integer indexing.");
    }
    // Each prune root contributes a baseline and one alternate quartet context.
    ensure_local_spr_batch_scoring_workspace(
        ctx,
        static_cast<int>(compact_node_count),
        static_cast<int>(root_count),
        batch_workspace);
    std::vector<DirectNNIContextOp> context_ops;
    std::vector<NodeOpInfo> placement_ops;
    std::vector<fp_t> branch_lengths(
        static_cast<size_t>(batch_workspace.node_capacity), fp_t(0));
    std::vector<int> compact_to_root(
        static_cast<size_t>(batch_workspace.node_capacity), -1);
    std::vector<int> compact_to_source_target(
        static_cast<size_t>(batch_workspace.node_capacity), -1);
    const size_t operation_capacity = mlipper::util::checked_product(
        "direct NNI operation reserve", root_count, size_t{2});
    context_ops.reserve(operation_capacity);
    placement_ops.reserve(operation_capacity);

    for (size_t local_idx = 0; local_idx < root_count; ++local_idx) {
        const LocalSPRPruneRootWorkItem& item = work_items[range_begin + local_idx];
        LocalSPRPruneRootEvaluationResult& result = results[local_idx];
        result.enumerated_candidates = item.legal_inner_candidate_edges.size();
        const int a = item.prune_root_id;
        if (a < 0 || a >= static_cast<int>(base_tree.nodes.size())) continue;
        const int u = base_tree.nodes[static_cast<size_t>(a)].parent;
        if (u < 0) continue;
        const TreeNode& central_child = base_tree.nodes[static_cast<size_t>(u)];
        const int b = central_child.left == a ? central_child.right :
            (central_child.right == a ? central_child.left : -1);
        const int p = central_child.parent;
        if (b < 0 || p < 0 || p == base_tree.root_id) continue;
        const TreeNode& central_parent = base_tree.nodes[static_cast<size_t>(p)];
        const int d = central_parent.left == u ? central_parent.right :
            (central_parent.right == u ? central_parent.left : -1);
        if (d < 0) continue;

        copy_unscaled_up_clv_to_query_slot(
            source_dev,
            a,
            batch_workspace.subtree_workspace.dev,
            static_cast<int>(local_idx),
            ctx.stream);

        const int baseline_dst = static_cast<int>(context_ops.size()) + 1;
        context_ops.push_back(DirectNNIContextOp{b, u, -1, -1, -1, baseline_dst});
        branch_lengths[static_cast<size_t>(baseline_dst)] =
            base_tree.nodes[static_cast<size_t>(b)].branch_length_to_parent +
            central_child.branch_length_to_parent;
        compact_to_root[static_cast<size_t>(baseline_dst)] = static_cast<int>(local_idx);
        compact_to_source_target[static_cast<size_t>(baseline_dst)] = b;

        const int alternative_dst = static_cast<int>(context_ops.size()) + 1;
        context_ops.push_back(DirectNNIContextOp{d, -1, p, b, u, alternative_dst});
        branch_lengths[static_cast<size_t>(alternative_dst)] =
            base_tree.nodes[static_cast<size_t>(d)].branch_length_to_parent;
        compact_to_root[static_cast<size_t>(alternative_dst)] = static_cast<int>(local_idx);
        compact_to_source_target[static_cast<size_t>(alternative_dst)] = d;

        NodeOpInfo op{};
        op.dir_tag = static_cast<uint8_t>(CLV_DIR_DOWN_LEFT);
        placement_ops.push_back(make_compact_local_spr_candidate_op(
            op, baseline_dst, static_cast<int>(local_idx)));
        placement_ops.push_back(make_compact_local_spr_candidate_op(
            op, alternative_dst, static_cast<int>(local_idx)));
    }
    if (placement_ops.empty()) return results;

    build_direct_nni_target_contexts(
        source_dev,
        batch_workspace.subtree_workspace.dev,
        context_ops,
        ctx.stream);
    CUDA_CHECK(cudaMemcpyAsync(
        batch_workspace.subtree_workspace.dev.d_blen,
        branch_lengths.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            branch_lengths.size(), "direct NNI branch-length upload"),
        cudaMemcpyHostToDevice,
        ctx.stream));
    PlacementTuningConfig tuning = ctx.placement_tuning;
    tuning.export_placement_topk = static_cast<int>(placement_ops.size());
    batch_workspace.candidate_ops.tuning = tuning;
    UploadPlacementOps(batch_workspace.candidate_ops, placement_ops, ctx.stream);
    ensure_local_spr_placement_scratch_capacity(
        batch_workspace.subtree_workspace.dev,
        placement_ops.size(),
        batch_workspace.placement_scratch,
        ctx.stream);
    RawPlacementResult batch_result = EvaluatePlacementCandidates(
        batch_workspace.subtree_workspace.dev,
        batch_workspace.candidate_ops.d_ops,
        batch_workspace.candidate_ops.num_ops,
        1,
        ctx.stream,
        false,
        tuning,
        &batch_workspace.placement_scratch);

    std::vector<RawPlacementResult> per_root(root_count);
    for (auto placement : batch_result.top_placements) {
        const int compact = placement.target_id;
        if (compact <= 0 || compact >= static_cast<int>(compact_to_root.size())) continue;
        const int local_idx = compact_to_root[static_cast<size_t>(compact)];
        if (local_idx < 0) continue;
        placement.target_id = compact_to_source_target[static_cast<size_t>(compact)];
        per_root[static_cast<size_t>(local_idx)].top_placements.push_back(placement);
    }

    for (size_t local_idx = 0; local_idx < root_count; ++local_idx) {
        const auto& item = work_items[range_begin + local_idx];
        const int u = base_tree.nodes[static_cast<size_t>(item.prune_root_id)].parent;
        const TreeNode& central_child = base_tree.nodes[static_cast<size_t>(u)];
        const int b = central_child.left == item.prune_root_id
            ? central_child.right : central_child.left;
        const double baseline = find_unchanged_topology_loglikelihood(
            per_root[local_idx], base_tree, b);
        if (!std::isfinite(baseline)) continue;
        append_positive_gain_local_spr_candidates(
            per_root[local_idx],
            base_tree,
            base_tree,
            unit,
            item.prune_root_id,
            item.subtree_nodes,
            baseline,
            b,
            candidate_pool_limit,
            results[local_idx].candidates);
    }
    return results;
}

// Dispatch one contiguous work range through the specialized direct-NNI path
// or the general per-lane SPR path. Both branches merge results through the
// same deterministic per-unit top-k and accounting contract.
void evaluate_local_spr_work_range_batched_inner(
    const TopologyRefinementSearchContext& ctx,
    const TreeBuildResult& base_tree,
    const LocalSPRRepairUnit& unit,
    const std::vector<LocalSPRPruneRootWorkItem>& work_items,
    size_t range_begin,
    size_t range_end,
    int inner_search_radius,
    bool staged_radius_expansion,
    int outer_seed_limit,
    int candidate_pool_limit,
    std::vector<LocalSPRScoringLane>& scoring_lanes,
    tbb::task_arena& task_arena,
    LocalSPRBatchScoringWorkspace& batch_workspace,
    std::vector<LocalSPRCandidateMove>& unit_topk,
    LocalSPRSearchSummary& search_summary)
{
    if (ctx.move_type == mlipper::TopologyMoveType::NNI &&
        ctx.per_rate_scaling) {
        std::vector<LocalSPRPruneRootEvaluationResult> batch_results =
            score_direct_nni_batch(
                ctx,
                base_tree,
                unit,
                work_items,
                range_begin,
                range_end,
                candidate_pool_limit,
                batch_workspace);
        for (LocalSPRPruneRootEvaluationResult& evaluation : batch_results) {
            merge_local_spr_prune_root_evaluation(
                std::move(evaluation),
                candidate_pool_limit,
                unit_topk,
                search_summary);
        }
        return;
    }

    local_spr_assert(
        !scoring_lanes.empty(),
        "local SPR batched scoring requested with zero scoring lanes");

    for (size_t batch_begin = range_begin;
         batch_begin < range_end;
         batch_begin += scoring_lanes.size()) {
        const size_t batch_size =
            std::min(scoring_lanes.size(), range_end - batch_begin);
        if (batch_size == 1) {
            TopologyRefinementSearchContext lane_ctx = ctx;
            lane_ctx.stream = scoring_lanes[0].stream;
            merge_local_spr_prune_root_evaluation(
                evaluate_local_spr_prune_root(
                    lane_ctx,
                    base_tree,
                    unit,
                    work_items[batch_begin],
                    inner_search_radius,
                    staged_radius_expansion,
                    outer_seed_limit,
                    candidate_pool_limit,
                    scoring_lanes[0].workspace),
                candidate_pool_limit,
                unit_topk,
                search_summary);
            continue;
        }

        std::vector<LocalSPRPreparedPruneRoot> prepared_roots(batch_size);
        task_arena.execute([&]() {
            tbb::parallel_for(size_t{0}, batch_size, [&](size_t lane_idx) {
                const size_t item_idx = batch_begin + lane_idx;
                LocalSPRScoringLane& lane = scoring_lanes[lane_idx];
                ensure_local_spr_scoring_lane_current_device(
                    lane,
                    "evaluate_local_spr_work_range_batched_inner");
                TopologyRefinementSearchContext lane_ctx = ctx;
                lane_ctx.stream = lane.stream;
                prepared_roots[lane_idx] =
                    prepare_local_spr_prune_root_inner_context(
                        lane_ctx,
                        base_tree,
                        work_items[item_idx],
                        static_cast<int>(lane_idx),
                        lane.workspace);
                if (prepared_roots[lane_idx].valid) {
                    CUDA_CHECK(cudaEventRecord(lane.ready_event, lane.stream));
                }
            });
        });

        std::vector<RawPlacementResult> inner_results =
            score_local_spr_prepared_inner_batch(
                ctx,
                prepared_roots,
                scoring_lanes,
                batch_workspace);

        for (size_t batch_idx = 0; batch_idx < prepared_roots.size(); ++batch_idx) {
            LocalSPRPreparedPruneRoot& prepared = prepared_roots[batch_idx];
            LocalSPRPruneRootEvaluationResult evaluation;
            evaluation.enumerated_candidates = prepared.enumerated_candidates;
            if (!prepared.valid) {
                merge_local_spr_prune_root_evaluation(
                    std::move(evaluation),
                    candidate_pool_limit,
                    unit_topk,
                    search_summary);
                continue;
            }

            RawPlacementResult placement_result =
                std::move(inner_results[batch_idx]);
            const double baseline_logL = find_unchanged_topology_loglikelihood(
                placement_result,
                prepared.pruned_tree,
                prepared.prune_info.sibling_id);
            if (!std::isfinite(baseline_logL)) {
                merge_local_spr_prune_root_evaluation(
                    std::move(evaluation),
                    candidate_pool_limit,
                    unit_topk,
                    search_summary);
                continue;
            }

            if (staged_radius_expansion) {
                LocalSPRScoringLane& lane =
                    scoring_lanes[static_cast<size_t>(prepared.lane_index)];
                TopologyRefinementSearchContext lane_ctx = ctx;
                lane_ctx.stream = lane.stream;
                append_staged_radius_expansion(
                    lane_ctx,
                    unit,
                    prepared.pruned_tree,
                    prepared.pruned_host,
                    prepared.node_to_tip,
                    prepared.subtree_mask,
                    prepared.prune_info,
                    prepared.prune_root_id,
                    inner_search_radius,
                    outer_seed_limit,
                    baseline_logL,
                    lane.workspace,
                    prepared.inner_required_update_ops,
                    placement_result,
                    evaluation.enumerated_candidates);
            }

            append_positive_gain_local_spr_candidates(
                placement_result,
                base_tree,
                prepared.pruned_tree,
                unit,
                prepared.prune_root_id,
                prepared.subtree_nodes,
                baseline_logL,
                prepared.prune_info.sibling_id,
                candidate_pool_limit,
                evaluation.candidates);
            merge_local_spr_prune_root_evaluation(
                std::move(evaluation),
                candidate_pool_limit,
                unit_topk,
                search_summary);
        }
    }
}

// Execute one complete candidate-search stage. CPU work first enumerates legal
// prune roots and inner-radius edges per repair unit. GPU lanes then score those
// edges in batches; promising seeds may expand to the outer radius. Only each
// unit's deterministic top-k survives into the globally ranked result.
std::vector<LocalSPRCandidateMove> rank_local_spr_candidates(
    const TopologyRefinementSearchContext& ctx,
    const TreeBuildResult& base_tree,
    const std::vector<LocalSPRRepairUnit>& repair_units,
    LocalSPRRankingWorkspace& ranking_workspace,
    tbb::task_arena& task_arena,
    int scoring_lane_count,
    LocalSPRSearchSummary& search_summary)
{
    std::vector<LocalSPRCandidateMove> ranked_candidates;

    ensure_local_spr_ranking_workspace(
        ctx,
        base_tree,
        scoring_lane_count,
        ranking_workspace);
    std::vector<LocalSPRScoringLane>& scoring_lanes =
        ranking_workspace.scoring_lanes;
    LocalSPRBatchScoringWorkspace& batch_scoring_workspace =
        ranking_workspace.batch_scoring_workspace;

    const int inner_search_radius = std::min(ctx.radius, 2);
    const bool staged_radius_expansion =
        ctx.radius > inner_search_radius;
    const int local_spr_candidate_pool_limit = ctx.topk_per_unit;
    const int local_spr_outer_seed_limit = ctx.topk_per_unit;
    std::vector<std::vector<LocalSPRPruneRootWorkItem>>
        prepared_unit_work_items(repair_units.size());
    task_arena.execute([&]() {
        tbb::parallel_for(
            size_t{0}, repair_units.size(), [&](size_t unit_idx) {
                prepared_unit_work_items[unit_idx] =
                    prepare_local_spr_unit_work_items(
                        base_tree,
                        repair_units[unit_idx],
                        inner_search_radius,
                        ctx);
            });
    });
    if (ctx.move_type != mlipper::TopologyMoveType::NNI &&
        !ctx.forbidden_regraft_node_ids.empty()) {
        for (auto& unit_items : prepared_unit_work_items) {
            unit_items.erase(
                std::remove_if(
                    unit_items.begin(),
                    unit_items.end(),
                    [&](const LocalSPRPruneRootWorkItem& item) {
                        return std::any_of(
                            item.subtree_nodes.begin(),
                            item.subtree_nodes.end(),
                            [&](int node_id) {
                                return std::find(
                                    ctx.forbidden_regraft_node_ids.begin(),
                                    ctx.forbidden_regraft_node_ids.end(),
                                    node_id) !=
                                    ctx.forbidden_regraft_node_ids.end();
                            });
                    }),
                unit_items.end());
            for (LocalSPRPruneRootWorkItem& item : unit_items) {
                auto& edges = item.legal_inner_candidate_edges;
                edges.erase(
                    std::remove_if(
                        edges.begin(),
                        edges.end(),
                        [&](int edge_child) {
                            return local_spr_edge_touches_forbidden_node(
                                ctx, base_tree, edge_child);
                        }),
                    edges.end());
            }
            unit_items.erase(
                std::remove_if(
                    unit_items.begin(),
                    unit_items.end(),
                    [](const LocalSPRPruneRootWorkItem& item) {
                        return item.legal_inner_candidate_edges.empty();
                    }),
                unit_items.end());
        }
    }
    ranked_candidates.reserve(mlipper::util::checked_product(
        "local SPR ranked candidate reserve",
        repair_units.size(),
        static_cast<size_t>(
            std::max(1, local_spr_candidate_pool_limit))));
    for (size_t unit_idx = 0; unit_idx < repair_units.size(); ++unit_idx) {
        const LocalSPRRepairUnit& unit = repair_units[unit_idx];
        std::vector<LocalSPRCandidateMove> unit_topk;
        std::vector<LocalSPRPruneRootWorkItem>& prune_root_work_items =
            prepared_unit_work_items[unit_idx];
        if (prune_root_work_items.empty()) {
            continue;
        }

        auto evaluate_range = [&](size_t range_begin, size_t range_end) {
            evaluate_local_spr_work_range_batched_inner(
                ctx,
                base_tree,
                unit,
                prune_root_work_items,
                range_begin,
                range_end,
                inner_search_radius,
                staged_radius_expansion,
                local_spr_outer_seed_limit,
                local_spr_candidate_pool_limit,
                scoring_lanes,
                task_arena,
                batch_scoring_workspace,
                unit_topk,
                search_summary);
        };

        evaluate_range(0, prune_root_work_items.size());
        search_summary.retained_candidates = mlipper::util::checked_add_size(
            search_summary.retained_candidates,
            unit_topk.size(),
            "local SPR retained candidate count");
        ranked_candidates.insert(
            ranked_candidates.end(),
            std::make_move_iterator(unit_topk.begin()),
            std::make_move_iterator(unit_topk.end()));
    }
    return ranked_candidates;
}
