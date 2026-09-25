#pragma once

#include <string>
#include <vector>

#include "io/parse_file.hpp"
#include "placement/placement.cuh"
#include "tree.hpp"

namespace mlipper {
namespace divide_and_conquer {

// The two node IDs that identify one undirected tree edge.
struct TreeEdgeEndpoints {
    int endpoint_a = -1;
    int endpoint_b = -1;
};

struct SubtreePlan {
    // Global node IDs selected for the bounded problem. The local/global maps
    // become valid only after build_subtree_from_plan materializes the tree.
    std::vector<int> skeleton_nodes;
    std::vector<int> selected_tip_ids;
    std::vector<int> subtree_node_ids;
    std::vector<int> local_to_global;
    std::vector<int> global_to_local;
};

struct DirectionalBoundaryPort {
    // local_node_id is a virtual tip; outside_global_node_id identifies the
    // full-tree directional CLV copied into that tip slot.
    int local_node_id = -1;
    int outside_global_node_id = -1;
};

struct PreparedDirectionalSubtree {
    // These values are one coherent local representation. Repacking subtree
    // independently would invalidate boundary_ports and local_to_global.
    TreeBuildResult subtree;
    parse::Alignment alignment;
    HostPacking host_packing;
    std::vector<int> local_to_global;
    std::vector<DirectionalBoundaryPort> boundary_ports;
};

// Return the inclusive simple path between two current node IDs.
std::vector<int> build_tree_path_nodes(
    const TreeBuildResult& tree,
    int start,
    int end);

// Compute undirected edge distance to the nearest marked node. Unreachable
// entries remain -1; source_mask must be indexed like tree.nodes.
std::vector<int> multi_source_bfs_distances(
    const TreeBuildResult& tree,
    const std::vector<char>& source_mask);

// Connect anchor edges into a skeleton, then expand outward by complete BFS
// levels until tip_budget is reached. This function plans but does not renumber.
SubtreePlan build_subtree_plan_from_anchors(
    const TreeBuildResult& tree,
    const std::vector<TreeEdgeEndpoints>& anchors,
    int tip_budget);

// Build the exact local sector surrounding central internal edges for NNI.
SubtreePlan build_nni_sector_plan(
    const TreeBuildResult& tree,
    const std::vector<TreeEdgeEndpoints>& central_edges);

// Materialize a dense local tree and populate both local/global ID maps.
TreeBuildResult build_subtree_from_plan(
    const TreeBuildResult& tree,
    const SubtreePlan& plan);

// Build a bounded local likelihood problem. Cut full-tree edges become virtual
// boundary tips whose CLVs represent the omitted side of each edge.
PreparedDirectionalSubtree prepare_anchor_subtree_with_directional_boundaries(
    const TreeBuildResult& tree,
    const std::vector<TreeEdgeEndpoints>& anchors,
    int tip_budget,
    const std::vector<std::string>& full_msa_tip_names,
    const std::vector<std::string>& full_msa_rows,
    size_t sites,
    int states,
    const EigResult& eig,
    const std::vector<double>& rate_multipliers,
    int rate_cats);

PreparedDirectionalSubtree prepare_nni_sector_with_directional_boundaries(
    const TreeBuildResult& tree,
    const std::vector<TreeEdgeEndpoints>& central_edges,
    const std::vector<std::string>& full_msa_tip_names,
    const std::vector<std::string>& full_msa_rows,
    size_t sites,
    int states,
    const EigResult& eig,
    const std::vector<double>& rate_multipliers,
    int rate_cats);

// Encode query characters into the same state representation used by DeviceTree.
PlacementQueryBatch build_query_batch_from_placement_queries(
    const std::vector<mlipper::SequenceRecord>& placement_queries,
    size_t sites,
    int states);

} // namespace divide_and_conquer
} // namespace mlipper
