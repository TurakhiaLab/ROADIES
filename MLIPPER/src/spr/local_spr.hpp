#pragma once

#include <functional>
#include <string>
#include <unordered_set>
#include <vector>

#include "tree/divide_and_conquer.hpp"
#include "tree/tree.hpp"

namespace mlipper {
struct PlacementResult;

enum class TopologyMoveType {
    SPR,
    NNI,
};
}

struct LocalSPRSessionWorkspaceSet;

// Non-owning inputs and mutable session state for one refinement invocation.
// The caller must keep every referenced object alive until the run completes.
struct TopologyRefinementRunContext {
    TopologyRefinementRunContext(
        TopologyRefinementState& state_arg,
        PlacementOpBuffer& placement_ops_arg,
        cudaStream_t stream_arg,
        std::vector<std::string>& current_names_arg,
        std::vector<std::string>& current_rows_arg,
        std::string& current_tree_newick_arg,
        const std::vector<unsigned>& pattern_weights_arg,
        const std::vector<double>& rate_weights_arg,
        const std::vector<double>& rate_multipliers_arg,
        const std::vector<double>& pi_arg,
        size_t sites_arg,
        int states_arg,
        int rate_cats_arg,
        bool per_rate_scaling_arg,
        const std::vector<
            mlipper::divide_and_conquer::TreeEdgeEndpoints>& anchors_arg,
        LocalSPRSessionWorkspaceSet& session_workspaces_arg)
        : state(state_arg),
          placement_ops(placement_ops_arg),
          stream(stream_arg),
          current_names(current_names_arg),
          current_rows(current_rows_arg),
          current_tree_newick(current_tree_newick_arg),
          pattern_weights_arg(pattern_weights_arg),
          rate_weights(rate_weights_arg),
          rate_multipliers(rate_multipliers_arg),
          pi(pi_arg),
          sites(sites_arg),
          states(states_arg),
          rate_cats(rate_cats_arg),
          per_rate_scaling(per_rate_scaling_arg),
          anchors(anchors_arg),
          session_workspaces(session_workspaces_arg)
    {
    }

    TopologyRefinementState& state;
    PlacementOpBuffer& placement_ops;
    cudaStream_t stream = nullptr;
    std::vector<std::string>& current_names;
    std::vector<std::string>& current_rows;
    std::string& current_tree_newick;
    const std::vector<unsigned>& pattern_weights_arg;
    const std::vector<double>& rate_weights;
    const std::vector<double>& rate_multipliers;
    const std::vector<double>& pi;
    size_t sites = 0;
    int states = 0, rate_cats = 0;
    bool per_rate_scaling = false;
    int radius = 4;
    int cluster_threshold = 0;
    int topk_per_unit = 0;
    int rounds = 0;
    const std::vector<mlipper::divide_and_conquer::TreeEdgeEndpoints>& anchors;
    LocalSPRSessionWorkspaceSet& session_workspaces;
    // Invoked before an accepted candidate replaces state.tree. This permits
    // workflow-specific bookkeeping against the pre-move topology.
    std::function<void(
        const TreeBuildResult&,
        int,
        int,
        double,
        double)> on_accept_move;
    int accepted_move_count = 0;
    mlipper::TopologyMoveType move_type = mlipper::TopologyMoveType::SPR;
    std::vector<int> directional_boundary_node_ids;
    std::unordered_set<int> directional_boundary_tip_node_ids;
    std::unordered_set<int> allowed_nni_central_edge_child_ids;
    int symmetric_branch_sweeps = 0;
};

void run_topology_refinement(TopologyRefinementRunContext& ctx);
