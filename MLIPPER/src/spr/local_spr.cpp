#include "local_spr.hpp"
#include "local_spr_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <unordered_set>
#include <vector>

#include "io/tree_newick.hpp"
#include "placement/placement.cuh"
#include "likelihood/root_loglikelihood.cuh"
#include "util/mlipper_util.h"
#include "util/checked_size.hpp"

namespace mltreeio = mlipper::treeio;

namespace {

struct LocalSPRLikelihoodEvaluation {
    double total = -std::numeric_limits<double>::infinity();
    TreeBuildResult optimized_tree;
};

void install_local_spr_virtual_boundary_messages(
    const TopologyRefinementRunContext& ctx,
    DeviceTree& destination,
    cudaStream_t stream)
{
    if (ctx.directional_boundary_node_ids.empty()) return;
    const DeviceTree& source = ctx.state.device;
    const size_t per_node = source.per_node_elems();
    const size_t scaler_span = source.scaler_elems();
    for (int node_id : ctx.directional_boundary_node_ids) {
        if (node_id < 0 || node_id >= source.N || node_id >= destination.N) {
            throw std::runtime_error(
                "local SPR virtual boundary node id is out of range.");
        }
        if (source.device_id == destination.device_id) {
            CUDA_CHECK(cudaMemcpyAsync(
                destination.d_clv_up + static_cast<size_t>(node_id) * per_node,
                source.d_clv_up + static_cast<size_t>(node_id) * per_node,
                sizeof(fp_t) * per_node,
                cudaMemcpyDeviceToDevice,
                stream));
            CUDA_CHECK(cudaMemcpyAsync(
                destination.d_site_scaler_up +
                    static_cast<size_t>(node_id) * scaler_span,
                source.d_site_scaler_up +
                    static_cast<size_t>(node_id) * scaler_span,
                sizeof(unsigned) * scaler_span,
                cudaMemcpyDeviceToDevice,
                stream));
        } else {
            CUDA_CHECK(cudaMemcpyPeerAsync(
                destination.d_clv_up + static_cast<size_t>(node_id) * per_node,
                destination.device_id,
                source.d_clv_up + static_cast<size_t>(node_id) * per_node,
                source.device_id,
                sizeof(fp_t) * per_node,
                stream));
            CUDA_CHECK(cudaMemcpyPeerAsync(
                destination.d_site_scaler_up +
                    static_cast<size_t>(node_id) * scaler_span,
                destination.device_id,
                source.d_site_scaler_up +
                    static_cast<size_t>(node_id) * scaler_span,
                source.device_id,
                sizeof(unsigned) * scaler_span,
                stream));
        }
    }
}

std::string write_local_spr_context_tree_newick(
    const TopologyRefinementRunContext& ctx,
    const TreeBuildResult& tree)
{
    if (ctx.directional_boundary_node_ids.empty()) {
        return mltreeio::write_tree_to_newick_string(tree);
    }
    TreeBuildResult output_tree = tree;
    for (int node_id : ctx.directional_boundary_node_ids) {
        if (node_id >= 0 &&
            node_id < static_cast<int>(output_tree.nodes.size())) {
            output_tree.nodes[static_cast<size_t>(node_id)].is_tip = true;
        }
    }
    return mltreeio::write_tree_to_newick_string(output_tree);
}

template <typename T>
void cuda_free_if_set(T*& ptr)
{
    if (ptr == nullptr) {
        return;
    }
    CUDA_CHECK(cudaFree(ptr));
    ptr = nullptr;
}

void release_local_spr_eval_workspace(
    LocalSPREvalWorkspace& workspace,
    cudaStream_t fallback_stream)
{
    const cudaStream_t stream = workspace.stream != nullptr
        ? workspace.stream
        : fallback_stream;
    if (workspace.device_id >= 0) {
        CUDA_CHECK(cudaSetDevice(workspace.device_id));
    }
    const bool has_workspace_state =
        workspace.initialized ||
        workspace.ops.d_ops != nullptr ||
        workspace.subtree_workspace.dev.device_id >= 0 ||
        workspace.subtree_workspace.dev.d_tipchars != nullptr;
    if (!has_workspace_state && workspace.stream == nullptr) {
        return;
    }
    if (workspace.ops.d_ops != nullptr) {
        free_placement_op_buffer(workspace.ops, stream);
    }
    CUDA_CHECK(cudaStreamSynchronize(stream));
    release_subtree_workspace(workspace.subtree_workspace);
    if (workspace.stream != nullptr) {
        CUDA_CHECK(cudaStreamDestroy(workspace.stream));
    }
    workspace = LocalSPREvalWorkspace{};
}

LocalSPRLikelihoodEvaluation evaluate_local_spr_tree_loglikelihood(
    const TopologyRefinementRunContext& ctx,
    const TreeBuildResult& candidate_tree,
    const std::vector<int>& local_smoothing_edges,
    LocalSPREvalWorkspace& workspace,
    cudaStream_t stream)
{
    if (workspace.device_id >= 0) {
        CUDA_CHECK(cudaSetDevice(workspace.device_id));
    }
    HostPacking candidate_host_pack =
        build_local_spr_tree_host_packing(ctx, candidate_tree);
    workspace.tree = candidate_tree;
    workspace.host_pack = std::move(candidate_host_pack);
    load_subtree_workspace(
        workspace.subtree_workspace,
        SubtreeWorkspaceLoadConfig{
            workspace.tree,
            workspace.host_pack,
            ctx.state.eig,
            ctx.rate_weights,
            ctx.rate_multipliers,
            ctx.pi,
            ctx.sites,
            ctx.states,
            ctx.rate_cats,
            ctx.per_rate_scaling,
            nullptr,
            false,
            0,
            0,
        },
        stream,
        "evaluate_local_spr_tree_loglikelihood");
    DeviceTree& eval_dev = workspace.subtree_workspace.dev;
    if (eval_dev.d_query_pmat != nullptr ||
        eval_dev.d_pmat_mid_prox != nullptr ||
        eval_dev.d_pmat_mid_dist != nullptr) {
        ensure_device_tree_current_device(
            eval_dev,
            "evaluate_local_spr_tree_loglikelihood");
        cuda_free_if_set(eval_dev.d_query_pmat);
        cuda_free_if_set(eval_dev.d_pmat_mid_prox);
        cuda_free_if_set(eval_dev.d_pmat_mid_dist);
    }
    workspace.initialized = true;
    if (workspace.tree.nodes.empty() || eval_dev.N == 0) {
        throw std::runtime_error("Local SPR eval produced empty tree/device structures.");
    }
    if (workspace.tree.root_id < 0) {
        throw std::runtime_error("Local SPR eval produced tree with invalid root_id.");
    }
    if (ctx.directional_boundary_node_ids.empty()) {
        UpdateTreeClvs(
            eval_dev,
            workspace.tree,
            workspace.host_pack,
            workspace.ops,
            stream);
    } else {
        install_local_spr_virtual_boundary_messages(
            ctx, eval_dev, stream);
        UpdateTreeClvsPreservingTipClvs(
            eval_dev,
            workspace.tree,
            workspace.host_pack,
            workspace.ops,
            stream);
    }
    double logl = mlipper::likelihood::root::compute_root_loglikelihood(
        eval_dev,
        workspace.tree.root_id,
        nullptr,
        0.0,
        stream);
    bool branch_optimization_accepted = true;
    if (!local_smoothing_edges.empty() &&
        ctx.symmetric_branch_sweeps > 0) {
        mlipper::BranchOptimizationOptions branch_options;
        branch_options.sweeps = ctx.symmetric_branch_sweeps;
        branch_options.newton_iterations = 5;
        branch_options.likelihood_tolerance = 1.0e-3;
        branch_options.update_scheme = mlipper::EdgeUpdateScheme::Jacobi;
        branch_options.clv_retention =
            ctx.directional_boundary_node_ids.empty()
            ? mlipper::ClvRetention::RebuildAll
            : mlipper::ClvRetention::PreserveTips;
        branch_options.acceptance_scope =
            mlipper::AcceptanceScope::LocalSubtree;
        const mlipper::BranchLengthOptimizationResult optimization =
            RunSelectedTreeEdgeJacobiBranchLengthOptimization(
                eval_dev,
                workspace.tree,
                workspace.host_pack,
                workspace.ops,
                ctx.state.eig,
                ctx.rate_multipliers,
                local_smoothing_edges,
                stream,
                branch_options);
        branch_optimization_accepted = optimization.accepted;
        if (branch_optimization_accepted) {
            // The Jacobi optimizer's final sweep already reloads the accepted
            // branch lengths, rebuilds all CLVs, and evaluates this likelihood.
            logl = optimization.log_likelihood_after;
        }
        if (!branch_optimization_accepted) {
            logl = optimization.log_likelihood_before;
        }
    }
    LocalSPRLikelihoodEvaluation result;
    result.total = logl;
    result.optimized_tree = branch_optimization_accepted
        ? workspace.tree
        : candidate_tree;
    return result;
}

static PlacementTuningConfig make_local_spr_placement_tuning(
    mlipper::TopologyMoveType move_type,
    size_t node_count)
{
    PlacementTuningConfig tuning;
    tuning.full_opt_passes = 4;
    if (node_count >
        static_cast<size_t>(std::numeric_limits<int>::max() / 2)) {
        throw std::length_error(
            "Local SPR node count exceeds placement top-k indexing.");
    }
    tuning.export_placement_topk = std::max(
        1, static_cast<int>(node_count) * 2);
#if !defined(MLIPPER_USE_DOUBLE)
    tuning.enable_double_rerank = false;
#endif
    if (move_type == mlipper::TopologyMoveType::NNI) {
        tuning.pendant_branch_min = TOPOLOGY_INTERNAL_BRANCH_LEN_MIN;
        tuning.split_branch_min = TOPOLOGY_INTERNAL_BRANCH_LEN_MIN;
    }
    return tuning;
}

static void append_edge_neighborhood(
    const TreeBuildResult& tree,
    int node_id,
    std::vector<int>& edges)
{
    if (node_id < 0 || node_id >= static_cast<int>(tree.nodes.size())) {
        return;
    }
    edges.push_back(node_id);
    const TreeNode& node = tree.nodes[static_cast<size_t>(node_id)];
    edges.push_back(node.parent);
    edges.push_back(node.left);
    edges.push_back(node.right);
    if (node.parent >= 0) {
        const TreeNode& parent =
            tree.nodes[static_cast<size_t>(node.parent)];
        edges.push_back(parent.left);
        edges.push_back(parent.right);
    }
}

int validate_local_spr_candidates(
    std::vector<LocalSPRCandidateMove> validation_candidates,
    const std::vector<LocalSPRRepairUnit>& repair_units,
    TopologyRefinementRunContext& ctx,
    LocalSPREvalWorkspace& eval_workspace,
    TreeBuildResult& base_tree,
    double& current_logL)
{
    int accepted_candidates = 0;
    for (size_t candidate_idx = 0;
         candidate_idx < validation_candidates.size();
         ++candidate_idx) {
        const LocalSPRCandidateMove candidate =
            validation_candidates[candidate_idx];
        if (candidate.repair_unit_id < 0 ||
            candidate.repair_unit_id >= static_cast<int>(repair_units.size())) {
            continue;
        }

        const LocalSPRRepairUnit& unit =
            repair_units[(size_t)candidate.repair_unit_id];
        if (ctx.move_type != mlipper::TopologyMoveType::NNI &&
            (std::find(
                 ctx.directional_boundary_node_ids.begin(),
                 ctx.directional_boundary_node_ids.end(),
                 candidate.regraft_child_id) !=
                 ctx.directional_boundary_node_ids.end() ||
             std::find(
                 ctx.directional_boundary_node_ids.begin(),
                 ctx.directional_boundary_node_ids.end(),
                 candidate.regraft_parent_id) !=
                 ctx.directional_boundary_node_ids.end())) {
            continue;
        }
        if (!local_spr_candidate_still_legal(
                base_tree,
                unit,
                candidate,
                nullptr)) {
            continue;
        }

        TreeBuildResult candidate_tree = base_tree;
        PruneInfo prune_info;
        if (!prune_subtree_for_spr(
                candidate_tree,
                candidate.prune_root_id,
                prune_info)) {
            continue;
        }

        std::vector<int> current_subtree_nodes;
        collect_subtree_node_ids(
            candidate_tree,
            prune_info.pruned_id,
            current_subtree_nodes);
        const std::vector<char> current_subtree_mask =
            build_local_spr_subtree_mask(candidate_tree, current_subtree_nodes);
        if (candidate.regraft_child_id < 0 ||
            candidate.regraft_child_id >=
                static_cast<int>(candidate_tree.nodes.size())) {
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

        const TreeNode& pruned_node = candidate_tree.nodes[
            static_cast<size_t>(prune_info.pruned_id)];
        const TreeNode& target_node = candidate_tree.nodes[
            static_cast<size_t>(candidate.regraft_child_id)];
        const double pruned_branch_min =
            pruned_node.is_tip ||
            ctx.directional_boundary_tip_node_ids.count(prune_info.pruned_id)
                ? OPT_BRANCH_LEN_MIN
                : TOPOLOGY_INTERNAL_BRANCH_LEN_MIN;
        const double target_branch_min =
            target_node.is_tip ||
            ctx.directional_boundary_tip_node_ids.count(
                candidate.regraft_child_id)
                ? OPT_BRANCH_LEN_MIN
                : TOPOLOGY_INTERNAL_BRANCH_LEN_MIN;
        regraft_subtree_for_spr(
            candidate_tree,
            prune_info,
            candidate.regraft_child_id,
            candidate.pendant_length,
            candidate.proximal_length,
            pruned_branch_min,
            target_branch_min);
        local_spr_assert(
            subtree_fully_inside_mask(
                candidate_tree,
                prune_info.pruned_id,
                unit.envelope_mask),
            "committed subtree escaped repair envelope after regraft");
        std::vector<int> smoothing_edges = candidate.regraft_path_nodes;
        append_edge_neighborhood(
            base_tree,
            candidate.prune_root_id,
            smoothing_edges);
        append_edge_neighborhood(
            base_tree,
            candidate.old_parent_id,
            smoothing_edges);
        append_edge_neighborhood(
            base_tree,
            candidate.regraft_child_id,
            smoothing_edges);
        smoothing_edges.erase(
            std::remove_if(
                smoothing_edges.begin(), smoothing_edges.end(),
                [&](int node_id) {
                    return node_id < 0 ||
                        node_id >= static_cast<int>(base_tree.nodes.size()) ||
                        node_id == base_tree.root_id ||
                        std::find(
                            ctx.directional_boundary_node_ids.begin(),
                            ctx.directional_boundary_node_ids.end(),
                            node_id) !=
                            ctx.directional_boundary_node_ids.end();
                }),
            smoothing_edges.end());
        std::sort(smoothing_edges.begin(), smoothing_edges.end());
        smoothing_edges.erase(
            std::unique(smoothing_edges.begin(), smoothing_edges.end()),
            smoothing_edges.end());

        LocalSPRLikelihoodEvaluation baseline_evaluation;
        if (ctx.symmetric_branch_sweeps > 0) {
            baseline_evaluation = evaluate_local_spr_tree_loglikelihood(
                ctx,
                base_tree,
                smoothing_edges,
                eval_workspace,
                ctx.stream);
        } else {
            baseline_evaluation.total = current_logL;
        }
        LocalSPRLikelihoodEvaluation candidate_evaluation =
            evaluate_local_spr_tree_loglikelihood(
                ctx,
                candidate_tree,
                smoothing_edges,
                eval_workspace,
                ctx.stream);
        const double candidate_logL = candidate_evaluation.total;
        const double comparison_logL = baseline_evaluation.total;
        const double exact_gain = candidate_logL - comparison_logL;
        // Symmetric smoothing evaluates a separately optimized baseline copy.
        // That copy is not the accepted state and can occasionally score below
        // it. A topology move must beat both comparisons; otherwise a locally
        // positive delta can still lower the rebuilt sector likelihood.
        const double accepted_tree_gain = candidate_logL - current_logL;
        const bool positive_move = exact_gain > 0.0 &&
            accepted_tree_gain > 0.0;
        if (positive_move) {
            if (ctx.on_accept_move) {
                ctx.on_accept_move(
                    base_tree,
                    candidate.prune_root_id,
                    candidate.regraft_child_id,
                    candidate.pendant_length,
                    candidate.proximal_length);
            }
            current_logL = candidate_logL;
            if (ctx.symmetric_branch_sweeps > 0) {
                candidate_tree =
                    std::move(candidate_evaluation.optimized_tree);
            }
            base_tree = std::move(candidate_tree);
            ++accepted_candidates;
            // All scores came from the pre-commit topology. Remove candidates
            // whose prune context or destination became stale before trying
            // the next move.
            if (candidate_idx + 1 < validation_candidates.size()) {
                auto keep_begin =
                    validation_candidates.begin() +
                    static_cast<std::ptrdiff_t>(candidate_idx + 1);
                keep_begin = std::remove_if(
                    keep_begin,
                    validation_candidates.end(),
                    [&](const LocalSPRCandidateMove& pending_candidate) {
                        if (pending_candidate.repair_unit_id < 0 ||
                            pending_candidate.repair_unit_id >=
                                static_cast<int>(repair_units.size())) {
                            return true;
                        }
                        // A post-placement re-place sweep may rank several
                        // destinations for the same query tip. Once one is
                        // accepted, every remaining score for that prune root
                        // is stale because it was computed on the pre-move
                        // topology. A later sweep can move the query again
                        // after rebuilding its exact context.
                        if (pending_candidate.prune_root_id ==
                            candidate.prune_root_id) {
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

void rebuild_after_local_spr_round(
    TopologyRefinementRunContext& ctx,
    TreeBuildResult& base_tree)
{
    const size_t boundary_per_node = ctx.state.device.per_node_elems();
    const size_t boundary_scaler_span = ctx.state.device.scaler_elems();
    std::vector<fp_t> boundary_clvs;
    std::vector<unsigned> boundary_scalers;
    if (!ctx.directional_boundary_node_ids.empty()) {
        const size_t boundary_clv_count = mlipper::util::checked_product(
            "local SPR boundary CLVs",
            ctx.directional_boundary_node_ids.size(),
            boundary_per_node);
        const size_t boundary_scaler_count = mlipper::util::checked_product(
            "local SPR boundary scalers",
            ctx.directional_boundary_node_ids.size(),
            boundary_scaler_span);
        boundary_clvs.resize(boundary_clv_count);
        boundary_scalers.resize(boundary_scaler_count);
        for (size_t boundary_idx = 0;
             boundary_idx < ctx.directional_boundary_node_ids.size();
             ++boundary_idx) {
            const int node_id =
                ctx.directional_boundary_node_ids[boundary_idx];
            CUDA_CHECK(cudaMemcpyAsync(
                boundary_clvs.data() + boundary_idx * boundary_per_node,
                ctx.state.device.d_clv_up +
                    static_cast<size_t>(node_id) * boundary_per_node,
                sizeof(fp_t) * boundary_per_node,
                cudaMemcpyDeviceToHost,
                ctx.stream));
            CUDA_CHECK(cudaMemcpyAsync(
                boundary_scalers.data() +
                    boundary_idx * boundary_scaler_span,
                ctx.state.device.d_site_scaler_up +
                    static_cast<size_t>(node_id) * boundary_scaler_span,
                sizeof(unsigned) * boundary_scaler_span,
                cudaMemcpyDeviceToHost,
                ctx.stream));
        }
        CUDA_CHECK(cudaStreamSynchronize(ctx.stream));
    }
    free_placement_op_buffer(ctx.placement_ops, ctx.stream);
    CUDA_CHECK(cudaStreamSynchronize(ctx.stream));

    const PlacementTuningConfig placement_tuning = ctx.placement_ops.tuning;
    ctx.placement_ops = PlacementOpBuffer{};
    ctx.placement_ops.tuning = placement_tuning;
    const PlacementQueryBatch empty_queries;
    HostPacking rebuilt_host_pack =
        build_local_spr_tree_host_packing(ctx, base_tree);
    ctx.state.tree = base_tree;
    ctx.state.host_packing = std::move(rebuilt_host_pack);
    ctx.state.queries = PlacementQueryBatch{};
    reload_device_tree_live_data(
        ctx.state.device,
        ctx.state.tree,
        ctx.state.host_packing,
        &empty_queries,
        ctx.stream);
    if (ctx.directional_boundary_node_ids.empty()) {
        UpdateTreeClvs(
            ctx.state.device,
            ctx.state.tree,
            ctx.state.host_packing,
            ctx.placement_ops,
            ctx.stream);
    } else {
        for (size_t boundary_idx = 0;
             boundary_idx < ctx.directional_boundary_node_ids.size();
             ++boundary_idx) {
            const int node_id =
                ctx.directional_boundary_node_ids[boundary_idx];
            CUDA_CHECK(cudaMemcpyAsync(
                ctx.state.device.d_clv_up +
                    static_cast<size_t>(node_id) * boundary_per_node,
                boundary_clvs.data() + boundary_idx * boundary_per_node,
                sizeof(fp_t) * boundary_per_node,
                cudaMemcpyHostToDevice,
                ctx.stream));
            CUDA_CHECK(cudaMemcpyAsync(
                ctx.state.device.d_site_scaler_up +
                    static_cast<size_t>(node_id) * boundary_scaler_span,
                boundary_scalers.data() +
                    boundary_idx * boundary_scaler_span,
                sizeof(unsigned) * boundary_scaler_span,
                cudaMemcpyHostToDevice,
                ctx.stream));
        }
        UpdateTreeClvsPreservingTipClvs(
            ctx.state.device,
            ctx.state.tree,
            ctx.state.host_packing,
            ctx.placement_ops,
            ctx.stream);
    }

    base_tree = ctx.state.tree;
}

} // namespace

void release_local_spr_session_workspace_set(
    LocalSPRSessionWorkspaceSet& workspace_set,
    cudaStream_t stream)
{
    release_local_spr_eval_workspace(workspace_set.eval_workspace, stream);
    release_local_spr_ranking_workspace(workspace_set.ranking_workspace);
}

// Run iterative SPR/NNI refinement with approximate candidate scoring followed
// by exact validation. Each accepted round rebuilds host/device representations
// and audits the rebuilt likelihood against the validated score; this is the
// synchronization boundary that makes cached node-indexed search state stale.
void run_topology_refinement(TopologyRefinementRunContext& ctx) {
    const char* refinement_name =
        ctx.move_type == mlipper::TopologyMoveType::NNI
            ? "NNI"
            : "Local SPR";
    if (ctx.anchors.empty() &&
        ctx.move_type != mlipper::TopologyMoveType::NNI) {
        std::cout << refinement_name
                  << " skipped: no anchors were provided.\n";
        return;
    }

    int rounds_executed = 0;
    const bool uses_direct_nni_scoring =
        ctx.move_type == mlipper::TopologyMoveType::NNI &&
        ctx.per_rate_scaling;
    const int scoring_lane_count = uses_direct_nni_scoring
        ? 1
        : choose_local_spr_scoring_lane_count(ctx.state.device);
    tbb::task_arena task_arena(scoring_lane_count);
    LocalSPREvalWorkspace& eval_workspace =
        ctx.session_workspaces.eval_workspace;
    LocalSPRRankingWorkspace& ranking_workspace =
        ctx.session_workspaces.ranking_workspace;
    const PlacementTuningConfig local_spr_tuning =
        make_local_spr_placement_tuning(
            ctx.move_type,
            ctx.state.tree.nodes.size());

    const TopologyRefinementSearchContext search_ctx{
        ctx.state,
        ctx.stream,
        ctx.pattern_weights_arg,
        ctx.rate_weights,
        ctx.rate_multipliers,
        ctx.pi,
        ctx.sites,
        ctx.states,
        ctx.rate_cats,
        ctx.per_rate_scaling,
        local_spr_tuning,
        ctx.radius,
        ctx.topk_per_unit,
        ctx.move_type,
        ctx.directional_boundary_node_ids,
        ctx.allowed_nni_central_edge_child_ids,
    };

    ctx.accepted_move_count = 0;
    for (int round = 1; round <= ctx.rounds; ++round) {
        ++rounds_executed;
        double current_logL =
            mlipper::likelihood::root::compute_root_loglikelihood(
                ctx.state.device,
                ctx.state.tree.root_id,
                nullptr,
                0.0,
                ctx.stream);
        const double round_start_logL = current_logL;
        TreeBuildResult base_tree = ctx.state.tree;

        std::vector<LocalSPRRepairUnit> repair_units;
        if (ctx.move_type == mlipper::TopologyMoveType::NNI) {
            repair_units = build_nni_repair_units(
                base_tree,
                ctx.allowed_nni_central_edge_child_ids);
        } else {
            const std::vector<LocalSPRInsertionAnchor> current_anchors =
                refresh_local_spr_anchor_parents(base_tree, ctx.anchors);
            repair_units = build_local_spr_repair_units(
                base_tree,
                current_anchors,
                ctx.cluster_threshold,
                ctx.radius);
        }
        if (repair_units.empty()) {
            std::cout << refinement_name << " round " << round
                      << "/" << ctx.rounds
                      << " skipped: no repair units were generated.\n";
            break;
        }

        LocalSPRSearchSummary search_summary;
        search_summary.unit_count = repair_units.size();
        std::vector<LocalSPRCandidateMove> ranked_candidates =
            rank_local_spr_candidates(
                search_ctx,
                base_tree,
                repair_units,
                ranking_workspace,
                task_arena,
                scoring_lane_count,
                search_summary);

        std::sort(
            ranked_candidates.begin(),
            ranked_candidates.end(),
            local_spr_candidate_better);

        std::vector<LocalSPRCandidateMove> validation_candidates;
        const bool validate_full_ranked_pool =
            ctx.move_type == mlipper::TopologyMoveType::NNI;
        if (validate_full_ranked_pool) {
            validation_candidates = ranked_candidates;
        } else {
            validation_candidates = select_local_spr_candidates(
                ranked_candidates,
                static_cast<int>(base_tree.nodes.size()));
        }
        search_summary.selected_candidates = validation_candidates.size();

        std::cout << refinement_name << " round " << round
                  << "/" << ctx.rounds
                  << " subtree repair units: " << search_summary.unit_count
                  << ", enumerated candidates: "
                  << search_summary.enumerated_candidates
                  << ", retained candidates: "
                  << search_summary.retained_candidates
                  << ", selected: " << search_summary.selected_candidates
                  << " (radius=" << ctx.radius
                  << ", cluster=" << ctx.cluster_threshold
                  << ", topk=" << ctx.topk_per_unit
                  << ", selection="
                  << (validate_full_ranked_pool
                          ? "dynamic-validation"
                          : "one-per-unit")
                  << ", pool=per-unit-topk"
                  << ")\n";

        // Candidate ranking can occupy several full local-tree GPU lanes.
        // Validation needs only the retained host-side candidates and a
        // single exact-likelihood workspace.  Keeping all ranking lanes
        // alive while allocating that workspace caused severe allocator
        // pressure (about 50 seconds for a 500-tip validation on the
        // 200k fixture).  Release them before exact validation; this does
        // not alter candidate scores, ordering, or acceptance semantics.
        release_local_spr_ranking_workspace(ranking_workspace);

        const int accepted_this_round = validate_local_spr_candidates(
            std::move(validation_candidates),
            repair_units,
            ctx,
            eval_workspace,
            base_tree,
            current_logL);
        ctx.accepted_move_count += accepted_this_round;
        if (accepted_this_round <= 0) {
            std::cout << refinement_name << " round " << round
                      << " accepted no candidates; stopping.\n";
            break;
        }
        std::cout << refinement_name << " round " << round
                  << " accepted " << accepted_this_round
                  << " candidate(s).\n";
        rebuild_after_local_spr_round(
            ctx,
            base_tree);
        const double rebuilt_logL =
            mlipper::likelihood::root::compute_root_loglikelihood(
                ctx.state.device,
                ctx.state.tree.root_id,
                nullptr,
                0.0,
                ctx.stream);
        const double rebuild_logL_diff = rebuilt_logL - current_logL;
        std::cout << refinement_name << " likelihood audit:"
                  << " before=" << round_start_logL
                  << " validated=" << current_logL
                  << " rebuilt=" << rebuilt_logL
                  << " accepted_delta="
                  << (rebuilt_logL - round_start_logL)
                  << " diff=" << rebuild_logL_diff
                  << "\n";
        const double rebuild_logL_tolerance =
            1e-10 * std::max(1.0, std::abs(current_logL));
        if (!std::isfinite(rebuilt_logL) ||
            std::abs(rebuild_logL_diff) > rebuild_logL_tolerance) {
            throw std::runtime_error(
                std::string(refinement_name) +
                " rebuilt topology likelihood differs from the validated "
                "accepted likelihood.");
        }
        current_logL = rebuilt_logL;
        ctx.current_tree_newick =
            write_local_spr_context_tree_newick(ctx, base_tree);
    }
    std::cout << refinement_name << " joint refinement done (radius="
              << ctx.radius
              << ", rounds=" << rounds_executed
              << "/" << ctx.rounds << ").\n";
}
