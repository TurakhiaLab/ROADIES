#pragma once

#include <cuda_runtime.h>

#include <vector>

#include "optimize/optimization_types.hpp"
#include "tree/tree.hpp"

struct RawPlacementResult {
    int target_id = -1;
    double loglikelihood = 0.0;
    double distal_length = 0.0;
    double pendant_length = 0.0;
    struct RankedPlacement {
        int target_id = -1;
        double loglikelihood = 0.0;
        double distal_length = 0.0;
        double pendant_length = 0.0;
        double like_weight_ratio = 0.0;
    };
    std::vector<RankedPlacement> top_placements;
};

struct PlacementScratchOverride {
    fp_t* d_sumtable = nullptr;
    fp_t* d_likelihoods = nullptr;
    fp_t* d_query_pmat = nullptr;
    fp_t* d_prev_loglk = nullptr;
    int* d_active_ops = nullptr;
    fp_t* d_pmat_mid_prox = nullptr;
    fp_t* d_pmat_mid_dist = nullptr;
    fp_t* d_edge_midpoint_clv = nullptr;
    unsigned* d_edge_midpoint_scaler = nullptr;
    fp_t* d_new_pendant_length = nullptr;
    fp_t* d_new_proximal_length = nullptr;
    fp_t* d_prev_pendant_length = nullptr;
    fp_t* d_prev_proximal_length = nullptr;
    size_t sumtable_capacity_ops = 0;
    size_t likelihood_capacity_ops = 0;
    size_t query_pmat_capacity_ops = 0;
    size_t ranking_state_capacity_ops = 0;
    size_t midpoint_pmat_capacity_nodes = 0;
    size_t midpoint_clv_capacity_nodes = 0;
    size_t branch_length_capacity_nodes = 0;
};

// Scores one query against d_ops without changing the tree. All override
// pointers are borrowed for the duration of the call; their matching capacity
// fields are validated before kernels launch. d_ops use child-endpoint target
// IDs and must describe the same resident topology as D. smoothing bounds the
// branch-coordinate iterations; returned target IDs remain in that topology.
RawPlacementResult EvaluatePlacementCandidates(
    const DeviceTree& D,
    const NodeOpInfo* d_ops,
    int num_ops,
    int smoothing,
    cudaStream_t stream,
    bool enable_local_child_refine,
    const PlacementTuningConfig& tuning = PlacementTuningConfig{},
    const PlacementScratchOverride* scratch_override = nullptr);

// Resident full-tree Gauss-Seidel branch optimization. Each edge update is
// committed before the next edge is scored; no site-batch workspace is used.
mlipper::BranchLengthOptimizationResult
RunAcceptedFullTreeSequentialBranchLengthOptimization(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    const PlacementQueryBatch* queries,
    cudaStream_t stream,
    const mlipper::BranchOptimizationOptions& options);

// Jacobi optimization proposes every selected edge from one shared CLV state,
// then accepts or rejects the edge group according to options.
mlipper::BranchLengthOptimizationResult
RunSelectedTreeEdgeJacobiBranchLengthOptimization(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    const std::vector<int>& edge_child_ids,
    cudaStream_t stream,
    const mlipper::BranchOptimizationOptions& options);
