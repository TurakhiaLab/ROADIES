#pragma once
#ifndef MLIPPER_PLACEMENT_CUH
#define MLIPPER_PLACEMENT_CUH

#include <limits>
#include <vector>
#include <cuda_runtime.h>
#include "tree/tree.hpp"
#include "likelihood/partial_likelihood.cuh"
#include "likelihood/root_likelihood.cuh"

// Host entry points for evaluating placement likelihoods while keeping ops on device.
void InsertLikelihoodEvaluationKernel(
    const DeviceTree& D,
    PlacementQueryBatch& Q,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    const NodeOpInfo* d_ops,
    int op_idx,
    fp_t* d_likelihoods,
    void* d_temp_reduce,
    size_t temp_bytes_reduce,
    cudaStream_t stream);

struct PlacementResult {
    int target_id = -1;
    double loglikelihood = 0.0;
    // This stores the jplace distal coordinate for the chosen edge.
    double proximal_length = 0.0;
    double pendant_length = 0.0;
    struct RankedPlacement {
        int target_id = -1;
        double loglikelihood = 0.0;
        // Same convention as PlacementResult::proximal_length above.
        double proximal_length = 0.0;
        double pendant_length = 0.0;
        double like_weight_ratio = 0.0;
    };
    std::vector<RankedPlacement> top_placements;
};

PlacementResult PlacementEvaluationKernel(
    const DeviceTree& D,
    const NodeOpInfo* d_ops,
    int num_ops,
    int smoothing,
    cudaStream_t stream,
    bool enable_local_child_refine);

#endif // TREE_GENERATION_TREE_PLACEMENT_CUH
