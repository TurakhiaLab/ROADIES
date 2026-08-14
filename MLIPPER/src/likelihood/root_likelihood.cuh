#pragma once
#ifndef TREE_GENERATION_ROOT_LIKELIHOOD_CUH
#define TREE_GENERATION_ROOT_LIKELIHOOD_CUH

#include <cuda_runtime.h>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <vector>
#include "tree/tree.hpp"
#include "partial_likelihood.cuh"

namespace root_likelihood {

double compute_root_loglikelihood_total(
    const DeviceTree& D,
    int root_id,
    const unsigned* d_pattern_w,
    const int* d_invar_indices,
    double invar_proportion,
    cudaStream_t stream = 0);

// Compute placement root log-likelihood per placement op into device buffer.
void Placement_Root_Loglk(
    const DeviceTree& D,
    const NodeOpInfo* d_ops,
    const int* d_op_indices,
    int num_ops,
    const fp_t* d_pendant_pmats, // [num_ops * rate_cats * states * states]
    const fp_t* d_distal_pmats,  // [N * rate_cats * states * states]
    const fp_t* d_proximal_pmats,// [N * rate_cats * states * states]
    fp_t* d_out,                  // [num_ops]
    cudaStream_t stream = 0);

} // namespace root_likelihood

#endif // TREE_GENERATION_ROOT_LIKELIHOOD_CUH
