#pragma once

#include <cuda_runtime.h>

#include <cstddef>

#include "tree/tree.hpp"

// Pendant-side derivative kernel. It builds midpoint state and derivative
// sumtable rows directly inside the kernel.
__global__ void LikelihoodDerivativePendantKernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int op_idx,
    const int* op_indices,
    const int*    __restrict__ invariant_site,
    const fp_t* __restrict__ invar_proportion,
    fp_t invar_scalar,
    fp_t* __restrict__ sumtable,
    const unsigned* __restrict__ pattern_weights,
    int max_iter,
    fp_t* new_branch_length,
    size_t sumtable_stride,
    const fp_t* prev_branch_lengths,
    const int* active_ops,
    double branch_min);

// Proximal-side derivative kernel. It builds midpoint state and derivative
// sumtable rows directly inside the kernel.
__global__ void LikelihoodDerivativeProximalKernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int op_idx,
    const int* op_indices,
    const int*    __restrict__ invariant_site,
    const fp_t* __restrict__ invar_proportion,
    fp_t invar_scalar,
    fp_t* __restrict__ sumtable,
    const unsigned* __restrict__ pattern_weights,
    int max_iter,
    fp_t* new_branch_length,
    size_t sumtable_stride,
    const fp_t* prev_branch_lengths,
    const int* active_ops,
    double branch_min);

// Optimize one existing tree edge from the two directional CLVs surrounding
// it. The edge sumtable is built once and reused throughout the local Newton
// loop, matching the RAxML/corax branch-coordinate contract.
void OptimizeSingleTreeEdgeFromCurrentClvs(
    const DeviceTree& D,
    int target_id,
    fp_t* d_sumtable,
    double* d_partial_gradient,
    double* d_partial_hessian,
    double* d_newton_state,
    int* d_newton_failure,
    int max_iter,
    cudaStream_t stream = 0);

// High-occupancy DNA+G4 path used by sequential tree-edge optimization. A
// half-warp cooperates on each site; other model shapes use the generic path.
void OptimizeSingleTreeEdgeFromCurrentClvsWarpSite(
    const DeviceTree& D,
    int target_id,
    fp_t* d_sumtable,
    double* d_partial_gradient,
    double* d_partial_hessian,
    double* d_newton_state,
    int* d_newton_failure,
    int max_iter,
    cudaStream_t stream = 0);
