/**
 * Derivative kernel declarations.
 */
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
    const int* active_ops);

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
    const int* active_ops);
