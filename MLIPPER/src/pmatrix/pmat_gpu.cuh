#pragma once

#include <cuda_runtime.h>

#include "util/precision.hpp"

// Build a transition matrix from the eigendecomposition uploaded to GPU.
// Assumes a small state space (for example DNA with 4 states) and supports
// up to 16 states.
__device__ void pmatrix_from_triple_gpu(
    const fp_t* Vinv,
    const fp_t* V,
    const fp_t* rate_eigenvalues,
    fp_t rate_scale,
    fp_t branch_length,
    fp_t pinv,
    fp_t* out_pmat,
    int state_count);

// Rebuild [node][rate][row][column] PMATs from device-resident branch lengths
// and eigendecomposition data. Invalid host-side arguments throw.
void build_all_branch_pmats_gpu(
    int node_count,
    int state_count,
    int rate_categories,
    const fp_t* d_branch_lengths,
    const fp_t* d_eigenvectors,
    const fp_t* d_inverse_eigenvectors,
    const fp_t* d_rate_eigenvalues,
    fp_t* d_pmats,
    cudaStream_t stream = 0);

// Rebuild one branch PMAT after a sequential coordinate update. The branch
// length remains device-resident, so a BLO traversal does not synchronize a
// scalar through the host between adjacent edges. Invalid arguments throw.
void build_single_branch_pmat_gpu(
    int node_id,
    int state_count,
    int rate_categories,
    const fp_t* d_branch_lengths,
    const fp_t* d_eigenvectors,
    const fp_t* d_inverse_eigenvectors,
    const fp_t* d_rate_eigenvalues,
    fp_t* d_pmats,
    cudaStream_t stream = 0);
