#pragma once

#include <cuda_runtime.h>
#include "util/precision.hpp"

// Build a transition matrix from the eigendecomposition uploaded to GPU.
// Assumes a small state space (for example DNA with 4 states) and supports
// up to 16 states.
__device__ void pmatrix_from_triple_device(
    const fp_t* Vinv,
    const fp_t* V,
    const fp_t* rate_eigenvalues,
    fp_t rate_scale,
    fp_t branch_length,
    fp_t pinv,
    fp_t* out_pmat,
    int state_count);
