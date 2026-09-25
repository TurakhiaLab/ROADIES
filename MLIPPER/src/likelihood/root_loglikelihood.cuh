#pragma once

#include <cuda_runtime.h>

#include "tree/tree.hpp"

namespace mlipper::likelihood::root {

// Computes the likelihood at a selected upward CLV, matching the role of
// libpll's pll_compute_root_loglikelihood. This call synchronizes stream before
// returning because the reduced likelihood is copied back to the host.
double compute_root_loglikelihood(
    const DeviceTree& D,
    int root_id,
    const int* d_invar_indices,
    double invar_proportion,
    cudaStream_t stream = 0);

} // namespace mlipper::likelihood::root
