#pragma once

#include <cuda_runtime.h>

#include "tree/tree.hpp"

namespace mlipper::likelihood::placement {

// Builds the three-way midpoint partials used to score query placements. The
// traversal operations identify candidate edges; the tree topology is unchanged.
void update_midpoint_partials(
    const DeviceTree& D,
    const NodeOpInfo* d_ops,
    int num_ops,
    cudaStream_t stream = 0);

// Batched placement counterpart of libpll's pll_compute_edge_loglikelihood.
// Writes one log likelihood per operation to borrowed device storage d_out.
void compute_edge_loglikelihoods(
    const DeviceTree& D,
    const NodeOpInfo* d_ops,
    int num_ops,
    const fp_t* d_pendant_pmats,
    const fp_t* d_distal_pmats,
    const fp_t* d_proximal_pmats,
    fp_t* d_out,
    cudaStream_t stream = 0);

} // namespace mlipper::likelihood::placement
