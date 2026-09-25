#pragma once

#include <cuda_runtime.h>

#include "tree/tree.hpp"

namespace mlipper::likelihood::partials {

// Serial traversal form: every site walks the ordered postorder operation list.
// Use this when operations may depend on earlier operations in the same launch.
__global__ void UpdatePartialsUpwardKernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops);

// Level form: blockIdx.y selects an operation. All operations supplied in one
// launch must be independent members of the same traversal level.
__global__ void UpdatePartialsUpwardLevelKernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops);

__global__ void InitializeTipPartialsKernel(const DeviceTree D);

// Downward equivalents follow preorder dependencies and the same serial/level
// scheduling contract as their upward counterparts.
__global__ void UpdatePartialsDownwardKernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops);

__global__ void UpdatePartialsDownwardLevelKernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops);

// DNA+G4 single-operation kernels. Sixteen lanes cooperate on each site,
// exposing one lane per rate/output-state component.
__global__ void BuildTreeEdgeOutsideWarpSiteKernel(
    const DeviceTree D,
    const NodeOpInfo* op);
__global__ void UpdateTreeUpwardWarpSiteKernel(
    const DeviceTree D,
    const NodeOpInfo* op);
__global__ void RefreshTreeChildDownWarpSiteKernel(
    const DeviceTree D,
    int target_id);

} // namespace mlipper::likelihood::partials
