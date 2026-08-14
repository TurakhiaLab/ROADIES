#pragma once
#ifndef PARTIAL_LIKELIHOOD_CUH
#define PARTIAL_LIKELIHOOD_CUH
#include "tree/tree.hpp"
#include <cuda_runtime.h>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include "util/mlipper_util.h"
#if defined(MLIPPER_USE_DOUBLE)
#define SCALE_THRESHOLD_EXPONENT -256
#else
#define SCALE_THRESHOLD_EXPONENT -32
#endif

// Downward specializations for states=4 (ratecat-specific).
template<int RATE_CATS>
__device__ __forceinline__ void compute_downward_inner_inner_ratecat(const DeviceTree& D, const NodeOpInfo& op, unsigned int site);

template<int RATE_CATS>
__device__ __forceinline__ void compute_downward_inner_tip_ratecat(const DeviceTree& D, const NodeOpInfo& op, unsigned int site);

template<int RATE_CATS>
__device__ __forceinline__ void compute_downward_tip_inner_ratecat(const DeviceTree& D, const NodeOpInfo& op, unsigned int site);

// Midpoint helper (states=4) used by placement.
template<int RATE_CATS>
__device__ void compute_midpoint_inner_inner_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site,
    bool proximal_mode = false,
    int op_pmat_idx = 0,
    bool active_thread = true,
    fp_t* shared_target_mat = nullptr,
    fp_t* shared_parent_mat = nullptr);

__device__ void compute_midpoint_inner_inner_states4_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site,
    bool proximal_mode = false,
    int op_pmat_idx = 0,
    bool active_thread = true,
    fp_t* shared_target_mat = nullptr,
    fp_t* shared_parent_mat = nullptr);

// Generic downward helpers.
__device__ void compute_downward_inner_inner_generic(const DeviceTree& D, const NodeOpInfo& op, unsigned int site);
__device__ void compute_downward_inner_tip_generic(const DeviceTree& D, const NodeOpInfo& op, unsigned int site);
__device__ void compute_downward_tip_inner_generic(const DeviceTree& D, const NodeOpInfo& op, unsigned int site);

__global__ void Rtree_Likelihood_Site_Parallel_Upward_Kernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops);

__global__ void InitializeTipClvUpKernel(const DeviceTree D);

// Downward (preorder) site-parallel kernel.
__global__ void Rtree_Likelihood_Site_Parallel_Downward_Kernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops);

__global__ void Rtree_Likelihood_Site_Parallel_Downward_Level_Kernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops);
#endif // PARTIAL_LIKELIHOOD_CUH
