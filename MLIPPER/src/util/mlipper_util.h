#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <cuda_runtime.h>

#include "optimize/optimization_types.hpp"
#include "tree/tree.hpp"

#if defined(MLIPPER_USE_DOUBLE)
constexpr int kClvScaleThresholdExponent = -256;
#else
constexpr int kClvScaleThresholdExponent = -32;
#endif

// Branch length defaults to match epa-ng constants.
constexpr double DEFAULT_BRANCH_LENGTH =
    mlipper::optimization::branch_lengths::kDefaultLength;
constexpr double OPT_BRANCH_LEN_MIN =
    mlipper::optimization::branch_lengths::kPlacementMinimumLength;
// Topology moves must be able to collapse weak internal resolutions without
// also allowing newly placed tip branches to collapse.
constexpr double TOPOLOGY_INTERNAL_BRANCH_LEN_MIN =
    mlipper::optimization::branch_lengths::kTreeMinimumLength;
constexpr double OPT_BRANCH_LEN_MAX =
    mlipper::optimization::branch_lengths::kMaximumLength;
constexpr double OPT_BRANCH_XTOL =
    mlipper::optimization::branch_lengths::kPlacementNewtonTolerance;

template <typename T>
__host__ __device__ inline T scalar_min(T lhs, T rhs) {
    return lhs < rhs ? lhs : rhs;
}

template <typename T>
__host__ __device__ inline T scalar_max(T lhs, T rhs) {
    return lhs > rhs ? lhs : rhs;
}

template <typename T>
__host__ __device__ inline T clamp_scalar(T value, T lower, T upper) {
    if (upper < lower) upper = lower;
    if (value < lower) return lower;
    if (value > upper) return upper;
    return value;
}

#if defined(__CUDACC__)
__device__ __forceinline__ double warp_reduce_sum_double(double value) {
    for (int offset = 16; offset > 0; offset >>= 1) {
        value += __shfl_down_sync(0xffffffff, value, offset);
    }
    return value;
}

__device__ __forceinline__ double block_reduce_sum_double(double value) {
    __shared__ double warp_totals[32];
    const int lane = threadIdx.x & 31;
    const int warp = threadIdx.x >> 5;
    const int active_warps = (blockDim.x + 31) >> 5;
    value = warp_reduce_sum_double(value);
    if (lane == 0) warp_totals[warp] = value;
    __syncthreads();
    if (warp == 0) {
        value = lane < active_warps ? warp_totals[lane] : 0.0;
        return warp_reduce_sum_double(value);
    }
    return 0.0;
}
#endif


__host__ __device__ inline double effective_split_branch_min(
    double total_branch_length,
    double min_branch_length = OPT_BRANCH_LEN_MIN)
{
    if (total_branch_length <= 0.0) return min_branch_length;
    return scalar_min(min_branch_length, 0.5 * total_branch_length);
}

__host__ __device__ inline void normalize_split_branch_lengths(
    double total_branch_length,
    double proposed_proximal_length,
    double min_branch_length,
    double& proximal_length_out,
    double& distal_length_out)
{
    if (total_branch_length <= 0.0) {
        proximal_length_out = min_branch_length;
        distal_length_out = min_branch_length;
        return;
    }

    const double lower_bound = effective_split_branch_min(total_branch_length, min_branch_length);
    const double upper_bound = scalar_max(lower_bound, total_branch_length - lower_bound);
    proximal_length_out = clamp_scalar(proposed_proximal_length, lower_bound, upper_bound);
    distal_length_out = total_branch_length - proximal_length_out;
}

__host__ __device__ inline double sanitize_branch_length(
    double branch_length,
    double min_branch_length = OPT_BRANCH_LEN_MIN,
    double max_branch_length = OPT_BRANCH_LEN_MAX,
    double default_branch_length = DEFAULT_BRANCH_LENGTH)
{
    if (!(branch_length > 0.0)) branch_length = default_branch_length;
    return clamp_scalar(branch_length, min_branch_length, max_branch_length);
}

enum NodeOpType : int {
    OP_TIP_TIP = 0,
    OP_TIP_INNER = 1,
    OP_INNER_INNER = 2,
    OP_DOWN_INNER_INNER = 3,
    OP_DOWN_INNER_TIP = 4,
    OP_DOWN_TIP_INNER = 5,
    OP_DOWN_TIP_TIP   = 6
};

enum ClvPool : uint8_t {
    CLV_POOL_UP = 0,
    CLV_POOL_DOWN = 1
};

enum ClvDir : uint8_t {
    CLV_DIR_UNSET      = 0,
    CLV_DIR_UP         = 1, // child -> parent
    CLV_DIR_DOWN_LEFT  = 2, // parent -> left child
    CLV_DIR_DOWN_RIGHT = 3  // parent -> right child
};

struct NodeOpInfo {
    int parent_id = -1;
    int left_id = -1;
    int right_id = -1;
    int left_tip_index = -1;
    int right_tip_index = -1;
    int op_type = OP_TIP_TIP;
    uint8_t clv_pool = static_cast<uint8_t>(CLV_POOL_UP);
    uint8_t dir_tag = static_cast<uint8_t>(CLV_DIR_UP);
};

// Downward placement operations repurpose right_tip_index as the query row.
// Keeping that encoding here avoids duplicating this subtle dispatch rule in
// every likelihood kernel.
__host__ __device__ inline int node_op_query_idx(const NodeOpInfo& op) {
    return (op.op_type == OP_DOWN_INNER_INNER && op.right_tip_index >= 0)
        ? op.right_tip_index
        : 0;
}

__host__ __device__ inline int node_op_target_id(const NodeOpInfo& op) {
    if (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT)) {
        return op.left_id;
    }
    if (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_RIGHT)) {
        return op.right_id;
    }
    return op.parent_id;
}

inline unsigned int ceil_log2_u32(unsigned int x) {
    if (x <= 1u) return 0u;
    unsigned int v = x - 1u;
    unsigned int r = 0u;
    while (v) { v >>= 1u; ++r; }
    return r;
}

__device__ __forceinline__ size_t per_node_span(const DeviceTree& D) {
    return D.sites * static_cast<size_t>(D.rate_cats) * static_cast<size_t>(D.states);
}

__device__ __forceinline__ int clv_node_slot(
    const DeviceTree& D,
    int node_id)
{
    if (node_id < 0 || node_id >= D.N) return -1;
    if (node_id >= D.capacity_N) return -1;
    return node_id;
}

__device__ __forceinline__ size_t scaler_span(const DeviceTree& D) {
    if (D.per_rate_scaling) {
        return D.sites * static_cast<size_t>(D.rate_cats);
    }
    return D.sites;
}

__device__ __forceinline__ size_t scaler_site_offset(
    const DeviceTree& D,
    size_t site)
{
    if (D.per_rate_scaling) {
        return site * static_cast<size_t>(D.rate_cats);
    }
    return site;
}

// A scaler pointer passed to likelihood kernels already points at one site.
// Select the rate-category entry when per-rate scaling is enabled.
__device__ __forceinline__ unsigned int clv_scaler_rate_slot(
    const DeviceTree& D,
    unsigned int rate_idx)
{
    return D.per_rate_scaling ? rate_idx : 0u;
}

__device__ __forceinline__ unsigned int read_clv_scaler_shift(
    const DeviceTree& D,
    const unsigned int* scaler,
    unsigned int rate_idx)
{
    return scaler ? scaler[clv_scaler_rate_slot(D, rate_idx)] : 0u;
}

__device__ __forceinline__ void write_clv_scaler_shift(
    const DeviceTree& D,
    unsigned int* scaler,
    unsigned int rate_idx,
    unsigned int value)
{
    if (scaler) scaler[clv_scaler_rate_slot(D, rate_idx)] = value;
}

__device__ __forceinline__ void add_clv_scaler_shift(
    const DeviceTree& D,
    unsigned int* scaler,
    unsigned int rate_idx,
    unsigned int shift)
{
    if (scaler && shift) scaler[clv_scaler_rate_slot(D, rate_idx)] += shift;
}

__device__ __forceinline__ unsigned int clv_scale_shift(fp_t maximum)
{
    const fp_t threshold = fp_ldexp(fp_t(1), kClvScaleThresholdExponent);
    return maximum < threshold
        ? static_cast<unsigned int>(-kClvScaleThresholdExponent)
        : 0u;
}

__device__ __forceinline__ void scale_clv_states4_if_needed(
    const DeviceTree& D,
    unsigned int* scaler,
    unsigned int rate_idx,
    fp_t* values,
    fp_t maximum)
{
    const unsigned int shift = clv_scale_shift(maximum);
    if (!shift) return;
    add_clv_scaler_shift(D, scaler, rate_idx, shift);
#pragma unroll
    for (int state = 0; state < 4; ++state) {
        fp_scale_pow2(values[state], static_cast<int>(shift));
    }
}

__device__ __forceinline__ void scale_clv_states4_if_needed(
    const DeviceTree& D,
    unsigned int* scaler,
    unsigned int rate_idx,
    fp_t* values)
{
    scale_clv_states4_if_needed(
        D,
        scaler,
        rate_idx,
        values,
        fp_hmax4(values[0], values[1], values[2], values[3]));
}

__device__ __forceinline__ unsigned int* scaler_ptr_for_node(
    unsigned int* base,
    const DeviceTree& D,
    int node_id,
    size_t site)
{
    if (!base) return nullptr;
    const int slot = clv_node_slot(D, node_id);
    if (slot < 0) return nullptr;
    return base + static_cast<size_t>(slot) * scaler_span(D) + scaler_site_offset(D, site);
}

__device__ __forceinline__ unsigned int* up_scaler_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site)
{
    return scaler_ptr_for_node(D.d_site_scaler_up, D, node_id, site);
}

__device__ __forceinline__ unsigned int* down_scaler_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site)
{
    return scaler_ptr_for_node(D.d_site_scaler_down, D, node_id, site);
}

__device__ __forceinline__ unsigned int* edge_midpoint_scaler_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site)
{
    return scaler_ptr_for_node(D.d_edge_midpoint_scaler, D, node_id, site);
}

__device__ __forceinline__ unsigned int* edge_outside_scaler_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site)
{
    return scaler_ptr_for_node(D.d_edge_outside_scaler, D, node_id, site);
}

template <typename T>
__device__ __forceinline__ T* clv_write_pool_base(
    const DeviceTree& D,
    const NodeOpInfo& op)
{
    return (op.clv_pool == static_cast<uint8_t>(CLV_POOL_DOWN))
        ? reinterpret_cast<T*>(D.d_clv_down)
        : reinterpret_cast<T*>(D.d_clv_up);
}

template <typename T>
__device__ __forceinline__ T* clv_read_pool_base(
    const DeviceTree& D,
    const NodeOpInfo& op)
{
    return (op.clv_pool == static_cast<uint8_t>(CLV_POOL_DOWN))
        ? reinterpret_cast<T*>(D.d_clv_down)
        : reinterpret_cast<T*>(D.d_clv_up);
}

template <typename T>
__device__ __forceinline__ T* clv_pool_ptr_for_node(
    T* base,
    const DeviceTree& D,
    int node_id,
    size_t site_offset = 0)
{
    if (!base) return nullptr;
    const int slot = clv_node_slot(D, node_id);
    if (slot < 0) return nullptr;
    return base + static_cast<size_t>(slot) * per_node_span(D) + site_offset;
}

template <typename T>
__device__ __forceinline__ T* clv_write_ptr_for_node(
    const DeviceTree& D,
    const NodeOpInfo& op,
    int node_id)
{
    T* base = clv_write_pool_base<T>(D, op);
    return clv_pool_ptr_for_node(base, D, node_id);
}

template <typename T>
__device__ __forceinline__ T* clv_read_ptr_for_node(
    const DeviceTree& D,
    const NodeOpInfo& op,
    int node_id)
{
    T* base = clv_read_pool_base<T>(D, op);
    return clv_pool_ptr_for_node(base, D, node_id);
}

template <typename T>
__device__ __forceinline__ T* clv_read_ptr_for_node(
    const DeviceTree& D,
    int node_id)
{
    return clv_pool_ptr_for_node(reinterpret_cast<T*>(D.d_clv_up), D, node_id);
}

template <typename T>
__device__ __forceinline__ T* up_clv_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site_offset = 0)
{
    return clv_pool_ptr_for_node(reinterpret_cast<T*>(D.d_clv_up), D, node_id, site_offset);
}

template <typename T>
__device__ __forceinline__ T* down_clv_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site_offset = 0)
{
    return clv_pool_ptr_for_node(reinterpret_cast<T*>(D.d_clv_down), D, node_id, site_offset);
}

template <typename T>
__device__ __forceinline__ T* edge_midpoint_clv_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site_offset = 0)
{
    return clv_pool_ptr_for_node(reinterpret_cast<T*>(D.d_edge_midpoint_clv), D, node_id, site_offset);
}

template <typename T>
__device__ __forceinline__ T* edge_outside_clv_ptr(
    const DeviceTree& D,
    int node_id,
    size_t site_offset = 0)
{
    return clv_pool_ptr_for_node(reinterpret_cast<T*>(D.d_edge_outside_clv), D, node_id, site_offset);
}

__device__ __forceinline__ unsigned int* site_scaler_ptr_base(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site,
    unsigned int rate_cats)
{
    (void)rate_cats;
    if (op.clv_pool == static_cast<uint8_t>(CLV_POOL_DOWN)) {
        return down_scaler_ptr(D, op.parent_id, site);
    }
    return up_scaler_ptr(D, op.parent_id, site);
}
