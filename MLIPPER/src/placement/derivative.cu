#include <cmath>
#include <chrono>
#include <cstdio>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <cooperative_groups.h>
#include <cuda_runtime.h>

#include "derivative.cuh"
#include "gpu/device_buffer.hpp"
#include "likelihood/root_loglikelihood.cuh"
#include "optimize/optimization_types.hpp"
#include "placement.cuh"
#include "pmatrix/pmat_gpu.cuh"
#include "tree/tree.hpp"
#include "tree/tree_topology_utils.hpp"
#include "util/mlipper_util.h"

constexpr int kMaxRateCats = 8;
using mlipper::optimization::branch_lengths::kMaximumLength;
using mlipper::optimization::branch_lengths::kNewtonTolerance;
using mlipper::optimization::branch_lengths::kTreeMinimumLength;

// Branch optimization builds a sumtable from the two likelihood messages on
// an edge. That table is independent of the trial branch length, so Newton
// iterations vary only exp(lambda*t) and its first two derivatives. Placement
// kernels apply the same mechanism separately to pendant and proximal
// coordinates; full-tree kernels update existing child edges.

static __device__ __forceinline__ unsigned int threshold_scale_shift(fp_t max_val)
{
    const fp_t scale_threshold = fp_ldexp(fp_t(1), kClvScaleThresholdExponent);
    return (max_val < scale_threshold)
        ? static_cast<unsigned int>(-kClvScaleThresholdExponent)
        : 0u;
}

// Each eigenvalue occupies four adjacent entries: exp(lambda*t), its first and
// second branch-length derivatives, and one reserved padding slot.
static __device__ __forceinline__ void build_diagtable_for_branch(
    const DeviceTree& D,
    fp_t branch_length,
    fp_t* diag_out)
{
    if (!diag_out || !D.d_lambdas) return;
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t total = rate_count * state_count;
    const size_t idx = threadIdx.x;
    if (idx >= total) return;
    const size_t rc = idx / state_count;
    const size_t st = idx - rc * state_count;
    const fp_t lambda = D.d_lambdas[rc * state_count + st];
    const fp_t e = fp_exp(lambda * branch_length);
    const size_t base = idx * 4;
    diag_out[base + 0] = e;
    diag_out[base + 1] = lambda * e;
    diag_out[base + 2] = lambda * lambda * e;
    diag_out[base + 3] = fp_t(0);
}

template<int RATE_CATS>
static __device__ __forceinline__ void build_diagtable_states4(
    const DeviceTree& D,
    fp_t branch_length,
    fp_t* diag_out)
{
    if (!diag_out || !D.d_lambdas) return;
    const size_t idx = threadIdx.x;
    if (idx >= static_cast<size_t>(RATE_CATS * 4)) return;
    const fp_t lambda = D.d_lambdas[idx];
    const fp_t e = fp_exp(lambda * branch_length);
    const size_t base = idx * 4;
    diag_out[base + 0] = e;
    diag_out[base + 1] = lambda * e;
    diag_out[base + 2] = lambda * lambda * e;
    diag_out[base + 3] = fp_t(0);
}

static __device__ __forceinline__ bool build_diagtable_states4_dispatch(
    const DeviceTree& D,
    fp_t branch_length,
    fp_t* diag_out)
{
    switch (D.rate_cats) {
        case 1:
            build_diagtable_states4<1>(D, branch_length, diag_out);
            return true;
        case 4:
            build_diagtable_states4<4>(D, branch_length, diag_out);
            return true;
        case 8:
            build_diagtable_states4<8>(D, branch_length, diag_out);
            return true;
        default:
            return false;
    }
}

static __device__ __forceinline__
unsigned int scaler_shift_at_site(
    const DeviceTree& D,
    const unsigned* __restrict__ scaler_base,
    size_t site_idx,
    size_t rate_idx);

static __device__ __forceinline__
unsigned int scaler_shift_from_site_ptr(
    const DeviceTree& D,
    const unsigned* __restrict__ scaler_site_ptr,
    size_t rate_idx);

// A block stages both half-branch PMATs before any thread evaluates sites.
// Every caller must participate because these helpers contain a block barrier.
template<int RATE_CATS>
__device__ __forceinline__ void load_midpoint_pmat_pair(
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat,
    const fp_t* target_mat,
    const fp_t* parent_mat)
{
    const int total_mat_elems = RATE_CATS * 16;
    for (int idx = threadIdx.x; idx < total_mat_elems; idx += blockDim.x) {
        shared_target_mat[idx] = target_mat[idx];
        shared_parent_mat[idx] = parent_mat[idx];
    }
    __syncthreads();
}

static __device__ __forceinline__ void load_midpoint_pmat_pair_generic(
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat,
    const fp_t* target_mat,
    const fp_t* parent_mat,
    int rate_cats)
{
    const int total_mat_elems = rate_cats * 16;
    for (int idx = threadIdx.x; idx < total_mat_elems; idx += blockDim.x) {
        shared_target_mat[idx] = target_mat[idx];
        shared_parent_mat[idx] = parent_mat[idx];
    }
    __syncthreads();
}

static __device__ __forceinline__
bool load_midpoint_pmat_pair_dispatch(
    const DeviceTree& D,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat,
    const fp_t* target_mat,
    const fp_t* parent_mat)
{
    switch (D.rate_cats) {
        case 1:
            load_midpoint_pmat_pair<1>(
                shared_target_mat, shared_parent_mat, target_mat, parent_mat);
            return true;
        case 4:
            load_midpoint_pmat_pair<4>(
                shared_target_mat, shared_parent_mat, target_mat, parent_mat);
            return true;
        case 8:
            load_midpoint_pmat_pair<8>(
                shared_target_mat, shared_parent_mat, target_mat, parent_mat);
            return true;
        default:
            if (D.rate_cats <= 0 || D.rate_cats > kMaxRateCats) return false;
            load_midpoint_pmat_pair_generic(
                shared_target_mat, shared_parent_mat, target_mat, parent_mat, D.rate_cats);
            return true;
    }
}

// The explicit row/column variants make the eigenbasis orientation visible at
// call sites; substituting one for the other changes the derivative projection.
static __device__ __forceinline__
fp4_t matvec4_rows(const fp_t* mat_rows, const fp4_t& vec)
{
    const fp4_t* rows = reinterpret_cast<const fp4_t*>(mat_rows);
    return make_fp4(
        fp_dot4(rows[0], vec),
        fp_dot4(rows[1], vec),
        fp_dot4(rows[2], vec),
        fp_dot4(rows[3], vec));
}

static __device__ __forceinline__
fp4_t matvec4_cols(const fp_t* mat_cols, const fp4_t& vec)
{
    return make_fp4(
        fp_fma(mat_cols[0], vec.x,
               fp_fma(mat_cols[4], vec.y,
                      fp_fma(mat_cols[8], vec.z, mat_cols[12] * vec.w))),
        fp_fma(mat_cols[1], vec.x,
               fp_fma(mat_cols[5], vec.y,
                      fp_fma(mat_cols[9], vec.z, mat_cols[13] * vec.w))),
        fp_fma(mat_cols[2], vec.x,
               fp_fma(mat_cols[6], vec.y,
                      fp_fma(mat_cols[10], vec.z, mat_cols[14] * vec.w))),
        fp_fma(mat_cols[3], vec.x,
               fp_fma(mat_cols[7], vec.y,
                      fp_fma(mat_cols[11], vec.z, mat_cols[15] * vec.w))));
}

// Build one pendant-side midpoint vector for a single site and rate category.
template<int RATE_CATS>
__device__ __forceinline__ fp4_t build_pendant_midpoint_site_rate(
    const DeviceTree& D,
    int target_id,
    size_t site_idx,
    size_t rate_idx,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat,
    unsigned int* total_shift_out,
    unsigned int* down_shift_out,
    unsigned int* up_shift_out,
    fp4_t* pbase_out,
    fp4_t* pup_out)
{
    const size_t rate_count = static_cast<size_t>(RATE_CATS);
    const size_t site_span = rate_count * 4;
    const size_t site_base = site_idx * site_span;
    const size_t rate_base = rate_idx * 4;

    const fp_t* edge_outside = edge_outside_clv_ptr<const fp_t>(D, target_id, site_base);
    const fp_t* target_up = up_clv_ptr<const fp_t>(D, target_id, site_base);
    if (!edge_outside || !target_up) return make_fp4(fp_t(0), fp_t(0), fp_t(0), fp_t(0));
    const unsigned* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site_idx);
    const unsigned* up_scaler = up_scaler_ptr(D, target_id, site_idx);

    const unsigned int down_shift = scaler_shift_from_site_ptr(D, edge_outside_scaler, rate_idx);
    const unsigned int up_shift = scaler_shift_from_site_ptr(D, up_scaler, rate_idx);
    unsigned int inherited_shift = down_shift + up_shift;

    const size_t mat_base = rate_idx * 16;
    const fp_t* Mtarget = shared_target_mat + mat_base;
    const fp_t* Mparent = shared_parent_mat + mat_base;
    const fp4_t Pup = reinterpret_cast<const fp4_t*>(target_up + rate_base)[0];
    const fp4_t Pbase = reinterpret_cast<const fp4_t*>(edge_outside + rate_base)[0];
    if (pbase_out) *pbase_out = Pbase;
    if (pup_out) *pup_out = Pup;
    if (down_shift_out) *down_shift_out = down_shift;
    if (up_shift_out) *up_shift_out = up_shift;

    fp4_t midpoint = matvec4_rows(Mparent, Pbase);
    const fp4_t target_proj = matvec4_rows(Mtarget, Pup);
    midpoint.x *= target_proj.x;
    midpoint.y *= target_proj.y;
    midpoint.z *= target_proj.z;
    midpoint.w *= target_proj.w;

    fp_t row_max = fp_hmax4(midpoint.x, midpoint.y, midpoint.z, midpoint.w);
    unsigned int total_shift = inherited_shift;
    {
        const unsigned int shift = threshold_scale_shift(row_max);
        if (shift) {
            total_shift += shift;
            midpoint.x = fp_ldexp(midpoint.x, shift);
            midpoint.y = fp_ldexp(midpoint.y, shift);
            midpoint.z = fp_ldexp(midpoint.z, shift);
            midpoint.w = fp_ldexp(midpoint.w, shift);
        }
    }

    if (total_shift_out) *total_shift_out = total_shift;
    return midpoint;
}

__device__ __forceinline__ fp4_t build_pendant_midpoint_site_rate_generic(
    const DeviceTree& D,
    int target_id,
    size_t site_idx,
    size_t rate_idx,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat,
    unsigned int* total_shift_out,
    unsigned int* down_shift_out,
    unsigned int* up_shift_out,
    fp4_t* pbase_out,
    fp4_t* pup_out)
{
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t site_span = rate_count * 4;
    const size_t site_base = site_idx * site_span;
    const size_t rate_base = rate_idx * 4;

    const fp_t* edge_outside = edge_outside_clv_ptr<const fp_t>(D, target_id, site_base);
    const fp_t* target_up = up_clv_ptr<const fp_t>(D, target_id, site_base);
    if (!edge_outside || !target_up) return make_fp4(fp_t(0), fp_t(0), fp_t(0), fp_t(0));
    const unsigned* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site_idx);
    const unsigned* up_scaler = up_scaler_ptr(D, target_id, site_idx);

    const unsigned int down_shift = scaler_shift_from_site_ptr(D, edge_outside_scaler, rate_idx);
    const unsigned int up_shift = scaler_shift_from_site_ptr(D, up_scaler, rate_idx);
    unsigned int inherited_shift = down_shift + up_shift;

    const size_t mat_base = rate_idx * 16;
    const fp_t* Mtarget = shared_target_mat + mat_base;
    const fp_t* Mparent = shared_parent_mat + mat_base;
    const fp4_t Pup = reinterpret_cast<const fp4_t*>(target_up + rate_base)[0];
    const fp4_t Pbase = reinterpret_cast<const fp4_t*>(edge_outside + rate_base)[0];
    if (pbase_out) *pbase_out = Pbase;
    if (pup_out) *pup_out = Pup;
    if (down_shift_out) *down_shift_out = down_shift;
    if (up_shift_out) *up_shift_out = up_shift;

    fp4_t midpoint = matvec4_rows(Mparent, Pbase);
    const fp4_t target_proj = matvec4_rows(Mtarget, Pup);
    midpoint.x *= target_proj.x;
    midpoint.y *= target_proj.y;
    midpoint.z *= target_proj.z;
    midpoint.w *= target_proj.w;

    fp_t row_max = fp_hmax4(midpoint.x, midpoint.y, midpoint.z, midpoint.w);
    unsigned int total_shift = inherited_shift;
    {
        const unsigned int shift = threshold_scale_shift(row_max);
        if (shift) {
            total_shift += shift;
            midpoint.x = fp_ldexp(midpoint.x, shift);
            midpoint.y = fp_ldexp(midpoint.y, shift);
            midpoint.z = fp_ldexp(midpoint.z, shift);
            midpoint.w = fp_ldexp(midpoint.w, shift);
        }
    }

    if (total_shift_out) *total_shift_out = total_shift;
    return midpoint;
}

// Build proximal-side midpoint vectors for a single site across all rate categories.
template<int RATE_CATS>
__device__ __forceinline__ void build_proximal_midpoint_site(
    const DeviceTree& D,
    int target_id,
    size_t site_idx,
    const fp_t* __restrict__ query_clv_base,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat,
    fp4_t* midpoint_rows,
    unsigned int* midpoint_shifts,
    fp4_t* edge_outside_rows,
    fp4_t* query_rows,
    fp4_t* parent_proj_rows,
    fp4_t* target_proj_rows)
{
    const size_t rate_count = static_cast<size_t>(RATE_CATS);
    const size_t site_span = rate_count * 4;
    const size_t site_base = site_idx * site_span;
    const fp_t* edge_outside = edge_outside_clv_ptr<const fp_t>(D, target_id, site_base);
    if (!edge_outside) return;
    const fp_t* query_clv = query_clv_base + site_base;
    const unsigned* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site_idx);

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        unsigned int inherited_shift = scaler_shift_from_site_ptr(D, edge_outside_scaler, r);

        const size_t rate_offset = static_cast<size_t>(r);
        const size_t rate_base = rate_offset * 4;
        const size_t mat_base = rate_offset * 16;
        const fp_t* Mtarget = shared_target_mat + mat_base;
        const fp_t* Mparent = shared_parent_mat + mat_base;
        const fp4_t Pup = reinterpret_cast<const fp4_t*>(query_clv + rate_base)[0];
        const fp4_t Pbase = reinterpret_cast<const fp4_t*>(edge_outside + rate_base)[0];
        const fp4_t parent_proj = matvec4_rows(Mparent, Pbase);
        const fp4_t target_proj = matvec4_rows(Mtarget, Pup);

        if (edge_outside_rows) edge_outside_rows[r] = Pbase;
        if (query_rows) query_rows[r] = Pup;
        if (parent_proj_rows) parent_proj_rows[r] = parent_proj;
        if (target_proj_rows) target_proj_rows[r] = target_proj;

        fp_t p0 = parent_proj.x * target_proj.x;
        fp_t p1 = parent_proj.y * target_proj.y;
        fp_t p2 = parent_proj.z * target_proj.z;
        fp_t p3 = parent_proj.w * target_proj.w;

        fp_t row_max = fp_hmax4(p0, p1, p2, p3);
        unsigned int total_shift = inherited_shift;
        {
            const unsigned int shift = threshold_scale_shift(row_max);
            if (shift) {
                total_shift += shift;
                p0 = fp_ldexp(p0, shift);
                p1 = fp_ldexp(p1, shift);
                p2 = fp_ldexp(p2, shift);
                p3 = fp_ldexp(p3, shift);
            }
        }

        midpoint_rows[r] = make_fp4(p0, p1, p2, p3);
        midpoint_shifts[r] = total_shift;
    }
}

__device__ __forceinline__ void build_proximal_midpoint_site_generic(
    const DeviceTree& D,
    int target_id,
    size_t site_idx,
    const fp_t* __restrict__ query_clv_base,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat,
    fp4_t* midpoint_rows,
    unsigned int* midpoint_shifts,
    fp4_t* edge_outside_rows,
    fp4_t* query_rows,
    fp4_t* parent_proj_rows,
    fp4_t* target_proj_rows)
{
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t site_span = rate_count * 4;
    const size_t site_base = site_idx * site_span;
    const fp_t* edge_outside = edge_outside_clv_ptr<const fp_t>(D, target_id, site_base);
    if (!edge_outside) return;
    const fp_t* query_clv = query_clv_base + site_base;
    const unsigned* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site_idx);

    for (int r = 0; r < D.rate_cats; ++r) {
        unsigned int inherited_shift = scaler_shift_from_site_ptr(D, edge_outside_scaler, r);

        const size_t rate_offset = static_cast<size_t>(r);
        const size_t rate_base = rate_offset * 4;
        const size_t mat_base = rate_offset * 16;
        const fp_t* Mtarget = shared_target_mat + mat_base;
        const fp_t* Mparent = shared_parent_mat + mat_base;
        const fp4_t Pup = reinterpret_cast<const fp4_t*>(query_clv + rate_base)[0];
        const fp4_t Pbase = reinterpret_cast<const fp4_t*>(edge_outside + rate_base)[0];
        const fp4_t parent_proj = matvec4_rows(Mparent, Pbase);
        const fp4_t target_proj = matvec4_rows(Mtarget, Pup);

        if (edge_outside_rows) edge_outside_rows[r] = Pbase;
        if (query_rows) query_rows[r] = Pup;
        if (parent_proj_rows) parent_proj_rows[r] = parent_proj;
        if (target_proj_rows) target_proj_rows[r] = target_proj;

        fp_t p0 = parent_proj.x * target_proj.x;
        fp_t p1 = parent_proj.y * target_proj.y;
        fp_t p2 = parent_proj.z * target_proj.z;
        fp_t p3 = parent_proj.w * target_proj.w;

        fp_t row_max = fp_hmax4(p0, p1, p2, p3);
        unsigned int total_shift = inherited_shift;
        {
            const unsigned int shift = threshold_scale_shift(row_max);
            if (shift) {
                total_shift += shift;
                p0 = fp_ldexp(p0, shift);
                p1 = fp_ldexp(p1, shift);
                p2 = fp_ldexp(p2, shift);
                p3 = fp_ldexp(p3, shift);
            }
        }

        midpoint_rows[r] = make_fp4(p0, p1, p2, p3);
        midpoint_shifts[r] = total_shift;
    }
}

// Build one site's derivative sumtable rows for the pendant-side branch update.
template<int RATE_CATS>
__device__ __forceinline__ void update_pendant_sumtable_site(
    DeviceTree D,
    int target_id,
    const fp_t* __restrict__ left_clv_base,
    size_t site_idx,
    fp_t* sumtable,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat)
{
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t site_span = rate_count * state_count;
    const size_t site_base = site_idx * site_span;
    const fp_t* left_clv = left_clv_base + site_base;
    fp_t* sumtable_ptr = sumtable + site_base;

    unsigned int midpoint_shifts[RATE_CATS];
    unsigned int active_rate_mask = 0u;
    unsigned int site_min_shift = 0u;
    bool have_signal = false;

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        const size_t rate_offset = static_cast<size_t>(r);
        const size_t rate_base = rate_offset * 4;
        const fp4_t Qclv = reinterpret_cast<const fp4_t*>(left_clv + rate_base)[0];
        const fp4_t Pclv = build_pendant_midpoint_site_rate<RATE_CATS>(
            D,
            target_id,
            site_idx,
            r,
            shared_target_mat,
            shared_parent_mat,
            &midpoint_shifts[r],
            nullptr,
            nullptr,
            nullptr,
            nullptr);
        fp_t* sumtable_row = sumtable_ptr + rate_base;

        const fp4_t piq = make_fp4(
            D.d_frequencies[0] * Qclv.x,
            D.d_frequencies[1] * Qclv.y,
            D.d_frequencies[2] * Qclv.z,
            D.d_frequencies[3] * Qclv.w);
        // MLIPPER stores PMAT factors in the standard order:
        //   P = V * diag(exp(lambda * t)) * Vinv
        // The sumtable projection must use the same basis convention:
        //   left[j]  = sum_k piq[k] * V[k,j]
        //   right[j] = sum_k Vinv[j,k] * clv[k]
        const fp4_t left_proj = matvec4_cols(D.d_V, piq);
        const fp4_t right_proj = matvec4_rows(D.d_Vinv, Pclv);

        sumtable_row[0] = left_proj.x * right_proj.x;
        sumtable_row[1] = left_proj.y * right_proj.y;
        sumtable_row[2] = left_proj.z * right_proj.z;
        sumtable_row[3] = left_proj.w * right_proj.w;

        const fp_t sum_row_max = fp_fmax(
            fp_fmax(sumtable_row[0], sumtable_row[1]),
            fp_fmax(sumtable_row[2], sumtable_row[3]));
        if (sum_row_max > fp_t(0)) {
            if (!have_signal || midpoint_shifts[r] < site_min_shift) {
                site_min_shift = midpoint_shifts[r];
            }
            active_rate_mask |= (1u << r);
            have_signal = true;
        }
    }

    if (!have_signal) return;

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        if ((active_rate_mask & (1u << r)) == 0u) continue;
        const int diff = static_cast<int>(midpoint_shifts[r]) - static_cast<int>(site_min_shift);
        if (diff <= 0) continue;
        const size_t rate_base = static_cast<size_t>(r) * 4;
        fp_t* sumtable_row = sumtable_ptr + rate_base;
        sumtable_row[0] = fp_ldexp(sumtable_row[0], -diff);
        sumtable_row[1] = fp_ldexp(sumtable_row[1], -diff);
        sumtable_row[2] = fp_ldexp(sumtable_row[2], -diff);
        sumtable_row[3] = fp_ldexp(sumtable_row[3], -diff);
    }
}

__device__ __forceinline__ void update_pendant_sumtable_site_generic(
    DeviceTree D,
    int target_id,
    const fp_t* __restrict__ left_clv_base,
    size_t site_idx,
    fp_t* sumtable,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat)
{
    if (D.rate_cats <= 0 || D.rate_cats > kMaxRateCats) return;

    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t site_span = rate_count * state_count;
    const size_t site_base = site_idx * site_span;
    const fp_t* left_clv = left_clv_base + site_base;
    fp_t* sumtable_ptr = sumtable + site_base;

    unsigned int midpoint_shifts[kMaxRateCats];
    bool rate_has_signal[kMaxRateCats];
    unsigned int site_min_shift = 0u;
    bool have_signal = false;

    for (int r = 0; r < D.rate_cats; ++r) {
        const size_t rate_offset = static_cast<size_t>(r);
        const size_t rate_base = rate_offset * 4;
        const fp4_t Qclv = reinterpret_cast<const fp4_t*>(left_clv + rate_base)[0];
        const fp4_t Pclv = build_pendant_midpoint_site_rate_generic(
            D,
            target_id,
            site_idx,
            r,
            shared_target_mat,
            shared_parent_mat,
            &midpoint_shifts[r],
            nullptr,
            nullptr,
            nullptr,
            nullptr);
        fp_t* sumtable_row = sumtable_ptr + rate_base;

        const fp4_t piq = make_fp4(
            D.d_frequencies[0] * Qclv.x,
            D.d_frequencies[1] * Qclv.y,
            D.d_frequencies[2] * Qclv.z,
            D.d_frequencies[3] * Qclv.w);
        const fp4_t left_proj = matvec4_cols(D.d_V, piq);
        const fp4_t right_proj = matvec4_rows(D.d_Vinv, Pclv);

        sumtable_row[0] = left_proj.x * right_proj.x;
        sumtable_row[1] = left_proj.y * right_proj.y;
        sumtable_row[2] = left_proj.z * right_proj.z;
        sumtable_row[3] = left_proj.w * right_proj.w;

        const fp_t sum_row_max = fp_fmax(
            fp_fmax(sumtable_row[0], sumtable_row[1]),
            fp_fmax(sumtable_row[2], sumtable_row[3]));
        rate_has_signal[r] = (sum_row_max > fp_t(0));
        if (rate_has_signal[r]) {
            if (!have_signal || midpoint_shifts[r] < site_min_shift) {
                site_min_shift = midpoint_shifts[r];
            }
            have_signal = true;
        }
    }

    if (!have_signal) return;

    for (int r = 0; r < D.rate_cats; ++r) {
        if (!rate_has_signal[r]) continue;
        const int diff = static_cast<int>(midpoint_shifts[r]) - static_cast<int>(site_min_shift);
        if (diff <= 0) continue;
        const size_t rate_base = static_cast<size_t>(r) * 4;
        fp_t* sumtable_row = sumtable_ptr + rate_base;
        sumtable_row[0] = fp_ldexp(sumtable_row[0], -diff);
        sumtable_row[1] = fp_ldexp(sumtable_row[1], -diff);
        sumtable_row[2] = fp_ldexp(sumtable_row[2], -diff);
        sumtable_row[3] = fp_ldexp(sumtable_row[3], -diff);
    }
}

// Build one site's derivative sumtable rows for the proximal-side branch update.
template<int RATE_CATS>
__device__ __forceinline__ void update_proximal_sumtable_site(
    DeviceTree D,
    int target_id,
    const fp_t* __restrict__ left_clv_base,
    const unsigned* __restrict__ left_scaler_base,
    const fp_t* __restrict__ query_clv_base,
    size_t site_idx,
    fp_t* sumtable,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat)
{
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t site_span = rate_count * state_count;
    const size_t site_base = site_idx * site_span;
    const fp_t* left_clv = left_clv_base + site_base;
    fp_t* sumtable_ptr = sumtable + site_base;
    fp4_t midpoint_rows[RATE_CATS];
    unsigned int midpoint_shifts[RATE_CATS];
    bool rate_has_signal[RATE_CATS];
    unsigned int site_min_shift = 0u;
    bool have_signal = false;

    build_proximal_midpoint_site<RATE_CATS>(
        D,
        target_id,
        site_idx,
        query_clv_base,
        shared_target_mat,
        shared_parent_mat,
        midpoint_rows,
        midpoint_shifts,
        nullptr,
        nullptr,
        nullptr,
        nullptr);

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        const size_t rate_offset = static_cast<size_t>(r);
        const size_t rate_base = rate_offset * 4;
        const fp4_t Qclv = reinterpret_cast<const fp4_t*>(left_clv + rate_base)[0];
        const fp4_t Pclv = midpoint_rows[r];
        fp_t* sumtable_row = sumtable_ptr + rate_base;

        const fp4_t piq = make_fp4(
            D.d_frequencies[0] * Qclv.x,
            D.d_frequencies[1] * Qclv.y,
            D.d_frequencies[2] * Qclv.z,
            D.d_frequencies[3] * Qclv.w);
        // Keep the derivative projection aligned with MLIPPER's PMAT basis:
        //   P = V * diag(exp(lambda * t)) * Vinv
        // so the sumtable uses V on the left and Vinv on the right.
        const fp4_t left_proj = matvec4_cols(D.d_V, piq);
        const fp4_t right_proj = matvec4_rows(D.d_Vinv, Pclv);

        sumtable_row[0] = left_proj.x * right_proj.x;
        sumtable_row[1] = left_proj.y * right_proj.y;
        sumtable_row[2] = left_proj.z * right_proj.z;
        sumtable_row[3] = left_proj.w * right_proj.w;

        const fp_t row_max = fp_fmax(
            fp_fmax(sumtable_row[0], sumtable_row[1]),
            fp_fmax(sumtable_row[2], sumtable_row[3]));
        const unsigned int total_shift =
            scaler_shift_at_site(D, left_scaler_base, site_idx, r) +
            midpoint_shifts[r];
        rate_has_signal[r] = (row_max > fp_t(0));
        if (rate_has_signal[r]) {
            if (!have_signal || total_shift < site_min_shift) {
                site_min_shift = total_shift;
            }
            midpoint_shifts[r] = total_shift;
            have_signal = true;
        }
    }

    if (!have_signal) return;

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        if (!rate_has_signal[r]) continue;
        const int diff = static_cast<int>(midpoint_shifts[r]) - static_cast<int>(site_min_shift);
        if (diff <= 0) continue;
        const size_t rate_base = static_cast<size_t>(r) * 4;
        fp_t* sumtable_row = sumtable_ptr + rate_base;
        sumtable_row[0] = fp_ldexp(sumtable_row[0], -diff);
        sumtable_row[1] = fp_ldexp(sumtable_row[1], -diff);
        sumtable_row[2] = fp_ldexp(sumtable_row[2], -diff);
        sumtable_row[3] = fp_ldexp(sumtable_row[3], -diff);
    }
}

// Reductions keep both derivatives in lockstep so the block leader observes a
// matched (df, ddf) pair before applying the shared Newton state transition.
static __device__ __forceinline__
void reduce_block_derivatives(
    double& local_df,
    double& local_ddf,
    double* warp_df_sums,
    double* warp_ddf_sums,
    double& block_df,
    double& block_ddf)
{
    const unsigned int lane = threadIdx.x & 31;
    const unsigned int warp = threadIdx.x >> 5;
    const unsigned int warp_count = (blockDim.x + 31) >> 5;
    for (int offset = 16; offset > 0; offset >>= 1) {
        local_df += __shfl_down_sync(0xffffffff, local_df, offset);
        local_ddf += __shfl_down_sync(0xffffffff, local_ddf, offset);
    }
    if (lane == 0) {
        warp_df_sums[warp] = local_df;
        warp_ddf_sums[warp] = local_ddf;
    }
    __syncthreads();
    if (warp == 0) {
        double warp_df_total = (lane < warp_count) ? warp_df_sums[lane] : 0.0;
        double warp_ddf_total = (lane < warp_count) ? warp_ddf_sums[lane] : 0.0;
        for (int offset = 16; offset > 0; offset >>= 1) {
            warp_df_total += __shfl_down_sync(0xffffffff, warp_df_total, offset);
            warp_ddf_total += __shfl_down_sync(0xffffffff, warp_ddf_total, offset);
        }
        if (lane == 0) {
            block_df = warp_df_total;
            block_ddf = warp_ddf_total;
        }
    }
    __syncthreads();
}

// Uniform safeguarded update used by the single- and multi-block tree-edge
// Newton kernels. status: 0=active, 1=converged, 2=numerical failure.
static __device__ __forceinline__ void apply_raxml_fast_newton_update(
    double df,
    double ddf,
    double tolerance,
    double maximum_step,
    double original_branch,
    double& branch,
    double& lower,
    double& upper,
    int& status)
{
    if (!isfinite(df) || !isfinite(ddf) || !isfinite(branch)) {
        branch = original_branch;
        status = 2;
        return;
    }

    double delta = 0.0;
    if (ddf > 0.0) {
        if (fabs(df) < tolerance) {
            status = 1;
            return;
        }
        if (df < 0.0) lower = branch;
        else upper = branch;
        delta = -df / ddf;
    } else if (ddf != 0.0) {
        delta = -df / fabs(ddf);
    } else if (fabs(df) < tolerance) {
        status = 1;
        return;
    } else {
        branch = original_branch;
        status = 2;
        return;
    }

    if (!isfinite(delta)) {
        branch = original_branch;
        status = 2;
        return;
    }
    delta = max(min(delta, maximum_step), -maximum_step);
    if (branch + delta < lower) delta = lower - branch;
    if (branch + delta > upper) delta = upper - branch;
    if (!isfinite(delta)) {
        branch = original_branch;
        status = 2;
    } else if (fabs(delta) < tolerance) {
        status = 1;
    } else {
        const double candidate = branch + delta;
        if (!isfinite(candidate)) {
            branch = original_branch;
            status = 2;
        } else {
            branch = clamp_scalar(candidate, kTreeMinimumLength, kMaximumLength);
        }
    }
}

static __device__ __forceinline__
void apply_newton_update(
    double& branch_value,
    int& stop_iterations,
    double& branch_lower_bound,
    double& branch_upper_bound,
    double max_step,
    double block_df,
    double block_ddf)
{
    const double tolerance = OPT_BRANCH_XTOL;

    if (!isfinite(branch_value) || !isfinite(block_df) ||
        !isfinite(block_ddf) || !isfinite(max_step)) {
        stop_iterations = 1;
        return;
    }

    double branch_delta = 0.0;
    double proposed_branch = branch_value;
    if (block_ddf > 0.0) {
        // Match EPA-ng's bracketed Newton step:
        // update the active bound based on the sign of the first derivative,
        // then use the raw Newton proposal for the step.
        if (fabs(block_df) < tolerance) {
            stop_iterations = 1;
            return;
        }
        if (block_df < 0.0) {
            branch_lower_bound = branch_value;
        } else {
            branch_upper_bound = branch_value;
        }
        const double newton_delta = -block_df / block_ddf;
        const double newton_branch = branch_value + newton_delta;
        const bool outside_bracket =
            (((branch_value - branch_upper_bound) * block_ddf - block_df) *
             ((branch_value - branch_lower_bound) * block_ddf - block_df)) >= 0.0;

        if (outside_bracket) {
            // Fall back to a conservative bracket midpoint when the Newton
            // proposal is not safely bracketed. This mirrors the "safe"
            // behavior of EPA-ng's branch-length optimizer and prevents the
            // step from repeatedly chasing the global cap.
            proposed_branch = branch_lower_bound + 0.5 * (branch_upper_bound - branch_lower_bound);
            branch_delta = proposed_branch - branch_value;
        } else {
            branch_delta = newton_delta;
            proposed_branch = newton_branch;
        }
    } else {
        // Match EPA-ng/libpll's fallback for negative curvature: preserve the
        // ascent direction from the first derivative, but damp the step using
        // the magnitude of the second derivative.
        if (fabs(block_ddf) < tolerance) {
            stop_iterations = 1;
            return;
        }
        branch_delta = -block_df / fabs(block_ddf);
        proposed_branch = branch_value + branch_delta;
    }

    double clipped_delta = max(min(branch_delta, max_step), -max_step);
    if (branch_value + clipped_delta < branch_lower_bound) {
        clipped_delta = branch_lower_bound - branch_value;
    }
    if (branch_value + clipped_delta > branch_upper_bound) {
        clipped_delta = branch_upper_bound - branch_value;
    }

    if (fabs(clipped_delta) < tolerance) {
        stop_iterations = 1;
        return;
    }

    double final_branch = branch_value + clipped_delta;
    const bool hit_lower_cap = (final_branch <= OPT_BRANCH_LEN_MIN);
    const bool hit_upper_cap = (final_branch >= OPT_BRANCH_LEN_MAX);
    if (hit_lower_cap) final_branch = OPT_BRANCH_LEN_MIN;
    if (hit_upper_cap) final_branch = OPT_BRANCH_LEN_MAX;
    branch_value = final_branch;

    // EPA-ng-style behavior: once the Newton update is forced against a hard
    // bound, stop iterating instead of repeatedly chasing the cap. This avoids
    // turning the bound itself into the apparent optimum.
    if (hit_lower_cap || hit_upper_cap) {
        stop_iterations = 1;
    }
}

static __device__ __forceinline__
unsigned int scaler_shift_at_site(
    const DeviceTree& D,
    const unsigned* __restrict__ scaler_base,
    size_t site_idx,
    size_t rate_idx)
{
    if (!scaler_base) return 0u;
    if (D.per_rate_scaling) {
        const size_t rate_count = static_cast<size_t>(D.rate_cats);
        const size_t site_base = site_idx * rate_count;
        return scaler_base[site_base + rate_idx];
    }
    return scaler_base[site_idx];
}

static __device__ __forceinline__
unsigned int scaler_shift_from_site_ptr(
    const DeviceTree& D,
    const unsigned* __restrict__ scaler_site_ptr,
    size_t rate_idx)
{
    if (!scaler_site_ptr) return 0u;
    if (D.per_rate_scaling) {
        return scaler_site_ptr[rate_idx];
    }
    return scaler_site_ptr[0];
}

static __device__ __forceinline__
fp_t load_pattern_weight_cached(
    const unsigned* __restrict__ pattern_weights,
    size_t site)
{
    return pattern_weights
        ? static_cast<fp_t>(pattern_weights[site])
        : fp_t(1);
}

// Evaluate one site's first and second derivatives of negative log likelihood
// with respect to the target edge length. Upward and outside messages are
// projected into eigen space, combined with exp(lambda*t) derivatives, aligned
// to a common scaler exponent across rate categories, and mixed with optional
// invariant-site mass. Outputs already include the compressed pattern weight.
static __device__ __forceinline__
void evaluate_tree_edge_site_derivatives_direct(
    const DeviceTree D,
    int target_id,
    const int* __restrict__ invariant_site,
    const fp_t* __restrict__ invar_proportion,
    fp_t invar_scalar,
    const fp_t* __restrict__ diag_shared,
    const unsigned* __restrict__ pattern_weights,
    size_t site,
    double& d1_out,
    double& d2_out)
{
    d1_out = 0.0;
    d2_out = 0.0;
    if (D.rate_cats <= 0 || D.rate_cats > kMaxRateCats || D.states <= 0) {
        return;
    }

    const size_t state_count = static_cast<size_t>(D.states);
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t site_span = rate_count * state_count;
    const size_t site_base = site * site_span;
    const fp_t* up_site = up_clv_ptr<const fp_t>(D, target_id, site_base);
    const fp_t* edge_outside_site = edge_outside_clv_ptr<const fp_t>(D, target_id, site_base);
    if (!up_site || !edge_outside_site || !D.d_V || !D.d_Vinv || !D.d_frequencies) {
        return;
    }

    const unsigned* up_scaler = up_scaler_ptr(D, target_id, site);
    const unsigned* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site);
    fp_t rate_lk0[kMaxRateCats];
    fp_t rate_lk1[kMaxRateCats];
    fp_t rate_lk2[kMaxRateCats];
    unsigned int rate_shifts[kMaxRateCats];
    bool rate_has_signal[kMaxRateCats];
    unsigned int site_min_shift = 0u;
    bool have_signal = false;
    const int inv = invariant_site ? invariant_site[site] : -1;

    for (int r = 0; r < D.rate_cats; ++r) {
        const size_t rate_offset = static_cast<size_t>(r) * state_count;
        const fp_t* up_row = up_site + rate_offset;
        const fp_t* edge_outside_row = edge_outside_site + rate_offset;
        const fp_t* diag_row = diag_shared + rate_offset * 4;
        fp_t cat0 = fp_t(0);
        fp_t cat1 = fp_t(0);
        fp_t cat2 = fp_t(0);

        for (int eig_idx = 0; eig_idx < D.states; ++eig_idx) {
            fp_t left_proj = fp_t(0);
            fp_t right_proj = fp_t(0);
            for (int state_idx = 0; state_idx < D.states; ++state_idx) {
                const fp_t weighted_up =
                    D.d_frequencies[state_idx] * up_row[static_cast<size_t>(state_idx)];
                left_proj = fp_fma(
                    D.d_V[static_cast<size_t>(state_idx) * state_count +
                          static_cast<size_t>(eig_idx)],
                    weighted_up,
                    left_proj);
                right_proj = fp_fma(
                    D.d_Vinv[static_cast<size_t>(eig_idx) * state_count +
                             static_cast<size_t>(state_idx)],
                    edge_outside_row[static_cast<size_t>(state_idx)],
                    right_proj);
            }

            const fp_t sum_value = left_proj * right_proj;
            const size_t diag_base = static_cast<size_t>(eig_idx) * 4;
            cat0 = fp_fma(sum_value, diag_row[diag_base + 0], cat0);
            cat1 = fp_fma(sum_value, diag_row[diag_base + 1], cat1);
            cat2 = fp_fma(
                sum_value, diag_row[diag_base + 2], cat2);
        }

        const fp_t pinv = invar_proportion ? invar_proportion[r] : invar_scalar;
        if (pinv > fp_t(0)) {
            const fp_t inv_site_lk = (inv < 0) ? fp_t(0) : (D.d_frequencies[inv] * pinv);
            const fp_t non_pinv = fp_t(1) - pinv;
            cat0 = cat0 * non_pinv + inv_site_lk;
            cat1 = cat1 * non_pinv;
            cat2 = cat2 * non_pinv;
        }

        rate_lk0[r] = cat0;
        rate_lk1[r] = cat1;
        rate_lk2[r] = cat2;
        rate_shifts[r] =
            scaler_shift_from_site_ptr(D, up_scaler, static_cast<size_t>(r)) +
            scaler_shift_from_site_ptr(D, edge_outside_scaler, static_cast<size_t>(r));
        rate_has_signal[r] = (cat0 > fp_t(0));
        if (rate_has_signal[r]) {
            if (!have_signal || rate_shifts[r] < site_min_shift) {
                site_min_shift = rate_shifts[r];
            }
            have_signal = true;
        }
    }

    if (!have_signal) {
        return;
    }

    fp_t site_lk0 = fp_t(0);
    fp_t site_lk1 = fp_t(0);
    fp_t site_lk2 = fp_t(0);
    for (int r = 0; r < D.rate_cats; ++r) {
        if (!rate_has_signal[r]) {
            continue;
        }
        const int shift_diff =
            static_cast<int>(rate_shifts[r]) - static_cast<int>(site_min_shift);
        fp_t scaled_cat0 = rate_lk0[r];
        fp_t scaled_cat1 = rate_lk1[r];
        fp_t scaled_cat2 = rate_lk2[r];
        if (shift_diff > 0) {
            scaled_cat0 = fp_ldexp(scaled_cat0, -shift_diff);
            scaled_cat1 = fp_ldexp(scaled_cat1, -shift_diff);
            scaled_cat2 = fp_ldexp(scaled_cat2, -shift_diff);
        }
        const fp_t rate_weight = D.d_rate_weights ? D.d_rate_weights[r] : fp_t(1);
        site_lk0 = fp_fma(scaled_cat0, rate_weight, site_lk0);
        site_lk1 = fp_fma(scaled_cat1, rate_weight, site_lk1);
        site_lk2 = fp_fma(
            scaled_cat2, rate_weight, site_lk2);
    }

    if (!(site_lk0 > fp_t(0))) {
        return;
    }

    const fp_t inv_lk0 = fp_t(1) / site_lk0;
    const fp_t d1 = -site_lk1 * inv_lk0;
    const fp_t d2 = d1 * d1 - (site_lk2 * inv_lk0);
    const fp_t weight = load_pattern_weight_cached(pattern_weights, site);
    d1_out = static_cast<double>(weight * d1);
    d2_out = static_cast<double>(weight * d2);
}

// The sumtable is independent of the proposed branch length. These helpers
// combine it with the current diagonal and optionally emit the weighted site
// log-likelihood used by placement scoring.
template<bool WRITE_PLACEMENT_CLV>
static __device__ __forceinline__
void evaluate_site_derivatives_general(
    const DeviceTree D,
    const int*    __restrict__ invariant_site,
    const fp_t* __restrict__ invar_proportion,
    fp_t invar_scalar,
    const fp_t* __restrict__ sumtable,
    const fp_t* __restrict__ lambdas,
    const unsigned* __restrict__ pattern_weights,
    size_t site,
    fp_t* placement_clv,
    double& d1_out,
    double& d2_out)
{
    fp_t site_lk0 = fp_t(0);
    fp_t site_lk1 = fp_t(0);
    fp_t site_lk2 = fp_t(0);

    const int inv = invariant_site ? invariant_site[site] : -1;
    const int states = D.states;
    const int rate_cats = D.rate_cats;
    const size_t state_count = static_cast<size_t>(states);
    const size_t rate_count = static_cast<size_t>(rate_cats);
    const size_t site_span = rate_count * state_count;
    const size_t site_base = site * site_span;
    const fp_t* site_sumtable = sumtable + site_base;
    const fp_t* rate_weights = D.d_rate_weights;
    const fp_t* frequencies = D.d_frequencies;

    for (int i = 0; i < rate_cats; ++i) {
        const size_t rate_base = static_cast<size_t>(i) * state_count;
        const fp_t* sum_row = site_sumtable + rate_base;
        const fp_t* diag_row = lambdas + rate_base * 4;

        fp_t cat0 = fp_t(0);
        fp_t cat1 = fp_t(0);
        fp_t cat2 = fp_t(0);

        for (int j = 0; j < states; ++j) {
            const fp_t sum_value = sum_row[j];
            cat0 = fp_fma(sum_value, diag_row[0], cat0);
            cat1 = fp_fma(sum_value, diag_row[1], cat1);
            cat2 = fp_fma(sum_value, diag_row[2], cat2);
            diag_row += 4;
        }
        const fp_t pinv = invar_proportion ? invar_proportion[i] : invar_scalar;
        if (pinv > fp_t(0)) {
            const fp_t inv_site_lk = (inv < 0) ? fp_t(0) : (frequencies[inv] * pinv);
            const fp_t non_pinv = fp_t(1) - pinv;
            cat0 = cat0 * non_pinv + inv_site_lk;
            cat1 = cat1 * non_pinv;
            cat2 = cat2 * non_pinv;
        }
        const fp_t w = rate_weights ? rate_weights[i] : fp_t(1);
        site_lk0 += cat0 * w;
        site_lk1 += cat1 * w;
        site_lk2 += cat2 * w;
    }

    const fp_t inv_lk0 = (site_lk0 != fp_t(0)) ? (fp_t(1) / site_lk0) : fp_t(0);
    const fp_t d1 = -site_lk1 * inv_lk0;
    const fp_t d2 = d1 * d1 - (site_lk2 * inv_lk0);

    const fp_t weight = load_pattern_weight_cached(pattern_weights, site);
    d1_out = static_cast<double>(weight * d1);
    d2_out = static_cast<double>(weight * d2);
    if constexpr (WRITE_PLACEMENT_CLV) {
        placement_clv[site] = fp_log(site_lk0 > FP_EPS ? site_lk0 : FP_EPS) * weight;
    }
}

template<int RATE_CATS, bool WRITE_PLACEMENT_CLV>
static __device__ __forceinline__
void evaluate_site_derivatives_states4_noinv(
    const DeviceTree D,
    const fp_t* __restrict__ sumtable,
    const fp_t* __restrict__ lambdas,
    const unsigned* __restrict__ pattern_weights,
    size_t site,
    fp_t* placement_clv,
    double& d1_out,
    double& d2_out)
{
    const size_t rate_count = static_cast<size_t>(RATE_CATS);
    const size_t site_span = rate_count * 4;
    const size_t site_base = site * site_span;
    const fp_t* site_sumtable = sumtable + site_base;
    const fp_t* rate_weights = D.d_rate_weights;
    fp_t site_lk0 = fp_t(0);
    fp_t site_lk1 = fp_t(0);
    fp_t site_lk2 = fp_t(0);
    bool eq_weights = true;
    fp_t eq_weight = fp_t(1);
    if (rate_weights && RATE_CATS > 1) {
        const fp_t w0 = rate_weights[0];
        eq_weight = w0;
#pragma unroll
        for (int r = 1; r < RATE_CATS; ++r) {
            if (rate_weights[r] != w0) {
                eq_weights = false;
                break;
            }
        }
    }

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        const size_t rate_offset = static_cast<size_t>(r);
        const size_t rate_base = rate_offset * 4;
        const fp_t* sum_row = site_sumtable + rate_base;
        const fp_t* diag_row = lambdas + rate_offset * 16;
        const fp_t rw = rate_weights ? rate_weights[r] : fp_t(1);

        const fp_t cat0 = fp_fma(
            sum_row[0], diag_row[0],
            fp_fma(sum_row[1], diag_row[4],
                   fp_fma(sum_row[2], diag_row[8],
                          sum_row[3] * diag_row[12])));
        const fp_t cat1 = fp_fma(
            sum_row[0], diag_row[1],
            fp_fma(sum_row[1], diag_row[5],
                   fp_fma(sum_row[2], diag_row[9],
                          sum_row[3] * diag_row[13])));
        const fp_t cat2 = fp_fma(
            sum_row[0], diag_row[2],
            fp_fma(sum_row[1], diag_row[6],
                   fp_fma(sum_row[2], diag_row[10],
                          sum_row[3] * diag_row[14])));

        if (eq_weights) {
            site_lk0 = fp_fma(cat0, eq_weight, site_lk0);
            site_lk1 = fp_fma(cat1, eq_weight, site_lk1);
            site_lk2 = fp_fma(cat2, eq_weight, site_lk2);
        } else {
            site_lk0 = fp_fma(cat0, rw, site_lk0);
            site_lk1 = fp_fma(cat1, rw, site_lk1);
            site_lk2 = fp_fma(cat2, rw, site_lk2);
        }
    }

    const fp_t inv_lk0 = (site_lk0 != fp_t(0)) ? (fp_t(1) / site_lk0) : fp_t(0);
    const fp_t d1 = -site_lk1 * inv_lk0;
    const fp_t d2 = d1 * d1 - (site_lk2 * inv_lk0);

    const fp_t weight = load_pattern_weight_cached(pattern_weights, site);
    d1_out = static_cast<double>(weight * d1);
    d2_out = static_cast<double>(weight * d2);
    if constexpr (WRITE_PLACEMENT_CLV) {
        placement_clv[site] = fp_log(site_lk0 > FP_EPS ? site_lk0 : FP_EPS) * weight;
    }
}

template<int RATE_CATS, bool WRITE_PLACEMENT_CLV>
static __device__ __forceinline__
void accumulate_site_derivatives_states4_noinv(
    const DeviceTree D,
    const fp_t* __restrict__ sumtable_op,
    const fp_t* __restrict__ diag_shared,
    const unsigned* __restrict__ pattern_weights,
    fp_t* placement_clv,
    unsigned int tid,
    unsigned int step,
    double& df_local,
    double& ddf_local)
{
    df_local = 0.0;
    ddf_local = 0.0;
    for (size_t site = tid; site < D.sites; site += step) {
        double d1 = 0.0;
        double d2 = 0.0;
        evaluate_site_derivatives_states4_noinv<RATE_CATS, WRITE_PLACEMENT_CLV>(
            D,
            sumtable_op,
            diag_shared,
            pattern_weights,
            site,
            placement_clv,
            d1,
            d2);
        df_local += d1;
        ddf_local += d2;
    }
}

template<bool WRITE_PLACEMENT_CLV>
static __device__ __forceinline__
void accumulate_site_derivatives(
    const DeviceTree D,
    const int* __restrict__ invariant_site,
    const fp_t* __restrict__ invar_proportion,
    fp_t invar_scalar,
    const fp_t* __restrict__ sumtable_op,
    const fp_t* __restrict__ diag_shared,
    const unsigned* __restrict__ pattern_weights,
    fp_t* placement_clv,
    unsigned int tid,
    unsigned int step,
    double& df_local,
    double& ddf_local)
{
    df_local = 0.0;
    ddf_local = 0.0;
    const bool specialized_states4_noinv =
        D.states == 4 &&
        invariant_site == nullptr &&
        invar_proportion == nullptr &&
        invar_scalar == fp_t(0);

    if (specialized_states4_noinv) {
        switch (D.rate_cats) {
            case 1:
                accumulate_site_derivatives_states4_noinv<1, WRITE_PLACEMENT_CLV>(
                    D, sumtable_op, diag_shared, pattern_weights, placement_clv, tid, step, df_local, ddf_local);
                return;
            case 4:
                accumulate_site_derivatives_states4_noinv<4, WRITE_PLACEMENT_CLV>(
                    D, sumtable_op, diag_shared, pattern_weights, placement_clv, tid, step, df_local, ddf_local);
                return;
            case 8:
                accumulate_site_derivatives_states4_noinv<8, WRITE_PLACEMENT_CLV>(
                    D, sumtable_op, diag_shared, pattern_weights, placement_clv, tid, step, df_local, ddf_local);
                return;
            default:
                break;
        }
    }

    for (size_t site = tid; site < D.sites; site += step) {
        double d1 = 0.0;
        double d2 = 0.0;
        evaluate_site_derivatives_general<WRITE_PLACEMENT_CLV>(
            D,
            invariant_site,
            invar_proportion,
            invar_scalar,
            sumtable_op,
            diag_shared,
            pattern_weights,
            site,
            placement_clv,
            d1,
            d2);
        df_local += d1;
        ddf_local += d2;
    }
}

__device__ __forceinline__ void update_proximal_sumtable_site_generic(
    DeviceTree D,
    int target_id,
    const fp_t* __restrict__ left_clv_base,
    const unsigned* __restrict__ left_scaler_base,
    const fp_t* __restrict__ query_clv_base,
    size_t site_idx,
    fp_t* sumtable,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat)
{
    if (D.rate_cats <= 0 || D.rate_cats > kMaxRateCats) return;

    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t site_span = rate_count * state_count;
    const size_t site_base = site_idx * site_span;
    const fp_t* left_clv = left_clv_base + site_base;
    fp_t* sumtable_ptr = sumtable + site_base;
    fp4_t midpoint_rows[kMaxRateCats];
    unsigned int midpoint_shifts[kMaxRateCats];
    bool rate_has_signal[kMaxRateCats];
    unsigned int site_min_shift = 0u;
    bool have_signal = false;

    build_proximal_midpoint_site_generic(
        D,
        target_id,
        site_idx,
        query_clv_base,
        shared_target_mat,
        shared_parent_mat,
        midpoint_rows,
        midpoint_shifts,
        nullptr,
        nullptr,
        nullptr,
        nullptr);

    for (int r = 0; r < D.rate_cats; ++r) {
        const size_t rate_offset = static_cast<size_t>(r);
        const size_t rate_base = rate_offset * 4;
        const fp4_t Qclv = reinterpret_cast<const fp4_t*>(left_clv + rate_base)[0];
        const fp4_t Pclv = midpoint_rows[r];
        fp_t* sumtable_row = sumtable_ptr + rate_base;

        const fp4_t piq = make_fp4(
            D.d_frequencies[0] * Qclv.x,
            D.d_frequencies[1] * Qclv.y,
            D.d_frequencies[2] * Qclv.z,
            D.d_frequencies[3] * Qclv.w);
        const fp4_t left_proj = matvec4_cols(D.d_V, piq);
        const fp4_t right_proj = matvec4_rows(D.d_Vinv, Pclv);

        sumtable_row[0] = left_proj.x * right_proj.x;
        sumtable_row[1] = left_proj.y * right_proj.y;
        sumtable_row[2] = left_proj.z * right_proj.z;
        sumtable_row[3] = left_proj.w * right_proj.w;

        const fp_t row_max = fp_fmax(
            fp_fmax(sumtable_row[0], sumtable_row[1]),
            fp_fmax(sumtable_row[2], sumtable_row[3]));
        const unsigned int total_shift =
            scaler_shift_at_site(D, left_scaler_base, site_idx, r) +
            midpoint_shifts[r];
        rate_has_signal[r] = (row_max > fp_t(0));
        if (rate_has_signal[r]) {
            if (!have_signal || total_shift < site_min_shift) {
                site_min_shift = total_shift;
            }
            midpoint_shifts[r] = total_shift;
            have_signal = true;
        }
    }

    if (!have_signal) return;

    for (int r = 0; r < D.rate_cats; ++r) {
        if (!rate_has_signal[r]) continue;
        const int diff = static_cast<int>(midpoint_shifts[r]) - static_cast<int>(site_min_shift);
        if (diff <= 0) continue;
        const size_t rate_base = static_cast<size_t>(r) * 4;
        fp_t* sumtable_row = sumtable_ptr + rate_base;
        sumtable_row[0] = fp_ldexp(sumtable_row[0], -diff);
        sumtable_row[1] = fp_ldexp(sumtable_row[1], -diff);
        sumtable_row[2] = fp_ldexp(sumtable_row[2], -diff);
        sumtable_row[3] = fp_ldexp(sumtable_row[3], -diff);
    }
}

// Each builder partitions one candidate's site range across the block. The
// resulting sumtable stays fixed while the kernel evaluates multiple Newton
// proposals for that branch coordinate.
template<int RATE_CATS>
static __device__ __forceinline__
void build_pendant_sumtable(
    DeviceTree D,
    int target_id,
    const fp_t* __restrict__ left_base,
    fp_t* sumtable_op,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat,
    unsigned int tid,
    unsigned int step)
{
    for (size_t site = tid; site < D.sites; site += step) {
        update_pendant_sumtable_site<RATE_CATS>(
            D,
            target_id,
            left_base,
            site,
            sumtable_op,
            shared_target_mat,
            shared_parent_mat);
    }
}

static __device__ __forceinline__
void build_pendant_sumtable_generic(
    DeviceTree D,
    int target_id,
    const fp_t* __restrict__ left_base,
    fp_t* sumtable_op,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat,
    unsigned int tid,
    unsigned int step)
{
    for (size_t site = tid; site < D.sites; site += step) {
        update_pendant_sumtable_site_generic(
            D,
            target_id,
            left_base,
            site,
            sumtable_op,
            shared_target_mat,
            shared_parent_mat);
    }
}

template<int RATE_CATS>
static __device__ __forceinline__
void build_proximal_sumtable(
    DeviceTree D,
    int target_id,
    const fp_t* __restrict__ left_base,
    const unsigned* __restrict__ left_scaler_base,
    const fp_t* __restrict__ query_clv_base,
    fp_t* sumtable_op,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat,
    unsigned int tid,
    unsigned int step)
{
    for (size_t site = tid; site < D.sites; site += step) {
        update_proximal_sumtable_site<RATE_CATS>(
            D,
            target_id,
            left_base,
            left_scaler_base,
            query_clv_base,
            site,
            sumtable_op,
            shared_target_mat,
            shared_parent_mat);
    }
}

static __device__ __forceinline__
void build_proximal_sumtable_generic(
    DeviceTree D,
    int target_id,
    const fp_t* __restrict__ left_base,
    const unsigned* __restrict__ left_scaler_base,
    const fp_t* __restrict__ query_clv_base,
    fp_t* sumtable_op,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat,
    unsigned int tid,
    unsigned int step)
{
    for (size_t site = tid; site < D.sites; site += step) {
        update_proximal_sumtable_site_generic(
            D,
            target_id,
            left_base,
            left_scaler_base,
            query_clv_base,
            site,
            sumtable_op,
            shared_target_mat,
            shared_parent_mat);
    }
}

// Pendant-side derivative kernel.
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
    double branch_min)
{
    if (!sumtable || !ops || op_idx < 0) {
        return;
    }

    const int op_local = op_idx + static_cast<int>(blockIdx.x);
    const int op_global = op_indices ? op_indices[op_local] : op_local;
    if (op_global < 0 || op_global >= D.N) return;
    if (active_ops && active_ops[op_local] == 0) return;
    const NodeOpInfo op = ops[op_global];
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_id = target_is_left ? op.left_id : op.right_id;
    if (target_id < 0 || target_id >= D.N) return;
    if (!D.d_query_clv || !D.d_edge_outside_clv || !D.d_clv_up || !D.d_pmat_mid_dist || !D.d_pmat_mid_prox) {
        return;
    }

    // Per-op output and branch state.
    const size_t op_base = static_cast<size_t>(op_local) * sumtable_stride;
    fp_t* sumtable_op = sumtable + op_base;
    __shared__ double branch_value_shared;
    __shared__ int stop_iterations;
    __shared__ double branch_lower_bound_shared;
    __shared__ double branch_upper_bound_shared;
    __shared__ double max_step_shared;
    double init_branch = branch_min;
    if (prev_branch_lengths) {
        init_branch = static_cast<double>(prev_branch_lengths[target_id]);
    } else {
        init_branch = DEFAULT_BRANCH_LENGTH;
    }
    if (threadIdx.x == 0) {
        if (init_branch < branch_min) init_branch = branch_min;
        if (init_branch > OPT_BRANCH_LEN_MAX) init_branch = OPT_BRANCH_LEN_MAX;
        branch_value_shared = init_branch;
        branch_lower_bound_shared = branch_min;
        branch_upper_bound_shared = OPT_BRANCH_LEN_MAX;
        max_step_shared = OPT_BRANCH_LEN_MAX / static_cast<double>(max_iter);
        stop_iterations = 0;
    }
    __syncthreads();

    // Shared-memory layout for branch diagonals and midpoint PMATs.
    extern __shared__ fp_t shmem[];
    fp_t* branch_diag_shared = shmem;
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    int query_idx = node_op_query_idx(op);
    if (query_idx < 0 || query_idx >= D.query_capacity) query_idx = 0;
    const fp_t* query_clv_base =
        D.d_query_clv +
        static_cast<size_t>(query_idx) * D.sites * rate_count * state_count;
    const size_t diag_span = rate_count * state_count * 4;
    const size_t midpoint_pmat_span = rate_count * 16;
    fp_t* target_midpoint_pmat_shared = branch_diag_shared + diag_span;
    fp_t* parent_midpoint_pmat_shared = target_midpoint_pmat_shared + midpoint_pmat_span;
    __shared__ double warp_df_sums[32];
    __shared__ double warp_ddf_sums[32];
    __shared__ double block_df_shared;
    __shared__ double block_ddf_shared;

    // Load midpoint PMATs for the current target edge.
    const size_t pmat_base = static_cast<size_t>(target_id) * midpoint_pmat_span;
    const fp_t* target_mat = D.d_pmat_mid_prox + pmat_base;
    const fp_t* parent_mat = D.d_pmat_mid_dist + pmat_base;
    if (!load_midpoint_pmat_pair_dispatch(
            D, target_midpoint_pmat_shared, parent_midpoint_pmat_shared, target_mat, parent_mat)) {
        return;
    }

    // Build the sumtable once for this operation.
    const unsigned int tid = threadIdx.x;
    const unsigned int step = blockDim.x;
    switch (D.rate_cats) {
        case 1:
            build_pendant_sumtable<1>(
                D,
                target_id,
                query_clv_base,
                sumtable_op,
                target_midpoint_pmat_shared,
                parent_midpoint_pmat_shared,
                tid,
                step);
            break;
        case 4:
            build_pendant_sumtable<4>(
                D,
                target_id,
                query_clv_base,
                sumtable_op,
                target_midpoint_pmat_shared,
                parent_midpoint_pmat_shared,
                tid,
                step);
            break;
        case 8:
            build_pendant_sumtable<8>(
                D,
                target_id,
                query_clv_base,
                sumtable_op,
                target_midpoint_pmat_shared,
                parent_midpoint_pmat_shared,
                tid,
                step);
            break;
        default:
            build_pendant_sumtable_generic(
                D,
                target_id,
                query_clv_base,
                sumtable_op,
                target_midpoint_pmat_shared,
                parent_midpoint_pmat_shared,
                tid,
                step);
            break;
    }
    __syncthreads();

    // Newton iterations over the already-built sumtable.
    double local_df = 0.0;
    double local_ddf = 0.0;
    for (int iter = 0; iter < max_iter; ++iter) {
        if (stop_iterations) break;
        const fp_t branch = static_cast<fp_t>(branch_value_shared);
        if (D.states == 4) {
            if (!build_diagtable_states4_dispatch(D, branch, branch_diag_shared)) {
                build_diagtable_for_branch(D, branch, branch_diag_shared);
            }
        } else {
            build_diagtable_for_branch(D, branch, branch_diag_shared);
        }
        __syncthreads();

        accumulate_site_derivatives<false>(
            D,
            invariant_site,
            invar_proportion,
            invar_scalar,
            sumtable_op,
            branch_diag_shared,
            pattern_weights,
            nullptr,
            tid,
            step,
            local_df,
            local_ddf);

        reduce_block_derivatives(local_df, local_ddf, warp_df_sums, warp_ddf_sums, block_df_shared, block_ddf_shared);

        if (threadIdx.x == 0) {
            apply_newton_update(
                branch_value_shared,
                stop_iterations,
                branch_lower_bound_shared,
                branch_upper_bound_shared,
                max_step_shared,
                block_df_shared,
                block_ddf_shared);
        }
        __syncthreads();
        if (stop_iterations) break;
    }

    if (threadIdx.x == 0) {
        new_branch_length[target_id] = static_cast<fp_t>(branch_value_shared);
    }
}

// Proximal-side derivative kernel.
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
    double branch_min)
{
    if (!sumtable || !ops || op_idx < 0) {
        return;
    }

    const int op_local = op_idx + static_cast<int>(blockIdx.x);
    const int op_global = op_indices ? op_indices[op_local] : op_local;
    if (op_global < 0 || op_global >= D.N) return;
    if (active_ops && active_ops[op_local] == 0) return;
    const NodeOpInfo op = ops[op_global];
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_id = target_is_left ? op.left_id : op.right_id;
    if (target_id < 0 || target_id >= D.N) return;
    if (!D.d_clv_up || !D.d_edge_outside_clv || !D.d_query_clv || !D.d_query_pmat || !D.d_pmat_mid_dist) {
        return;
    }

    // Per-op output and branch state.
    const size_t op_base = static_cast<size_t>(op_local) * sumtable_stride;
    fp_t* sumtable_op = sumtable + op_base;
    __shared__ double branch_value_shared;
    __shared__ int stop_iterations;
    __shared__ double branch_lower_bound_shared;
    __shared__ double branch_upper_bound_shared;
    __shared__ double max_step_shared;
    double init_branch = branch_min;
    const double total_branch = static_cast<double>(D.d_blen[target_id]);
    const double branch_lower_bound = effective_split_branch_min(total_branch, branch_min);
    const double branch_upper_bound = scalar_max(branch_lower_bound, total_branch - branch_lower_bound);
    if (prev_branch_lengths) {
        init_branch = static_cast<double>(prev_branch_lengths[target_id]);
    } else {
        init_branch = 0.5 * static_cast<double>(D.d_blen[target_id]);
    }
    if (threadIdx.x == 0) {
        init_branch = clamp_scalar(init_branch, branch_lower_bound, branch_upper_bound);
        branch_value_shared = init_branch;
        branch_lower_bound_shared = branch_lower_bound;
        branch_upper_bound_shared = branch_upper_bound;
        max_step_shared = scalar_max(
            OPT_BRANCH_XTOL,
            (branch_upper_bound_shared - branch_lower_bound_shared) / static_cast<double>(max_iter));
        stop_iterations = 0;
    }
    __syncthreads();

    // Shared-memory layout for branch diagonals and midpoint PMATs.
    extern __shared__ fp_t shmem[];
    fp_t* branch_diag_shared = shmem;
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t diag_span = rate_count * state_count * 4;
    const size_t midpoint_pmat_span = rate_count * 16;
    fp_t* target_midpoint_pmat_shared = branch_diag_shared + diag_span;
    fp_t* parent_midpoint_pmat_shared = target_midpoint_pmat_shared + midpoint_pmat_span;
    __shared__ double warp_df_sums[32];
    __shared__ double warp_ddf_sums[32];
    __shared__ double block_df_shared;
    __shared__ double block_ddf_shared;

    // Resolve the left-side CLV/scaler slice for this target edge.
    if (!D.d_site_scaler_up) return;
    const size_t node_site_span = D.sites * rate_count * state_count;
    const fp_t* left_base = up_clv_ptr<const fp_t>(D, target_id);
    int query_idx = node_op_query_idx(op);
    if (query_idx < 0 || query_idx >= D.query_capacity) query_idx = 0;
    const fp_t* query_clv_base =
        D.d_query_clv + static_cast<size_t>(query_idx) * node_site_span;
    const unsigned* left_scaler_base = up_scaler_ptr(D, target_id, 0);
    if (!left_base || !left_scaler_base) return;

    // Load midpoint PMATs for the current query/target pair.
    const size_t target_pmat_base = static_cast<size_t>(op_local) * midpoint_pmat_span;
    const size_t parent_pmat_base = static_cast<size_t>(target_id) * midpoint_pmat_span;
    const fp_t* target_mat = D.d_query_pmat + target_pmat_base;
    const fp_t* parent_mat = D.d_pmat_mid_dist + parent_pmat_base;
    if (!load_midpoint_pmat_pair_dispatch(
            D, target_midpoint_pmat_shared, parent_midpoint_pmat_shared, target_mat, parent_mat)) {
        return;
    }
    // Match EPA-ng's contract here: the proximal-side sumtable is prepared
    // once for the current smoothing step, then Newton iterates on that fixed
    // sumtable instead of rebuilding the distal-side PMAT inside the loop.
    const unsigned int tid = threadIdx.x;
    const unsigned int step = blockDim.x;
    switch (D.rate_cats) {
        case 1:
            build_proximal_sumtable<1>(
                D,
                target_id,
                left_base,
                left_scaler_base,
                query_clv_base,
                sumtable_op,
                target_midpoint_pmat_shared,
                parent_midpoint_pmat_shared,
                tid,
                step);
            break;
        case 4:
            build_proximal_sumtable<4>(
                D,
                target_id,
                left_base,
                left_scaler_base,
                query_clv_base,
                sumtable_op,
                target_midpoint_pmat_shared,
                parent_midpoint_pmat_shared,
                tid,
                step);
            break;
        case 8:
            build_proximal_sumtable<8>(
                D,
                target_id,
                left_base,
                left_scaler_base,
                query_clv_base,
                sumtable_op,
                target_midpoint_pmat_shared,
                parent_midpoint_pmat_shared,
                tid,
                step);
            break;
        default:
            build_proximal_sumtable_generic(
                D,
                target_id,
                left_base,
                left_scaler_base,
                query_clv_base,
                sumtable_op,
                target_midpoint_pmat_shared,
                parent_midpoint_pmat_shared,
                tid,
                step);
            break;
    }
    __syncthreads();

    // Newton iterations over the already-built sumtable.
    double local_df = 0.0;
    double local_ddf = 0.0;
    for (int iter = 0; iter < max_iter; ++iter) {
        if (stop_iterations) break;
        const fp_t branch = static_cast<fp_t>(branch_value_shared);
        if (D.states == 4) {
            if (!build_diagtable_states4_dispatch(D, branch, branch_diag_shared)) {
                build_diagtable_for_branch(D, branch, branch_diag_shared);
            }
        } else {
            build_diagtable_for_branch(D, branch, branch_diag_shared);
        }
        __syncthreads();

        accumulate_site_derivatives<false>(
            D,
            invariant_site,
            invar_proportion,
            invar_scalar,
            sumtable_op,
            branch_diag_shared,
            pattern_weights,
            nullptr,
            tid,
            step,
            local_df,
            local_ddf);

        reduce_block_derivatives(local_df, local_ddf, warp_df_sums, warp_ddf_sums, block_df_shared, block_ddf_shared);

        if (threadIdx.x == 0) {
            apply_newton_update(
                branch_value_shared,
                stop_iterations,
                branch_lower_bound_shared,
                branch_upper_bound_shared,
                max_step_shared,
                block_df_shared,
                block_ddf_shared);
        }
        __syncthreads();

        if (stop_iterations) break;
    }

    if (threadIdx.x == 0) {
        new_branch_length[target_id] = static_cast<fp_t>(branch_value_shared);
    }
}

// Propose every selected edge from the same frozen CLV/outside-message state.
// One block owns one child-endpoint edge and reduces weighted site derivatives;
// proposals are written separately so host code can accept or roll back the
// complete Jacobi group transactionally.
__global__ void TreeEdgeJacobiBranchLengthKernel(
    const DeviceTree D,
    const int* __restrict__ invariant_site,
    const fp_t* __restrict__ invar_proportion,
    fp_t invar_scalar,
    const unsigned* __restrict__ pattern_weights,
    int max_iter,
    fp_t* new_branch_length,
    const fp_t* prev_branch_lengths,
    const unsigned char* selected_edges)
{
    const int target_id = static_cast<int>(blockIdx.x);
    if (!new_branch_length || !prev_branch_lengths || target_id < 0 || target_id >= D.N) {
        return;
    }

    if ((selected_edges && !selected_edges[target_id]) ||
        target_id == D.root_id ||
        !D.d_clv_up ||
        !D.d_edge_outside_clv ||
        !D.d_site_scaler_up ||
        !D.d_edge_outside_scaler ||
        !D.d_lambdas ||
        !D.d_V ||
        !D.d_Vinv ||
        !D.d_frequencies ||
        D.rate_cats <= 0 ||
        D.rate_cats > kMaxRateCats ||
        D.states <= 0) {
        if (threadIdx.x == 0) {
            new_branch_length[target_id] =
                (target_id == D.root_id) ? fp_t(0) : prev_branch_lengths[target_id];
        }
        return;
    }

    __shared__ double branch_value_shared;
    __shared__ int stop_iterations;
    __shared__ double branch_lower_bound_shared;
    __shared__ double branch_upper_bound_shared;
    __shared__ double max_step_shared;
    double init_branch = static_cast<double>(prev_branch_lengths[target_id]);
    if (!(init_branch > 0.0)) {
        init_branch = DEFAULT_BRANCH_LENGTH;
    }
    if (threadIdx.x == 0) {
        init_branch = clamp_scalar(init_branch, OPT_BRANCH_LEN_MIN, OPT_BRANCH_LEN_MAX);
        branch_value_shared = init_branch;
        branch_lower_bound_shared = OPT_BRANCH_LEN_MIN;
        branch_upper_bound_shared = OPT_BRANCH_LEN_MAX;
        max_step_shared = OPT_BRANCH_LEN_MAX / static_cast<double>(scalar_max(max_iter, 1));
        stop_iterations = 0;
    }
    __syncthreads();

    extern __shared__ fp_t shmem[];
    fp_t* branch_diag_shared = shmem;
    __shared__ double warp_df_sums[32];
    __shared__ double warp_ddf_sums[32];
    __shared__ double block_df_shared;
    __shared__ double block_ddf_shared;
    const unsigned int tid = threadIdx.x;
    const unsigned int step = blockDim.x;

    for (int iter = 0; iter < max_iter; ++iter) {
        if (stop_iterations) {
            break;
        }
        const fp_t branch = static_cast<fp_t>(branch_value_shared);
        build_diagtable_for_branch(D, branch, branch_diag_shared);
        __syncthreads();

        double local_df = 0.0;
        double local_ddf = 0.0;
        for (size_t site = tid; site < D.sites; site += step) {
            double d1 = 0.0;
            double d2 = 0.0;
            evaluate_tree_edge_site_derivatives_direct(
                D,
                target_id,
                invariant_site,
                invar_proportion,
                invar_scalar,
                branch_diag_shared,
                pattern_weights,
                site,
                d1,
                d2);
            local_df += d1;
            local_ddf += d2;
        }

        reduce_block_derivatives(
            local_df,
            local_ddf,
            warp_df_sums,
            warp_ddf_sums,
            block_df_shared,
            block_ddf_shared);

        if (threadIdx.x == 0) {
            apply_newton_update(
                branch_value_shared,
                stop_iterations,
                branch_lower_bound_shared,
                branch_upper_bound_shared,
                max_step_shared,
                block_df_shared,
                block_ddf_shared);
        }
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        new_branch_length[target_id] = static_cast<fp_t>(branch_value_shared);
    }
}

__global__ void SingleTreeEdgeSumtableNewtonKernel(
    const DeviceTree D,
    int target_id,
    fp_t* __restrict__ sumtable,
    int* __restrict__ newton_failure,
    int max_iter)
{
    if (!sumtable || !newton_failure || target_id < 0 || target_id >= D.N ||
        target_id == D.root_id || !D.d_clv_up || !D.d_edge_outside_clv ||
        !D.d_site_scaler_up || !D.d_edge_outside_scaler ||
        !D.d_lambdas || !D.d_V || !D.d_Vinv || !D.d_frequencies ||
        !D.d_blen || D.rate_cats <= 0 || D.rate_cats > kMaxRateCats ||
        D.states <= 0) {
        return;
    }

    const unsigned int tid = threadIdx.x;
    const unsigned int step = blockDim.x;
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t site_span = state_count * rate_count;

    // The two directional CLVs exclude the optimized edge. Their eigenbasis
    // product is therefore constant for every Newton proposal on this edge.
    for (size_t site = tid; site < D.sites; site += step) {
        const size_t site_base = site * site_span;
        const fp_t* up_site =
            up_clv_ptr<const fp_t>(D, target_id, site_base);
        const fp_t* base_site =
            edge_outside_clv_ptr<const fp_t>(D, target_id, site_base);
        const unsigned* up_scaler = up_scaler_ptr(D, target_id, site);
        const unsigned* base_scaler =
            edge_outside_scaler_ptr(D, target_id, site);
        if (!up_site || !base_site) continue;

        unsigned int shifts[kMaxRateCats];
        unsigned int min_shift = 0u;
        for (int r = 0; r < D.rate_cats; ++r) {
            shifts[r] =
                scaler_shift_from_site_ptr(
                    D, up_scaler, static_cast<size_t>(r)) +
                scaler_shift_from_site_ptr(
                    D, base_scaler, static_cast<size_t>(r));
            if (r == 0 || shifts[r] < min_shift) min_shift = shifts[r];

            const size_t rate_base = static_cast<size_t>(r) * state_count;
            const fp_t* up_row = up_site + rate_base;
            const fp_t* base_row = base_site + rate_base;
            fp_t* sum_row = sumtable + site_base + rate_base;
            for (int eig = 0; eig < D.states; ++eig) {
                fp_t left = fp_t(0);
                fp_t right = fp_t(0);
                for (int state = 0; state < D.states; ++state) {
                    left = fp_fma(
                        D.d_V[static_cast<size_t>(state) * state_count + eig],
                        D.d_frequencies[state] * up_row[state], left);
                    right = fp_fma(
                        D.d_Vinv[static_cast<size_t>(eig) * state_count + state],
                        base_row[state], right);
                }
                sum_row[eig] = left * right;
            }
        }
        for (int r = 0; r < D.rate_cats; ++r) {
            const int diff = static_cast<int>(shifts[r]) -
                static_cast<int>(min_shift);
            if (diff <= 0) continue;
            fp_t* sum_row = sumtable + site_base +
                static_cast<size_t>(r) * state_count;
            for (int eig = 0; eig < D.states; ++eig) {
                sum_row[eig] = fp_ldexp(sum_row[eig], -diff);
            }
        }
    }
    __syncthreads();

    extern __shared__ fp_t branch_diag[];
    __shared__ double branch_value;
    __shared__ double lower_bound;
    __shared__ double upper_bound;
    __shared__ double max_step;
    __shared__ double original_branch;
    __shared__ int stop_iterations;
    __shared__ double warp_df[32], warp_ddf[32], block_df, block_ddf;
    if (threadIdx.x == 0) {
        const double input_branch = static_cast<double>(D.d_blen[target_id]);
        const bool valid_input_branch = isfinite(input_branch);
        original_branch = valid_input_branch
            ? clamp_scalar(input_branch, kTreeMinimumLength, kMaximumLength)
            : DEFAULT_BRANCH_LENGTH;
        branch_value = original_branch;
        lower_bound = kTreeMinimumLength;
        upper_bound = kMaximumLength;
        max_step = kMaximumLength /
            static_cast<double>(scalar_max(max_iter, 1));
        stop_iterations = valid_input_branch ? 0 : 2;
        if (!valid_input_branch) {
            atomicCAS(newton_failure, 0, target_id + 1);
        }
    }
    __syncthreads();

    for (int iteration = 0; iteration < max_iter; ++iteration) {
        if (stop_iterations) break;
        build_diagtable_for_branch(
            D, static_cast<fp_t>(branch_value), branch_diag);
        __syncthreads();
        double local_df = 0.0;
        double local_ddf = 0.0;
        accumulate_site_derivatives<false>(
            D, nullptr, nullptr, fp_t(0), sumtable, branch_diag,
            D.d_pattern_weights_u, nullptr, tid, step,
            local_df, local_ddf);
        reduce_block_derivatives(
            local_df, local_ddf, warp_df, warp_ddf, block_df, block_ddf);
        if (threadIdx.x == 0) {
            // Match corax_opt_minimize_newton_multi() FAST exactly. In
            // particular, do not use the placement SAFE midpoint fallback:
            // RAxML-NG 0.9.0 uses the raw clipped Newton displacement.
            apply_raxml_fast_newton_update(
                block_df, block_ddf, kNewtonTolerance, max_step,
                original_branch, branch_value, lower_bound, upper_bound,
                stop_iterations);
            if (stop_iterations == 2) {
                atomicCAS(newton_failure, 0, target_id + 1);
            }
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        D.d_blen[target_id] = static_cast<fp_t>(branch_value);
    }
}

// One branch remains one optimization coordinate, but its site work is spread
// across several resident blocks. Cooperative grid barriers keep every block
// on the same Newton proposal without returning partial derivatives to the
// host or launching a separate kernel for each Newton iteration.
__global__ void SingleTreeEdgeSumtableNewtonMultiBlockKernel(
    const DeviceTree D,
    int target_id,
    fp_t* __restrict__ sumtable,
    double* __restrict__ partial_gradient,
    double* __restrict__ partial_hessian,
    double* __restrict__ newton_state,
    int* __restrict__ newton_failure,
    int max_iter)
{
    namespace cg = cooperative_groups;
    cg::grid_group grid = cg::this_grid();
    if (!sumtable || !partial_gradient || !partial_hessian ||
        !newton_state || !newton_failure ||
        target_id < 0 || target_id >= D.N ||
        target_id == D.root_id || !D.d_clv_up || !D.d_edge_outside_clv ||
        !D.d_site_scaler_up || !D.d_edge_outside_scaler ||
        !D.d_lambdas || !D.d_V || !D.d_Vinv || !D.d_frequencies ||
        !D.d_blen || D.rate_cats <= 0 || D.rate_cats > kMaxRateCats ||
        D.states <= 0) {
        return;
    }

    const size_t global_tid =
        static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const size_t global_step =
        static_cast<size_t>(gridDim.x) * blockDim.x;
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t site_span = state_count * rate_count;

    for (size_t site = global_tid; site < D.sites; site += global_step) {
        const size_t site_base = site * site_span;
        const fp_t* up_site =
            up_clv_ptr<const fp_t>(D, target_id, site_base);
        const fp_t* base_site =
            edge_outside_clv_ptr<const fp_t>(D, target_id, site_base);
        const unsigned* up_scaler = up_scaler_ptr(D, target_id, site);
        const unsigned* base_scaler =
            edge_outside_scaler_ptr(D, target_id, site);
        if (!up_site || !base_site) continue;

        unsigned int shifts[kMaxRateCats];
        unsigned int min_shift = 0u;
        for (int r = 0; r < D.rate_cats; ++r) {
            shifts[r] =
                scaler_shift_from_site_ptr(
                    D, up_scaler, static_cast<size_t>(r)) +
                scaler_shift_from_site_ptr(
                    D, base_scaler, static_cast<size_t>(r));
            if (r == 0 || shifts[r] < min_shift) min_shift = shifts[r];

            const size_t rate_base = static_cast<size_t>(r) * state_count;
            const fp_t* up_row = up_site + rate_base;
            const fp_t* base_row = base_site + rate_base;
            fp_t* sum_row = sumtable + site_base + rate_base;
            for (int eig = 0; eig < D.states; ++eig) {
                fp_t left = fp_t(0);
                fp_t right = fp_t(0);
                for (int state = 0; state < D.states; ++state) {
                    left = fp_fma(
                        D.d_V[static_cast<size_t>(state) * state_count + eig],
                        D.d_frequencies[state] * up_row[state], left);
                    right = fp_fma(
                        D.d_Vinv[static_cast<size_t>(eig) * state_count + state],
                        base_row[state], right);
                }
                sum_row[eig] = left * right;
            }
        }
        for (int r = 0; r < D.rate_cats; ++r) {
            const int diff = static_cast<int>(shifts[r]) -
                static_cast<int>(min_shift);
            if (diff <= 0) continue;
            fp_t* sum_row = sumtable + site_base +
                static_cast<size_t>(r) * state_count;
            for (int eig = 0; eig < D.states; ++eig) {
                sum_row[eig] = fp_ldexp(sum_row[eig], -diff);
            }
        }
    }
    grid.sync();

    // state: branch value, lower bound, upper bound, maximum step, status,
    // and the original value used for numerical-failure rollback.
    if (blockIdx.x == 0 && threadIdx.x == 0) {
        const double input_branch = static_cast<double>(D.d_blen[target_id]);
        const bool valid_input_branch = isfinite(input_branch);
        newton_state[5] = valid_input_branch
            ? clamp_scalar(input_branch, kTreeMinimumLength, kMaximumLength)
            : DEFAULT_BRANCH_LENGTH;
        newton_state[0] = newton_state[5];
        newton_state[1] = kTreeMinimumLength;
        newton_state[2] = kMaximumLength;
        newton_state[3] = kMaximumLength /
            static_cast<double>(scalar_max(max_iter, 1));
        newton_state[4] = valid_input_branch ? 0.0 : 2.0;
        if (!valid_input_branch) {
            atomicCAS(newton_failure, 0, target_id + 1);
        }
    }
    grid.sync();

    extern __shared__ fp_t branch_diag[];
    __shared__ double warp_df[32], warp_ddf[32], block_df, block_ddf;
    for (int iteration = 0; iteration < max_iter; ++iteration) {
        // status is written by one grid leader before the preceding grid
        // barrier. Every thread therefore takes this exit together; no block
        // can skip a later grid.sync() while another block enters it.
        if (newton_state[4] != 0.0) break;
        build_diagtable_for_branch(
            D, static_cast<fp_t>(newton_state[0]), branch_diag);
        __syncthreads();
        double local_df = 0.0;
        double local_ddf = 0.0;
        accumulate_site_derivatives<false>(
            D, nullptr, nullptr, fp_t(0), sumtable, branch_diag,
            D.d_pattern_weights_u, nullptr, global_tid, global_step,
            local_df, local_ddf);
        reduce_block_derivatives(
            local_df, local_ddf, warp_df, warp_ddf, block_df, block_ddf);
        if (threadIdx.x == 0) {
            partial_gradient[blockIdx.x] = block_df;
            partial_hessian[blockIdx.x] = block_ddf;
        }
        grid.sync();

        if (blockIdx.x == 0 && threadIdx.x == 0) {
            double df = 0.0;
            double ddf = 0.0;
            for (unsigned int block = 0; block < gridDim.x; ++block) {
                df += partial_gradient[block];
                ddf += partial_hessian[block];
            }
            double& branch = newton_state[0];
            double& lower = newton_state[1];
            double& upper = newton_state[2];
            int status = static_cast<int>(newton_state[4]);
            apply_raxml_fast_newton_update(
                df, ddf, kNewtonTolerance, newton_state[3],
                newton_state[5], branch, lower, upper, status);
            newton_state[4] = static_cast<double>(status);
            if (status == 2) {
                atomicCAS(newton_failure, 0, target_id + 1);
            }
        }
        grid.sync();
    }
    if (blockIdx.x == 0 && threadIdx.x == 0) {
        D.d_blen[target_id] = static_cast<fp_t>(newton_state[0]);
    }
}

__global__ void SingleTreeEdgeSumtableNewtonWarpSiteKernel(
    const DeviceTree D,
    int target_id,
    fp_t* __restrict__ sumtable,
    double* __restrict__ partial_gradient,
    double* __restrict__ partial_hessian,
    double* __restrict__ newton_state,
    int* __restrict__ newton_failure,
    int max_iter)
{
    namespace cg = cooperative_groups;
    constexpr int kLanesPerSite = 16;
    constexpr int kStates = 4;
    constexpr int kRates = 4;
    cg::grid_group grid = cg::this_grid();
    if (!sumtable || !partial_gradient || !partial_hessian ||
        !newton_state || !newton_failure ||
        target_id < 0 || target_id >= D.N ||
        target_id == D.root_id || D.states != kStates ||
        D.rate_cats != kRates || !D.d_clv_up || !D.d_edge_outside_clv ||
        !D.d_site_scaler_up || !D.d_edge_outside_scaler ||
        !D.d_lambdas || !D.d_V || !D.d_Vinv || !D.d_frequencies ||
        !D.d_blen) {
        return;
    }

    const unsigned int lane = threadIdx.x & 31u;
    const unsigned int site_lane = lane & (kLanesPerSite - 1);
    const unsigned int halfwarp_base = lane & ~(kLanesPerSite - 1);
    const unsigned int halfwarp_mask =
        (halfwarp_base == 0) ? 0x0000ffffu : 0xffff0000u;
    const size_t global_thread =
        static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const size_t global_group = global_thread / kLanesPerSite;
    const size_t group_step =
        static_cast<size_t>(gridDim.x) * blockDim.x / kLanesPerSite;
    const int rate = static_cast<int>(site_lane / kStates);
    const int eigen = static_cast<int>(site_lane % kStates);

    for (size_t site = global_group; site < D.sites; site += group_step) {
        const size_t site_base = site * kRates * kStates;
        const size_t rate_base = static_cast<size_t>(rate) * kStates;
        const fp_t* up_row =
            up_clv_ptr<const fp_t>(D, target_id, site_base + rate_base);
        const fp_t* base_row =
            edge_outside_clv_ptr<const fp_t>(D, target_id, site_base + rate_base);
        const unsigned* up_scaler = up_scaler_ptr(D, target_id, site);
        const unsigned* base_scaler =
            edge_outside_scaler_ptr(D, target_id, site);
        if (!up_row || !base_row) continue;

        fp_t left = fp_t(0);
        fp_t right = fp_t(0);
#pragma unroll
        for (int state = 0; state < kStates; ++state) {
            left = fp_fma(
                D.d_V[static_cast<size_t>(state) * kStates + eigen],
                D.d_frequencies[state] * up_row[state], left);
            right = fp_fma(
                D.d_Vinv[static_cast<size_t>(eigen) * kStates + state],
                base_row[state], right);
        }

        const unsigned int shift =
            scaler_shift_from_site_ptr(D, up_scaler, rate) +
            scaler_shift_from_site_ptr(D, base_scaler, rate);
        unsigned int min_shift = shift;
#pragma unroll
        for (int offset = 8; offset > 0; offset >>= 1) {
            min_shift = min(
                min_shift,
                __shfl_down_sync(
                    halfwarp_mask, min_shift, offset, kLanesPerSite));
        }
        min_shift = __shfl_sync(
            halfwarp_mask, min_shift, halfwarp_base);
        const int shift_difference =
            static_cast<int>(shift) - static_cast<int>(min_shift);
        fp_t value = left * right;
        if (shift_difference > 0) {
            value = fp_ldexp(value, -shift_difference);
        }
        sumtable[site_base + rate_base + eigen] = value;
    }
    grid.sync();

    if (blockIdx.x == 0 && threadIdx.x == 0) {
        const double input_branch = static_cast<double>(D.d_blen[target_id]);
        const bool valid_input_branch = isfinite(input_branch);
        newton_state[5] = valid_input_branch
            ? clamp_scalar(input_branch, kTreeMinimumLength, kMaximumLength)
            : DEFAULT_BRANCH_LENGTH;
        newton_state[0] = newton_state[5];
        newton_state[1] = kTreeMinimumLength;
        newton_state[2] = kMaximumLength;
        newton_state[3] = kMaximumLength /
            static_cast<double>(scalar_max(max_iter, 1));
        newton_state[4] = valid_input_branch ? 0.0 : 2.0;
        if (!valid_input_branch) {
            atomicCAS(newton_failure, 0, target_id + 1);
        }
    }
    grid.sync();

    extern __shared__ fp_t branch_diag[];
    __shared__ double warp_df[32], warp_ddf[32], block_df, block_ddf;
    for (int iteration = 0; iteration < max_iter; ++iteration) {
        // Keep all cooperative blocks on the same control-flow path across
        // every grid barrier.
        if (newton_state[4] != 0.0) break;
        build_diagtable_for_branch(
            D, static_cast<fp_t>(newton_state[0]), branch_diag);
        __syncthreads();

        double local_df = 0.0;
        double local_ddf = 0.0;
        for (size_t site = global_group; site < D.sites; site += group_step) {
            const size_t component =
                (site * kRates + static_cast<size_t>(rate)) * kStates + eigen;
            const size_t diag_component =
                (static_cast<size_t>(rate) * kStates + eigen) * 4;
            const fp_t sum_value = sumtable[component];
            const fp_t rate_weight = D.d_rate_weights
                ? D.d_rate_weights[rate]
                : fp_t(1);
            fp_t likelihood =
                sum_value * branch_diag[diag_component + 0] * rate_weight;
            fp_t first =
                sum_value * branch_diag[diag_component + 1] * rate_weight;
            fp_t second =
                sum_value * branch_diag[diag_component + 2] * rate_weight;
#pragma unroll
            for (int offset = 8; offset > 0; offset >>= 1) {
                likelihood += __shfl_down_sync(
                    halfwarp_mask, likelihood, offset, kLanesPerSite);
                first += __shfl_down_sync(
                    halfwarp_mask, first, offset, kLanesPerSite);
                second += __shfl_down_sync(
                    halfwarp_mask, second, offset, kLanesPerSite);
            }
            if (site_lane == 0) {
                const fp_t inverse = likelihood != fp_t(0)
                    ? fp_t(1) / likelihood
                    : fp_t(0);
                const fp_t d1 = -first * inverse;
                const fp_t d2 = d1 * d1 - second * inverse;
                const fp_t weight = load_pattern_weight_cached(
                    D.d_pattern_weights_u, site);
                local_df += static_cast<double>(weight * d1);
                local_ddf += static_cast<double>(weight * d2);
            }
        }
        reduce_block_derivatives(
            local_df, local_ddf, warp_df, warp_ddf, block_df, block_ddf);
        if (threadIdx.x == 0) {
            partial_gradient[blockIdx.x] = block_df;
            partial_hessian[blockIdx.x] = block_ddf;
        }
        grid.sync();

        if (blockIdx.x == 0 && threadIdx.x == 0) {
            double df = 0.0;
            double ddf = 0.0;
            for (unsigned int block = 0; block < gridDim.x; ++block) {
                df += partial_gradient[block];
                ddf += partial_hessian[block];
            }
            double& branch = newton_state[0];
            double& lower = newton_state[1];
            double& upper = newton_state[2];
            int status = static_cast<int>(newton_state[4]);
            apply_raxml_fast_newton_update(
                df, ddf, kNewtonTolerance, newton_state[3],
                newton_state[5], branch, lower, upper, status);
            newton_state[4] = static_cast<double>(status);
            if (status == 2) {
                atomicCAS(newton_failure, 0, target_id + 1);
            }
        }
        grid.sync();
    }
    if (blockIdx.x == 0 && threadIdx.x == 0) {
        D.d_blen[target_id] = static_cast<fp_t>(newton_state[0]);
    }
}

void OptimizeSingleTreeEdgeFromCurrentClvs(
    const DeviceTree& D,
    int target_id,
    fp_t* d_sumtable,
    double* d_partial_gradient,
    double* d_partial_hessian,
    double* d_newton_state,
    int* d_newton_failure,
    int max_iter,
    cudaStream_t stream)
{
    if (max_iter <= 0) return;
    ensure_device_tree_current_device(
        D, "OptimizeSingleTreeEdgeFromCurrentClvs");
    if (!d_sumtable || !d_newton_failure ||
        target_id < 0 || target_id >= D.N || target_id == D.root_id ||
        !D.d_clv_up || !D.d_edge_outside_clv ||
        !D.d_site_scaler_up || !D.d_edge_outside_scaler ||
        !D.d_lambdas || !D.d_V || !D.d_Vinv || !D.d_frequencies ||
        !D.d_blen || D.rate_cats <= 0 || D.rate_cats > kMaxRateCats ||
        D.states <= 0) {
        throw std::invalid_argument(
            "OptimizeSingleTreeEdgeFromCurrentClvs: invalid tree edge or buffers");
    }
    constexpr int kBlockSize = 256;
    const size_t shared_bytes = sizeof(fp_t) *
        static_cast<size_t>(D.rate_cats) * D.states * 4;

    struct GenericLaunchCache {
        int device = -1;
        size_t shared_bytes = 0;
        int cooperative_launch = 0;
        int maximum_blocks = 1;
    };
    static thread_local GenericLaunchCache launch_cache;
    int current_device = -1;
    CUDA_CHECK(cudaGetDevice(&current_device));
    if (launch_cache.device != current_device ||
        launch_cache.shared_bytes != shared_bytes) {
        int multiprocessors = 0;
        int active_blocks_per_sm = 0;
        CUDA_CHECK(cudaDeviceGetAttribute(
            &launch_cache.cooperative_launch,
            cudaDevAttrCooperativeLaunch, current_device));
        CUDA_CHECK(cudaDeviceGetAttribute(
            &multiprocessors, cudaDevAttrMultiProcessorCount,
            current_device));
        CUDA_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &active_blocks_per_sm,
            SingleTreeEdgeSumtableNewtonMultiBlockKernel,
            kBlockSize,
            shared_bytes));
        launch_cache.device = current_device;
        launch_cache.shared_bytes = shared_bytes;
        launch_cache.maximum_blocks =
            active_blocks_per_sm * multiprocessors;
    }
    const int automatic_blocks = static_cast<int>(scalar_min<size_t>(
        32, (D.sites + 127) / 128));
    const int site_blocks = scalar_max(
        1, scalar_min(
            scalar_min(automatic_blocks, D.N),
            launch_cache.maximum_blocks));
    if (launch_cache.cooperative_launch && site_blocks > 1 &&
        d_partial_gradient &&
        d_partial_hessian && d_newton_state) {
        void* args[] = {
            const_cast<DeviceTree*>(&D), &target_id, &d_sumtable,
            &d_partial_gradient, &d_partial_hessian, &d_newton_state,
            &d_newton_failure, &max_iter};
        CUDA_CHECK(cudaLaunchCooperativeKernel(
            reinterpret_cast<void*>(
                SingleTreeEdgeSumtableNewtonMultiBlockKernel),
            dim3(site_blocks), dim3(kBlockSize), args,
            shared_bytes, stream));
        CUDA_CHECK(cudaGetLastError());
    } else {
        SingleTreeEdgeSumtableNewtonKernel
            <<<1, kBlockSize, shared_bytes, stream>>>(
                D, target_id, d_sumtable, d_newton_failure, max_iter);
        CUDA_CHECK(cudaGetLastError());
    }
}

void OptimizeSingleTreeEdgeFromCurrentClvsWarpSite(
    const DeviceTree& D,
    int target_id,
    fp_t* d_sumtable,
    double* d_partial_gradient,
    double* d_partial_hessian,
    double* d_newton_state,
    int* d_newton_failure,
    int max_iter,
    cudaStream_t stream)
{
    if (max_iter <= 0) return;
    ensure_device_tree_current_device(
        D, "OptimizeSingleTreeEdgeFromCurrentClvsWarpSite");
    if (D.states != 4 || D.rate_cats != 4) {
        OptimizeSingleTreeEdgeFromCurrentClvs(
            D, target_id, d_sumtable, d_partial_gradient,
            d_partial_hessian, d_newton_state, d_newton_failure,
            max_iter, stream);
        return;
    }
    if (!d_sumtable || !d_partial_gradient || !d_partial_hessian ||
        !d_newton_state || !d_newton_failure ||
        target_id < 0 || target_id >= D.N || target_id == D.root_id ||
        !D.d_clv_up || !D.d_edge_outside_clv ||
        !D.d_site_scaler_up || !D.d_edge_outside_scaler ||
        !D.d_lambdas || !D.d_V || !D.d_Vinv || !D.d_frequencies ||
        !D.d_blen) {
        throw std::invalid_argument(
            "OptimizeSingleTreeEdgeFromCurrentClvsWarpSite: invalid tree edge or buffers");
    }

    constexpr int kBlockSize = 256;
    constexpr int kLanesPerSite = 16;
    const size_t shared_bytes = sizeof(fp_t) *
        static_cast<size_t>(D.rate_cats) * D.states * 4;
    struct WarpSiteLaunchCache {
        int device = -1;
        size_t sites = 0;
        int cooperative_launch = 0;
        int maximum_blocks = 1;
    };
    static thread_local WarpSiteLaunchCache launch_cache;
    int current_device = -1;
    CUDA_CHECK(cudaGetDevice(&current_device));
    if (launch_cache.device != current_device ||
        launch_cache.sites != D.sites) {
        int multiprocessors = 0;
        int active_blocks_per_sm = 0;
        CUDA_CHECK(cudaDeviceGetAttribute(
            &launch_cache.cooperative_launch,
            cudaDevAttrCooperativeLaunch, current_device));
        CUDA_CHECK(cudaDeviceGetAttribute(
            &multiprocessors, cudaDevAttrMultiProcessorCount,
            current_device));
        CUDA_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &active_blocks_per_sm,
            SingleTreeEdgeSumtableNewtonWarpSiteKernel,
            kBlockSize,
            shared_bytes));
        launch_cache.device = current_device;
        launch_cache.sites = D.sites;
        launch_cache.maximum_blocks =
            active_blocks_per_sm * multiprocessors;
    }
    const size_t required_blocks =
        (D.sites * kLanesPerSite + kBlockSize - 1) / kBlockSize;
    const int requested_blocks = static_cast<int>(scalar_min<size_t>(
        required_blocks,
        static_cast<size_t>(std::numeric_limits<int>::max())));
    int site_blocks = scalar_max(
        1, scalar_min(
            scalar_min(requested_blocks, launch_cache.maximum_blocks), D.N));
    if (!launch_cache.cooperative_launch || site_blocks <= 1) {
        OptimizeSingleTreeEdgeFromCurrentClvs(
            D, target_id, d_sumtable, d_partial_gradient,
            d_partial_hessian, d_newton_state, d_newton_failure,
            max_iter, stream);
        return;
    }

    void* args[] = {
        const_cast<DeviceTree*>(&D), &target_id, &d_sumtable,
        &d_partial_gradient, &d_partial_hessian, &d_newton_state,
        &d_newton_failure, &max_iter};
    CUDA_CHECK(cudaLaunchCooperativeKernel(
        reinterpret_cast<void*>(SingleTreeEdgeSumtableNewtonWarpSiteKernel),
        dim3(site_blocks), dim3(kBlockSize), args,
        shared_bytes, stream));
    CUDA_CHECK(cudaGetLastError());
}

mlipper::BranchLengthOptimizationResult
RunSelectedTreeEdgeJacobiBranchLengthOptimization(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    const std::vector<int>& edge_child_ids,
    cudaStream_t stream,
    const mlipper::BranchOptimizationOptions& options)
{
    if (options.update_scheme != mlipper::EdgeUpdateScheme::Jacobi ||
        options.acceptance_scope != mlipper::AcceptanceScope::LocalSubtree) {
        throw std::invalid_argument(
            "selected-edge optimizer requires Jacobi/local-subtree options");
    }
    mlipper::BranchLengthOptimizationResult result;
    result.sweeps = options.sweeps;
    result.newton_iterations = scalar_max(options.newton_iterations, 1);
    if (options.sweeps <= 0 || edge_child_ids.empty()) return result;

    ensure_device_tree_current_device(
        D, "RunSelectedTreeEdgeJacobiBranchLengthOptimization");
    const TreeBuildResult original_tree = T;
    const HostPacking original_host = H;
    result.attempted = true;
    result.log_likelihood_before = mlipper::likelihood::root::compute_root_loglikelihood(
        D, T.root_id, nullptr, 0.0, stream);

    // Jacobi acceptance is transactional across the CPU tree, its packed host
    // representation, and all derived device PMAT/CLV state.
    auto restore_original_state = [&]() {
        T = original_tree;
        H = original_host;
        if (options.clv_retention == mlipper::ClvRetention::PreserveTips) {
            reload_device_tree_live_data_preserving_clvs(
                D, T, H, nullptr, stream);
            UpdateTreeClvsPreservingTipClvs(
                D, T, H, placement_ops, stream);
        } else {
            reload_device_tree_live_data(D, T, H, nullptr, stream);
            UpdateTreeClvs(D, T, H, placement_ops, stream);
        }
    };

    try {
        std::vector<unsigned char> selected(static_cast<size_t>(D.N), 0);
        for (int edge_child_id : edge_child_ids) {
            if (edge_child_id >= 0 && edge_child_id < D.N &&
                edge_child_id != T.root_id &&
                T.nodes[static_cast<size_t>(edge_child_id)].parent >= 0) {
                selected[static_cast<size_t>(edge_child_id)] = 1;
            }
        }
        mlipper::gpu::DeviceBuffer<unsigned char> selected_device;
        selected_device.ensureCapacity(selected.size());
        CUDA_CHECK(cudaMemcpyAsync(
            selected_device.get(), selected.data(), selected.size(),
            cudaMemcpyHostToDevice, stream));

        dim3 block(256);
        dim3 grid(static_cast<unsigned int>(D.N));
        const size_t shmem_bytes = sizeof(fp_t) *
            static_cast<size_t>(D.rate_cats) *
            static_cast<size_t>(D.states) * 4;
        std::vector<fp_t> host_branch_lengths(static_cast<size_t>(D.N));
        double previous = result.log_likelihood_before;
        for (int sweep = 0; sweep < options.sweeps; ++sweep) {
            TreeEdgeJacobiBranchLengthKernel
                <<<grid, block, shmem_bytes, stream>>>(
                    D, nullptr, nullptr, fp_t(0), D.d_pattern_weights_u,
                    result.newton_iterations, D.d_new_proximal_length,
                    D.d_blen, selected_device.get());
            CUDA_CHECK(cudaGetLastError());
            CUDA_CHECK(cudaMemcpyAsync(
                host_branch_lengths.data(), D.d_new_proximal_length,
                sizeof(fp_t) * static_cast<size_t>(D.N),
                cudaMemcpyDeviceToHost, stream));
            CUDA_CHECK(cudaStreamSynchronize(stream));

            bool any_change = false;
            for (int node_id : edge_child_ids) {
                if (node_id < 0 || node_id >= D.N || !selected[node_id]) {
                    continue;
                }
                const fp_t updated = static_cast<fp_t>(sanitize_branch_length(
                    static_cast<double>(host_branch_lengths[node_id]),
                    OPT_BRANCH_LEN_MIN, OPT_BRANCH_LEN_MAX,
                    DEFAULT_BRANCH_LENGTH));
                any_change = any_change ||
                    std::abs(static_cast<double>(
                        T.nodes[static_cast<size_t>(node_id)]
                            .branch_length_to_parent) - updated) >
                        mlipper::optimization::branch_lengths::
                            kLengthChangeTolerance;
                T.nodes[static_cast<size_t>(node_id)]
                    .branch_length_to_parent = updated;
                H.blen[static_cast<size_t>(node_id)] = updated;
            }
            if (!any_change) break;
            fill_pmats_in_host_packing(
                T, H, er, rate_multipliers, D.states, D.rate_cats);
            if (options.clv_retention == mlipper::ClvRetention::PreserveTips) {
                reload_device_tree_live_data_preserving_clvs(
                    D, T, H, nullptr, stream);
                UpdateTreeClvsPreservingTipClvs(
                    D, T, H, placement_ops, stream);
            } else {
                reload_device_tree_live_data(D, T, H, nullptr, stream);
                UpdateTreeClvs(D, T, H, placement_ops, stream);
            }
            const double current =
                mlipper::likelihood::root::compute_root_loglikelihood(
                    D, T.root_id, nullptr, 0.0, stream);
            previous = current;
            if (current >= result.log_likelihood_before &&
                current - result.log_likelihood_before <
                    options.likelihood_tolerance) {
                break;
            }
        }
        result.log_likelihood_after = previous;
        result.accepted = std::isfinite(result.log_likelihood_after) &&
            result.log_likelihood_after >= result.log_likelihood_before -
                mlipper::optimization::branch_lengths::
                    kCandidateAcceptanceTolerance;
        if (!result.accepted) {
            restore_original_state();
            result.log_likelihood_after =
                mlipper::likelihood::root::compute_root_loglikelihood(
                    D, T.root_id, nullptr, 0.0, stream);
        }
    } catch (...) {
        restore_original_state();
        throw;
    }
    return result;
}

static void release_sequential_branch_optimization_scratch(
    PlacementOpBuffer& prepared_upward,
    PlacementOpBuffer& prepared_downward,
    cudaStream_t stream)
{
    free_placement_op_buffer(prepared_upward, stream);
    free_placement_op_buffer(prepared_downward, stream);
}

mlipper::BranchLengthOptimizationResult
RunAcceptedFullTreeSequentialBranchLengthOptimization(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    const PlacementQueryBatch* queries,
    cudaStream_t stream,
    const mlipper::BranchOptimizationOptions& options)
{
    if (options.update_scheme != mlipper::EdgeUpdateScheme::Sequential ||
        options.clv_retention != mlipper::ClvRetention::RebuildAll ||
        options.acceptance_scope != mlipper::AcceptanceScope::FullTree) {
        throw std::invalid_argument(
            "full-tree optimizer requires sequential/rebuild/full-tree options");
    }
    mlipper::BranchLengthOptimizationResult result;
    result.sweeps = options.sweeps;
    result.newton_iterations = scalar_max(options.newton_iterations, 1);
    if (options.sweeps <= 0 || D.N <= 1) return result;

    ensure_device_tree_current_device(
        D, "RunAcceptedFullTreeSequentialBranchLengthOptimization");
    const int artificial_root_child =
        T.nodes[static_cast<size_t>(T.root_id)].left;
    const int optimized_root_child =
        T.nodes[static_cast<size_t>(T.root_id)].right;
    if (artificial_root_child < 0 || optimized_root_child < 0) {
        throw std::runtime_error(
            "resident sequential BLO requires a bifurcating root");
    }

    // Canonicalize the two rooted halves of the same unrooted edge. Only the
    // second half is optimized; the artificial half remains exactly zero.
    const fp_t combined_root_length = static_cast<fp_t>(
        static_cast<double>(T.nodes[static_cast<size_t>(artificial_root_child)]
                                .branch_length_to_parent) +
        static_cast<double>(T.nodes[static_cast<size_t>(optimized_root_child)]
                                .branch_length_to_parent));
    T.nodes[static_cast<size_t>(artificial_root_child)]
        .branch_length_to_parent = fp_t(0);
    T.nodes[static_cast<size_t>(optimized_root_child)]
        .branch_length_to_parent = combined_root_length;
    H.blen[static_cast<size_t>(artificial_root_child)] = fp_t(0);
    H.blen[static_cast<size_t>(optimized_root_child)] = combined_root_length;
    fill_pmats_in_host_packing(
        T, H, er, rate_multipliers, D.states, D.rate_cats);
    reload_device_tree_live_data(D, T, H, queries, stream);
    UpdateTreeClvs(D, T, H, placement_ops, stream);

    const TreeBuildResult original_tree = T;
    const HostPacking original_host = H;
    result.attempted = true;
    result.log_likelihood_before = mlipper::likelihood::root::compute_root_loglikelihood(
        D, T.root_id, nullptr, 0.0, stream);

    const auto start = std::chrono::steady_clock::now();
    DeviceTree optimization_device = D;
    optimization_device.downward_pmat_indexing =
        DownwardPmatIndexing::Rows;
    PlacementOpBuffer prepared_upward;
    PlacementOpBuffer prepared_downward;
    mlipper::gpu::DeviceBuffer<fp_t> sumtable;
    mlipper::gpu::DeviceBuffer<double> gradient;
    mlipper::gpu::DeviceBuffer<double> hessian;
    mlipper::gpu::DeviceBuffer<double> newton_state;
    mlipper::gpu::DeviceBuffer<int> newton_failure;
    try {
        PrepareTreeClvOperations(
            T, H, prepared_upward, prepared_downward, stream);
        UpdateTreeClvsPrepared(
            optimization_device,
            prepared_upward,
            prepared_downward,
            stream);

        std::vector<int> upward_index(
            static_cast<size_t>(optimization_device.N), -1);
        for (size_t index = 0;
             index < prepared_upward.upward_ops_host.size(); ++index) {
            const int node_id =
                prepared_upward.upward_ops_host[index].parent_id;
            if (node_id >= 0 && node_id < optimization_device.N) {
                upward_index[static_cast<size_t>(node_id)] =
                    static_cast<int>(index);
            }
        }
        std::vector<int> downward_index(
            static_cast<size_t>(optimization_device.N), -1);
        for (size_t index = 0;
             index < prepared_downward.downward_ops_host.size(); ++index) {
            const NodeOpInfo& op = prepared_downward.downward_ops_host[index];
            const bool target_is_left =
                op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT);
            const int target_id = target_is_left ? op.left_id : op.right_id;
            if (target_id >= 0 && target_id < optimization_device.N) {
                downward_index[static_cast<size_t>(target_id)] =
                    static_cast<int>(index);
            }
        }

        const size_t sumtable_elems =
            optimization_device.sites *
            static_cast<size_t>(optimization_device.rate_cats) *
            optimization_device.states;
        sumtable.ensureCapacity(sumtable_elems);
        gradient.ensureCapacity(optimization_device.N);
        hessian.ensureCapacity(optimization_device.N);
        newton_state.ensureCapacity(6);
        newton_failure.ensureCapacity(1);
        CUDA_CHECK(cudaMemsetAsync(
            newton_failure.get(), 0, sizeof(int), stream));

        const std::vector<SequentialBranchTraversalStep> traversal =
            build_sequential_branch_traversal_steps(T);
        for (int sweep = 0; sweep < options.sweeps; ++sweep) {
            if (sweep > 0) {
                UpdateTreeClvsDownwardOnlyPrepared(
                    optimization_device, prepared_downward, stream);
            }
            for (const SequentialBranchTraversalStep& step : traversal) {
                if (step.kind ==
                    SequentialBranchTraversalStepKind::RefreshParentUpward) {
                    const int up_op =
                        upward_index[static_cast<size_t>(step.node_id)];
                    if (up_op >= 0) {
                        UpdateSingleTreeClvUpwardWarpSitePrepared(
                            optimization_device,
                            prepared_upward,
                            up_op,
                            stream);
                    }
                    continue;
                }
                const int child_id = step.node_id;
                const int down_op =
                    downward_index[static_cast<size_t>(child_id)];
                if (down_op < 0) {
                    throw std::runtime_error(
                        "resident sequential BLO is missing a downward operation");
                }
                BuildSingleTreeEdgeOutsideWarpSitePrepared(
                    optimization_device, prepared_downward, down_op, stream);
                if (child_id != artificial_root_child) {
                    OptimizeSingleTreeEdgeFromCurrentClvsWarpSite(
                        optimization_device, child_id,
                        sumtable.get(), gradient.get(), hessian.get(),
                        newton_state.get(), newton_failure.get(),
                        result.newton_iterations, stream);
                    build_single_branch_pmat_gpu(
                        child_id,
                        optimization_device.states,
                        optimization_device.rate_cats,
                        optimization_device.d_blen,
                        optimization_device.d_V,
                        optimization_device.d_Vinv,
                        optimization_device.d_lambdas,
                        optimization_device.d_pmat,
                        stream);
                }
                if (!T.nodes[static_cast<size_t>(child_id)].is_tip) {
                    RefreshSingleTreeChildDownPrepared(
                        optimization_device,
                        prepared_downward,
                        down_op,
                        child_id,
                        stream);
                }
            }
        }

        // One independent rebuild audits all incremental directional updates.
        UpdateTreeClvsPrepared(
            optimization_device,
            prepared_upward,
            prepared_downward,
            stream);
        result.log_likelihood_after = mlipper::likelihood::root::compute_root_loglikelihood(
            optimization_device,
            T.root_id,
            nullptr, 0.0, stream);
        std::vector<fp_t> branch_lengths(
            static_cast<size_t>(optimization_device.N));
        int first_newton_failure = 0;
        CUDA_CHECK(cudaMemcpyAsync(
            branch_lengths.data(), optimization_device.d_blen,
            static_cast<size_t>(optimization_device.N) * sizeof(fp_t),
            cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(
            &first_newton_failure, newton_failure.get(), sizeof(int),
            cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));
        if (first_newton_failure != 0) {
            throw std::runtime_error(
                "resident sequential BLO encountered non-finite Newton state at node " +
                std::to_string(first_newton_failure - 1));
        }
        for (int node_id = 0; node_id < optimization_device.N; ++node_id) {
            T.nodes[static_cast<size_t>(node_id)].branch_length_to_parent =
                branch_lengths[static_cast<size_t>(node_id)];
            H.blen[static_cast<size_t>(node_id)] =
                branch_lengths[static_cast<size_t>(node_id)];
        }
        release_sequential_branch_optimization_scratch(
            prepared_upward,
            prepared_downward,
            stream);
    } catch (...) {
        release_sequential_branch_optimization_scratch(
            prepared_upward,
            prepared_downward,
            stream);
        T = original_tree;
        H = original_host;
        reload_device_tree_live_data(D, T, H, queries, stream);
        UpdateTreeClvs(D, T, H, placement_ops, stream);
        throw;
    }
    result.derivative_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    result.accepted = std::isfinite(result.log_likelihood_after) &&
        result.log_likelihood_after >= result.log_likelihood_before -
            mlipper::optimization::branch_lengths::
                kCandidateAcceptanceTolerance;
    if (result.accepted) {
        UpdateTreeClvs(D, T, H, placement_ops, stream);
        return result;
    }

    T = original_tree;
    H = original_host;
    reload_device_tree_live_data(D, T, H, queries, stream);
    UpdateTreeClvs(D, T, H, placement_ops, stream);
    result.log_likelihood_after = mlipper::likelihood::root::compute_root_loglikelihood(
        D, T.root_id, nullptr, 0.0, stream);
    return result;
}
