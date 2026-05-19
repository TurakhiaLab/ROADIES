#include "root_likelihood.cuh"
#include <cuda_runtime.h>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "tree/tree.hpp"
#include "partial_likelihood.cuh"

namespace root_likelihood {

constexpr double kLn2 = 0.69314718055994530942;
constexpr int kMaxRateCats = 8;

__device__ __forceinline__ double warp_reduce_sum(double val)
{
    for (int offset = 16; offset > 0; offset >>= 1) {
        val += __shfl_down_sync(0xffffffff, val, offset);
    }
    return val;
}

__device__ __forceinline__ double block_reduce_sum(double val)
{
    __shared__ double shared[32]; // One slot per warp for up to 1024 threads.
    const int lane = threadIdx.x & 31;
    const int wid = threadIdx.x >> 5;
    const int active_warps = (blockDim.x + 31) >> 5;

    val = warp_reduce_sum(val);
    if (lane == 0) shared[wid] = val;
    __syncthreads();

    if (wid == 0) {
        val = (lane < active_warps) ? shared[lane] : 0.0;
        val = warp_reduce_sum(val);
        return val;
    }
    return 0.0;
}

__device__ __forceinline__ int target_id_from_op(const NodeOpInfo& op)
{
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const bool target_is_right = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_RIGHT));
    return target_is_left ? op.left_id : (target_is_right ? op.right_id : op.parent_id);
}

__device__ __forceinline__ unsigned int root_scaler_shift_at(
    const DeviceTree& D,
    size_t site_idx,
    size_t rate_idx,
    size_t rate_cats)
{
    if (!D.d_site_scaler) return 0u;
    if (D.per_rate_scaling) {
        return D.d_site_scaler[site_idx * rate_cats + rate_idx];
    }
    return D.d_site_scaler[site_idx];
}

__device__ __forceinline__ double finalize_root_site_loglik(
    double sum_rate,
    const fp_t* freqs,
    size_t site_idx,
    const unsigned* pattern_w,
    const int* invar_indices,
    double invar_proportion)
{
    double site_sum = (1.0 - invar_proportion) * sum_rate;
    if (invar_indices) {
        const int inv_idx = invar_indices[site_idx];
        if (inv_idx >= 0) site_sum += invar_proportion * freqs[inv_idx];
    }
    double loglk = log(site_sum > 1e-300 ? site_sum : 1e-300);
    loglk *= static_cast<double>(pattern_w ? pattern_w[site_idx] : 1u);
    return loglk;
}

__device__ __forceinline__ unsigned int placement_scaler_shift_at(
    const unsigned* scaler_pool,
    const DeviceTree& D,
    int node_id,
    size_t site_idx,
    size_t rate_idx)
{
    if (!scaler_pool || node_id < 0) return 0u;
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t per_node = D.per_rate_scaling
        ? (D.sites * rate_count)
        : D.sites;
    const size_t base = static_cast<size_t>(node_id) * per_node;
    if (D.per_rate_scaling) {
        return scaler_pool[base + site_idx * rate_count + rate_idx];
    }
    return scaler_pool[base + site_idx];
}

// Root likelihood helpers.
// Specialized device helpers for common state/rate counts.
template<int RC>
__device__ __forceinline__ double compute_root_loglikelihood_states4(
    const DeviceTree& D,
    const fp_t* clv_site,
    const fp_t* freqs, const fp_t* rate_weights,
    size_t site_idx,
    const unsigned* pattern_w,
    const int* invar_indices,
    double invar_proportion)
{
    const double pi0 = static_cast<double>(freqs[0]);
    const double pi1 = static_cast<double>(freqs[1]);
    const double pi2 = static_cast<double>(freqs[2]);
    const double pi3 = static_cast<double>(freqs[3]);

    double sum_rate = 0.0;
    #pragma unroll
    for (int r = 0; r < RC; ++r) {
        const fp4_t a = reinterpret_cast<const fp4_t*>(clv_site)[r];
        double val = fma(static_cast<double>(a.x), pi0,
                     fma(static_cast<double>(a.y), pi1,
                     fma(static_cast<double>(a.z), pi2, static_cast<double>(a.w) * pi3)));
        const unsigned int shift = root_scaler_shift_at(
            D,
            site_idx,
            static_cast<size_t>(r),
            static_cast<size_t>(RC));
        if (shift) val = ldexp(val, -static_cast<int>(shift));
        sum_rate = fma(static_cast<double>(rate_weights[r]), val, sum_rate);
    }
    return finalize_root_site_loglik(
        sum_rate,
        freqs,
        site_idx,
        pattern_w,
        invar_indices,
        invar_proportion);
}

template<int RC>
__device__ __forceinline__ double compute_root_loglikelihood_states5(
    const DeviceTree& D,
    const fp_t* clv_site,
    const fp_t* freqs, const fp_t* rate_weights,
    size_t site_idx,
    const unsigned* pattern_w,
    const int* invar_indices,
    double invar_proportion)
{
    const double pi0 = static_cast<double>(freqs[0]);
    const double pi1 = static_cast<double>(freqs[1]);
    const double pi2 = static_cast<double>(freqs[2]);
    const double pi3 = static_cast<double>(freqs[3]);
    const double pi4 = static_cast<double>(freqs[4]);

    double sum_rate = 0.0;
    #pragma unroll
    for (int r = 0; r < RC; ++r) {
        const fp_t* cr = clv_site + static_cast<size_t>(r) * 5;
        double val = static_cast<double>(cr[0])*pi0 + static_cast<double>(cr[1])*pi1
                   + static_cast<double>(cr[2])*pi2 + static_cast<double>(cr[3])*pi3
                   + static_cast<double>(cr[4])*pi4;
        const unsigned int shift = root_scaler_shift_at(
            D,
            site_idx,
            static_cast<size_t>(r),
            static_cast<size_t>(RC));
        if (shift) val = ldexp(val, -static_cast<int>(shift));
        sum_rate = fma(static_cast<double>(rate_weights[r]), val, sum_rate);
    }
    return finalize_root_site_loglik(
        sum_rate,
        freqs,
        site_idx,
        pattern_w,
        invar_indices,
        invar_proportion);
}

// Generic device root log-likelihood for any state/rate counts (fallback).
__device__ __forceinline__ double compute_root_loglikelihood_generic(
    const DeviceTree& D,
    const fp_t* clv_site,
    const fp_t* freqs, const fp_t* rate_weights,
    size_t site_idx, const unsigned* pattern_w,
    const int* invar_indices, double invar_proportion)
{
    double sum_rate = 0.0;
    const int rate_cats = D.rate_cats;
    const int states = D.states;
    const size_t state_count = static_cast<size_t>(states);
    const size_t rate_count = static_cast<size_t>(rate_cats);
    for (int r = 0; r < rate_cats; ++r) {
        const fp_t* cr = clv_site + static_cast<size_t>(r) * state_count;
        double val = 0.0;
        for (int s = 0; s < states; ++s) {
            val = fma(static_cast<double>(cr[s]), static_cast<double>(freqs[s]), val);
        }
        const unsigned int shift = root_scaler_shift_at(D, site_idx, static_cast<size_t>(r), rate_count);
        if (shift) val = ldexp(val, -static_cast<int>(shift));
        sum_rate = fma(static_cast<double>(rate_weights[r]), val, sum_rate);
    }
    return finalize_root_site_loglik(
        sum_rate,
        freqs,
        site_idx,
        pattern_w,
        invar_indices,
        invar_proportion);
}

// Device helper that allows explicit site index (usable from arbitrary kernels).
__device__ double compute_root_loglikelihood_at_site(
    const DeviceTree& D,
    const NodeOpInfo& op,
    const fp_t* freqs,
    const fp_t* rate_weights,
    const unsigned* pattern_w,
    const int* invar_indices,
    double invar_proportion,
    size_t site_idx)
{
    if (site_idx >= D.sites) return 0.0;
    const int target_id = target_id_from_op(op);

    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t per_site = rate_count * state_count;
    const size_t per_node = D.sites * per_site;
    const fp_t* clv_pool = D.d_clv_mid;
    if (!clv_pool || target_id < 0 || target_id >= D.N) return 0.0;
    const fp_t* clv_site = clv_pool + static_cast<size_t>(target_id) * per_node + site_idx * per_site;
    if (D.states == 4) {
        switch (D.rate_cats) {
            case 1:
                return compute_root_loglikelihood_states4<1>(
                    D, clv_site, freqs, rate_weights, site_idx, pattern_w, invar_indices, invar_proportion);
            case 4:
                return compute_root_loglikelihood_states4<4>(
                    D, clv_site, freqs, rate_weights, site_idx, pattern_w, invar_indices, invar_proportion);
            case 8:
                return compute_root_loglikelihood_states4<8>(
                    D, clv_site, freqs, rate_weights, site_idx, pattern_w, invar_indices, invar_proportion);
            default:
                return compute_root_loglikelihood_generic(
                    D, clv_site, freqs, rate_weights, site_idx, pattern_w, invar_indices, invar_proportion);
        }
    }
    if (D.states == 5) {
        switch (D.rate_cats) {
            case 1:
                return compute_root_loglikelihood_states5<1>(
                    D, clv_site, freqs, rate_weights, site_idx, pattern_w, invar_indices, invar_proportion);
            case 4:
                return compute_root_loglikelihood_states5<4>(
                    D, clv_site, freqs, rate_weights, site_idx, pattern_w, invar_indices, invar_proportion);
            case 8:
                return compute_root_loglikelihood_states5<8>(
                    D, clv_site, freqs, rate_weights, site_idx, pattern_w, invar_indices, invar_proportion);
            default:
                return compute_root_loglikelihood_generic(
                    D, clv_site, freqs, rate_weights, site_idx, pattern_w, invar_indices, invar_proportion);
        }
    }
    return compute_root_loglikelihood_generic(D, clv_site, freqs, rate_weights, site_idx, pattern_w, invar_indices, invar_proportion);
}

__global__ void RootLikelihoodTotalKernel(
    DeviceTree D,
    NodeOpInfo op,
    const unsigned* __restrict__ d_pattern_w,
    const int* __restrict__ d_invar_indices,
    double invar_proportion)
{
    double local_sum = 0.0;
    const size_t block_start =
        static_cast<size_t>(blockIdx.x) * static_cast<size_t>(blockDim.x) +
        static_cast<size_t>(threadIdx.x);
    const size_t stride = static_cast<size_t>(gridDim.x) * static_cast<size_t>(blockDim.x);
    for (size_t site = block_start; site < D.sites; site += stride) {
        local_sum += compute_root_loglikelihood_at_site(
            D,
            op,
            D.d_frequencies,
            D.d_rate_weights,
            d_pattern_w,
            d_invar_indices,
            invar_proportion,
            site);
    }

    const double block_sum = block_reduce_sum(local_sum);
    if (threadIdx.x == 0) {
        atomicAdd(D.d_root_loglik_total, block_sum);
    }
}

double compute_root_loglikelihood_total(
    const DeviceTree& D,
    int root_id,
    const unsigned* d_pattern_w,
    const int* d_invar_indices,
    double invar_proportion,
    cudaStream_t stream)
{
    if (root_id < 0 || root_id >= D.N) {
        throw std::runtime_error("Invalid root id.");
    }
    if (!D.d_frequencies || !D.d_rate_weights) {
        throw std::runtime_error("Device frequencies or rate weights are not initialized.");
    }
    if (!D.d_clv_mid || !D.d_clv_up) {
        throw std::runtime_error("Device CLV buffers are not initialized.");
    }
    if (!D.d_root_loglik_total) {
        throw std::runtime_error("Device root log-likelihood accumulator is not initialized.");
    }

    DeviceTree D_root = D;
    const size_t root_node = static_cast<size_t>(root_id);
    if (D.d_site_scaler_up) {
        const size_t scaler_span = D.scaler_elems();
        // Root likelihood still uses the legacy flat site-scaler view, but it
        // should be a root-local slice, not a global alias to the UP pool.
        D_root.d_site_scaler = D.d_site_scaler_up + root_node * scaler_span;
    }

    const size_t per_node = D.per_node_elems();
    if (per_node > 0) {
        D_root.d_clv_mid = D.d_clv_up + root_node * per_node;
    }

    NodeOpInfo root_op{};
    root_op.parent_id = 0;
    root_op.dir_tag = static_cast<uint8_t>(CLV_DIR_UP);
    root_op.clv_pool = static_cast<uint8_t>(CLV_POOL_UP);

    if (D_root.sites == 0) {
        return 0.0;
    }

    dim3 block(256);
    const unsigned int grid_x = static_cast<unsigned int>((D_root.sites + block.x - 1) / block.x);
    dim3 grid(grid_x);
    double total = 0.0;
    CUDA_CHECK(cudaMemsetAsync(D_root.d_root_loglik_total, 0, sizeof(double), stream));
    RootLikelihoodTotalKernel<<<grid, block, 0, stream>>>(
        D_root,
        root_op,
        d_pattern_w,
        d_invar_indices,
        invar_proportion);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(
        &total,
        D_root.d_root_loglik_total,
        sizeof(double),
        cudaMemcpyDeviceToHost,
        stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    return total;
}

// Placement log-likelihood helpers.
template<int RATE_CATS>
__device__ __forceinline__ double placement_site_loglk_sum_ratecat(
    const fp_t* rate_vals,
    const unsigned int* rate_shifts,
    unsigned int site_min_shift,
    const fp_t* rate_weights,
    fp_t site_weight)
{
    fp_t site_lk = fp_t(0);
    #pragma unroll
    for (int rc = 0; rc < RATE_CATS; ++rc) {
        fp_t val = rate_vals[rc];
        if (val > fp_t(0)) {
            const unsigned int diff = rate_shifts[rc] - site_min_shift;
            if (diff) val = fp_ldexp(val, -static_cast<int>(diff));
            site_lk += rate_weights[rc] * val;
        }
    }
    return static_cast<double>(site_weight) *
        (static_cast<double>(fp_log(site_lk > FP_EPS ? site_lk : FP_EPS))
         - static_cast<double>(site_min_shift) * kLn2);
}

__device__ __forceinline__ double placement_site_loglk_sum_generic(
    const fp_t* rate_vals,
    const unsigned int* rate_shifts,
    int rate_cats,
    unsigned int site_min_shift,
    const fp_t* rate_weights,
    fp_t site_weight)
{
    fp_t site_lk = fp_t(0);
    for (int rc = 0; rc < rate_cats; ++rc) {
        fp_t val = rate_vals[rc];
        if (val > fp_t(0)) {
            const unsigned int diff = rate_shifts[rc] - site_min_shift;
            if (diff) val = fp_ldexp(val, -static_cast<int>(diff));
            site_lk += rate_weights[rc] * val;
        }
    }
    return static_cast<double>(site_weight) *
        (static_cast<double>(fp_log(site_lk > FP_EPS ? site_lk : FP_EPS))
         - static_cast<double>(site_min_shift) * kLn2);
}

// RateCat contract: describe CLV/PMAT layout and compute one unweighted,
// unscaled per-rate likelihood value for the common kernel.
template<int RATE_CATS>
struct PlacementRateCatStates4 {
    static constexpr bool kRuntimeRateCats = false;
    static constexpr int kRateCatsStorage = RATE_CATS;

    __device__ __forceinline__ static fp_t compute_rate_value(
        const DeviceTree& D,
        int rc,
        const fp_t* query_clv,
        const fp_t* distal_clv,
        const fp_t* prox_clv,
        const fp_t* pendant_pmat,
        const fp_t* distal_pmat,
        const fp_t* prox_pmat)
    {
        const size_t rate_offset = static_cast<size_t>(rc) * 4;
        const size_t matrix_offset = static_cast<size_t>(rc) * 16;
        const fp4_t q = reinterpret_cast<const fp4_t*>(query_clv  + rate_offset)[0];
        const fp4_t d = reinterpret_cast<const fp4_t*>(distal_clv + rate_offset)[0];
        const fp4_t p = reinterpret_cast<const fp4_t*>(prox_clv   + rate_offset)[0];
        const fp4_t* p_pendant = reinterpret_cast<const fp4_t*>(pendant_pmat + matrix_offset);
        const fp4_t* p_distal  = reinterpret_cast<const fp4_t*>(distal_pmat  + matrix_offset);
        const fp4_t* p_prox    = reinterpret_cast<const fp4_t*>(prox_pmat    + matrix_offset);
        const fp_t* freqs = D.d_frequencies;

        const fp_t acc_pend0 = fp_dot4(p_pendant[0], q);
        const fp_t acc_pend1 = fp_dot4(p_pendant[1], q);
        const fp_t acc_pend2 = fp_dot4(p_pendant[2], q);
        const fp_t acc_pend3 = fp_dot4(p_pendant[3], q);

        const fp_t acc_dist0 = fp_dot4(p_distal[0], d);
        const fp_t acc_dist1 = fp_dot4(p_distal[1], d);
        const fp_t acc_dist2 = fp_dot4(p_distal[2], d);
        const fp_t acc_dist3 = fp_dot4(p_distal[3], d);

        const fp_t acc_prox0 = fp_dot4(p_prox[0], p);
        const fp_t acc_prox1 = fp_dot4(p_prox[1], p);
        const fp_t acc_prox2 = fp_dot4(p_prox[2], p);
        const fp_t acc_prox3 = fp_dot4(p_prox[3], p);

        // Match libpll's edge-likelihood accumulation: sum all state
        // contributions directly, then reconcile per-rate scalers at the
        // site level. These terms are expected to be non-negative because
        // they are built from PMAT rows, CLVs, and stationary freqs.
        const fp_t v0 = acc_pend0 * acc_dist0 * acc_prox0 * freqs[0];
        const fp_t v1 = acc_pend1 * acc_dist1 * acc_prox1 * freqs[1];
        const fp_t v2 = acc_pend2 * acc_dist2 * acc_prox2 * freqs[2];
        const fp_t v3 = acc_pend3 * acc_dist3 * acc_prox3 * freqs[3];
        return ((v0 + v1) + (v2 + v3));
    }
};

struct PlacementRateCatGeneric {
    static constexpr bool kRuntimeRateCats = true;
    static constexpr int kRateCatsStorage = kMaxRateCats;

    __device__ __forceinline__ static fp_t compute_rate_value(
        const DeviceTree& D,
        int rc,
        const fp_t* query_clv,
        const fp_t* distal_clv,
        const fp_t* prox_clv,
        const fp_t* pendant_pmat,
        const fp_t* distal_pmat,
        const fp_t* prox_pmat)
    {
        const int states = D.states;
        const size_t state_count = static_cast<size_t>(states);
        const size_t matrix_elems = state_count * state_count;
        const size_t rate = static_cast<size_t>(rc);
        const fp_t* p_pendant = pendant_pmat + rate * matrix_elems;
        const fp_t* p_distal  = distal_pmat  + rate * matrix_elems;
        const fp_t* p_prox    = prox_pmat    + rate * matrix_elems;
        const fp_t* qrow = query_clv + rate * state_count;
        const fp_t* drow = distal_clv + rate * state_count;
        const fp_t* prow = prox_clv   + rate * state_count;

        fp_t rate_sum = fp_t(0);
        for (int s = 0; s < states; ++s) {
            const size_t row = static_cast<size_t>(s) * state_count;
            fp_t acc_pend = fp_t(0);
            fp_t acc_dist = fp_t(0);
            fp_t acc_prox = fp_t(0);
            for (int k = 0; k < states; ++k) {
                const size_t idx = row + static_cast<size_t>(k);
                acc_pend = fp_fma(p_pendant[idx], qrow[k], acc_pend);
                acc_dist = fp_fma(p_distal[idx], drow[k], acc_dist);
                acc_prox = fp_fma(p_prox[idx], prow[k], acc_prox);
            }
            rate_sum = fp_fma(acc_pend * acc_dist * acc_prox, D.d_frequencies[s], rate_sum);
        }
        return rate_sum;
    }
};

template<typename RateCat>
__global__ void PlacementLoglkKernel(
    DeviceTree D,
    const NodeOpInfo* __restrict__ d_ops,
    const int* __restrict__ d_op_indices,
    const fp_t* __restrict__ d_pendant_pmats,
    const fp_t* __restrict__ d_distal_pmats,
    const fp_t* __restrict__ d_proximal_pmats,
    size_t per_query,
    size_t per_node_pmat,
    fp_t* __restrict__ d_out)
{
    const int op_local = static_cast<int>(blockIdx.y);
    const int op_idx = d_op_indices ? d_op_indices[op_local] : op_local;
    if (!d_ops || op_idx < 0 || op_idx >= D.N) return;
    const NodeOpInfo op = d_ops[op_idx];
    const int target_id = target_id_from_op(op);
    if (target_id < 0 || target_id >= D.N) return;

    const size_t op_offset = static_cast<size_t>(op_idx);
    const size_t target_offset = static_cast<size_t>(target_id);
    const fp_t* pendant_pmat = d_pendant_pmats ? d_pendant_pmats + op_offset * per_query : nullptr;
    const fp_t* distal_pmat  = d_distal_pmats ? d_distal_pmats + target_offset * per_node_pmat : nullptr;
    const fp_t* prox_pmat    = d_proximal_pmats ? d_proximal_pmats + target_offset * per_node_pmat : nullptr;
    if (!pendant_pmat || !distal_pmat || !prox_pmat) return;

    if (!D.d_query_clv || !D.d_clv_mid_base || !D.d_clv_up || !D.d_rate_weights || !D.d_frequencies) return;
    constexpr int kRateCatsStorage = RateCat::kRateCatsStorage;
    const int rate_cat_count = RateCat::kRuntimeRateCats ? D.rate_cats : kRateCatsStorage;
    const size_t rate_count = static_cast<size_t>(rate_cat_count);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t per_site = RateCat::kRuntimeRateCats
        ? (rate_count * state_count)
        : (static_cast<size_t>(kRateCatsStorage) * 4);
    const size_t per_node = D.sites * per_site;

    double local_sum = 0.0;
    const size_t site_step = static_cast<size_t>(blockDim.x);
    for (size_t site = static_cast<size_t>(threadIdx.x); site < D.sites; site += site_step) {
        const fp_t site_weight = static_cast<fp_t>(
            D.d_pattern_weights_u ? D.d_pattern_weights_u[site] : 1u);
        const size_t site_offset = site * per_site;
        const fp_t* query_clv = D.d_query_clv + site_offset;
        const fp_t* distal_clv = D.d_clv_mid_base + target_offset * per_node + site_offset;
        const fp_t* prox_clv = D.d_clv_up + target_offset * per_node + site_offset;

        fp_t rate_vals[kRateCatsStorage];
        unsigned int rate_shifts[kRateCatsStorage];
        unsigned int site_min_shift = 0u;
        bool have_positive = false;
        for (int rc = 0; rc < rate_cat_count; ++rc) {
            const size_t rate_idx = static_cast<size_t>(rc);
            const unsigned int distal_shift =
                placement_scaler_shift_at(D.d_site_scaler_mid_base, D, target_id, site, rate_idx);
            const unsigned int prox_shift =
                placement_scaler_shift_at(D.d_site_scaler_up, D, target_id, site, rate_idx);
            const fp_t rate_sum = RateCat::compute_rate_value(
                D,
                rc,
                query_clv,
                distal_clv,
                prox_clv,
                pendant_pmat,
                distal_pmat,
                prox_pmat);
            rate_vals[rc] = rate_sum;
            rate_shifts[rc] = distal_shift + prox_shift;
            if (rate_sum > fp_t(0)) {
                if (!have_positive || rate_shifts[rc] < site_min_shift) {
                    site_min_shift = rate_shifts[rc];
                }
                have_positive = true;
            }
        }
        if constexpr (RateCat::kRuntimeRateCats) {
            local_sum += placement_site_loglk_sum_generic(
                rate_vals,
                rate_shifts,
                rate_cat_count,
                site_min_shift,
                D.d_rate_weights,
                site_weight);
        } else {
            local_sum += placement_site_loglk_sum_ratecat<kRateCatsStorage>(
                rate_vals,
                rate_shifts,
                site_min_shift,
                D.d_rate_weights,
                site_weight);
        }
    }

    const double block_sum = block_reduce_sum(local_sum);
    if (threadIdx.x == 0) d_out[op_idx] = static_cast<fp_t>(block_sum);
}

void Placement_Root_Loglk(
    const DeviceTree& D,
    const NodeOpInfo* d_ops,
    const int* d_op_indices,
    int num_ops,
    const fp_t* d_pendant_pmats,
    const fp_t* d_distal_pmats,
    const fp_t* d_proximal_pmats,
    fp_t* d_out,
    cudaStream_t stream)
{
    if (num_ops <= 0) return;
    if (!d_ops || !d_pendant_pmats || !d_distal_pmats || !d_proximal_pmats) {
        throw std::runtime_error("Missing PMAT or ops pointers for placement loglk.");
    }
    if (!d_out) {
        throw std::runtime_error("Missing output buffer for placement loglk.");
    }
    if (!D.d_query_clv || !D.d_clv_mid_base || !D.d_clv_up || !D.d_rate_weights || !D.d_frequencies) {
        throw std::runtime_error("Placement buffers not initialized.");
    }

    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t per_query = rate_count * state_count * state_count;
    const size_t per_node_pmat = per_query;
    auto check_launch = [&](const char* stage) {
        const cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string(stage) + ": " + cudaGetErrorString(err));
        }
    };

    dim3 block(256);
    dim3 grid(1, static_cast<unsigned int>(num_ops));
    if (D.states == 4) {
        switch (D.rate_cats) {
            case 1:
                PlacementLoglkKernel<PlacementRateCatStates4<1>><<<grid, block, 0, stream>>>(
                    D, d_ops, d_op_indices, d_pendant_pmats, d_distal_pmats, d_proximal_pmats,
                    per_query, per_node_pmat, d_out);
                check_launch("PlacementLoglkKernel<PlacementRateCatStates4<1>>");
                return;
            case 4:
                PlacementLoglkKernel<PlacementRateCatStates4<4>><<<grid, block, 0, stream>>>(
                    D, d_ops, d_op_indices, d_pendant_pmats, d_distal_pmats, d_proximal_pmats,
                    per_query, per_node_pmat, d_out);
                check_launch("PlacementLoglkKernel<PlacementRateCatStates4<4>>");
                return;
            case 8:
                PlacementLoglkKernel<PlacementRateCatStates4<8>><<<grid, block, 0, stream>>>(
                    D, d_ops, d_op_indices, d_pendant_pmats, d_distal_pmats, d_proximal_pmats,
                    per_query, per_node_pmat, d_out);
                check_launch("PlacementLoglkKernel<PlacementRateCatStates4<8>>");
                return;
            default:
                break;
        }
    }

    PlacementLoglkKernel<PlacementRateCatGeneric><<<grid, block, 0, stream>>>(
        D,
        d_ops,
        d_op_indices,
        d_pendant_pmats,
        d_distal_pmats,
        d_proximal_pmats,
        per_query,
        per_node_pmat,
        d_out);
    check_launch("PlacementLoglkKernel<PlacementRateCatGeneric>");
}



} // namespace root_likelihood
