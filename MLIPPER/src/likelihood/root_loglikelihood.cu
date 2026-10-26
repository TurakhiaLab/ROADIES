#include "root_loglikelihood.cuh"

#include <cuda_runtime.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "gpu/device_buffer.hpp"
#include "tree/tree.hpp"
#include "util/mlipper_util.h"

namespace mlipper::likelihood::root {

constexpr double kLn2 = 0.69314718055994530942;

__device__ __forceinline__ unsigned int root_scaler_shift_at(
    const DeviceTree& D,
    const unsigned* root_scaler,
    size_t site_idx,
    size_t rate_idx,
    size_t rate_cats)
{
    if (!root_scaler) return 0u;
    if (D.per_rate_scaling) {
        return root_scaler[site_idx * rate_cats + rate_idx];
    }
    return root_scaler[site_idx];
}

__device__ __forceinline__ double root_logaddexp(double lhs, double rhs)
{
    if (isnan(lhs) || isnan(rhs)) return NAN;
    if (lhs == INFINITY || rhs == INFINITY) return INFINITY;
    if (lhs == -INFINITY) return rhs;
    if (rhs == -INFINITY) return lhs;
    const double hi = fmax(lhs, rhs);
    const double lo = fmin(lhs, rhs);
    return hi + log1p(exp(lo - hi));
}

__device__ __forceinline__ double finalize_log_root_site_loglik(
    double variable_log,
    const fp_t* freqs,
    size_t site_idx,
    const DeviceTree& D,
    const int* invar_indices,
    double invar_proportion,
    int states)
{
    if (!isfinite(invar_proportion) || invar_proportion < 0.0 ||
        invar_proportion >= 1.0) {
        return NAN;
    }
    double loglk = variable_log;
    if (isfinite(loglk) && invar_proportion > 0.0)
        loglk += log1p(-invar_proportion);
    if (invar_indices) {
        const int inv_idx = invar_indices[site_idx];
        if (inv_idx >= states) return NAN;
        const double invariant = inv_idx >= 0
            ? invar_proportion * static_cast<double>(freqs[inv_idx]) : 0.0;
        if (invariant > 0.0) {
            const double invariant_log = log(invariant);
            loglk = root_logaddexp(loglk, invariant_log);
        }
    }
    if (isnan(loglk) || loglk == INFINITY) return NAN;
    if (loglk == -INFINITY) loglk = log(1e-300);
    return loglk * static_cast<double>(
        D.d_pattern_weights_u ? D.d_pattern_weights_u[site_idx] : 1u);
}

// Specialized device helper for common four-state/rate-category shapes.
template<int RC>
__device__ __forceinline__ double compute_root_loglikelihood_states4(
    const DeviceTree& D,
    const fp_t* clv_site,
    const unsigned* root_scaler,
    const fp_t* freqs,
    const fp_t* rate_weights,
    size_t site_idx,
    const int* invar_indices,
    double invar_proportion)
{
    const double pi0 = static_cast<double>(freqs[0]);
    const double pi1 = static_cast<double>(freqs[1]);
    const double pi2 = static_cast<double>(freqs[2]);
    const double pi3 = static_cast<double>(freqs[3]);

    double log_sum_rate = -INFINITY;
#pragma unroll
    for (int r = 0; r < RC; ++r) {
        const fp4_t a = reinterpret_cast<const fp4_t*>(clv_site)[r];
        double val = fma(static_cast<double>(a.x), pi0,
                     fma(static_cast<double>(a.y), pi1,
                     fma(static_cast<double>(a.z), pi2, static_cast<double>(a.w) * pi3)));
        const unsigned int shift = root_scaler_shift_at(
            D,
            root_scaler,
            site_idx,
            static_cast<size_t>(r),
            static_cast<size_t>(RC));
        const double weight = static_cast<double>(rate_weights[r]);
        if (val > 0.0 && weight > 0.0) {
        const double term =
            log(weight) + log(val) - static_cast<double>(shift) * kLn2;
            log_sum_rate = root_logaddexp(log_sum_rate, term);
        }
    }
    return finalize_log_root_site_loglik(
        log_sum_rate,
        freqs,
        site_idx,
        D,
        invar_indices,
        invar_proportion,
        D.states);
}

// Generic device root log-likelihood for any state/rate counts (fallback).
__device__ __forceinline__ double compute_root_loglikelihood_generic(
    const DeviceTree& D,
    const fp_t* clv_site,
    const unsigned* root_scaler,
    const fp_t* freqs,
    const fp_t* rate_weights,
    size_t site_idx,
    const int* invar_indices,
    double invar_proportion)
{
    const int rate_cats = D.rate_cats;
    const int states = D.states;
    const size_t state_count = static_cast<size_t>(states);
    const size_t rate_count = static_cast<size_t>(rate_cats);
    double log_sum_rate = -INFINITY;
    for (int r = 0; r < rate_cats; ++r) {
        const fp_t* cr = clv_site + static_cast<size_t>(r) * state_count;
        double val = 0.0;
        for (int s = 0; s < states; ++s) {
            val = fma(
                static_cast<double>(cr[s]),
                static_cast<double>(freqs[s]),
                val);
        }
        const unsigned int shift = root_scaler_shift_at(
            D, root_scaler, site_idx, static_cast<size_t>(r), rate_count);
        const double weight = static_cast<double>(rate_weights[r]);
        if (val > 0.0 && weight > 0.0) {
            const double term =
                log(weight) + log(val) - static_cast<double>(shift) * kLn2;
            log_sum_rate = root_logaddexp(log_sum_rate, term);
        }
    }
    return finalize_log_root_site_loglik(
        log_sum_rate,
        freqs,
        site_idx,
        D,
        invar_indices,
        invar_proportion,
        D.states);
}

// Explicit site indexing lets both the root kernel and other device workflows
// reuse the same per-site likelihood semantics.
__device__ double compute_root_loglikelihood_at_site(
    const DeviceTree& D,
    const fp_t* root_clv,
    const unsigned* root_scaler,
    const fp_t* freqs,
    const fp_t* rate_weights,
    const int* invar_indices,
    double invar_proportion,
    size_t site_idx)
{
    if (site_idx >= D.sites) return 0.0;

    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t per_site = rate_count * state_count;
    if (!root_clv) return 0.0;
    const fp_t* clv_site = root_clv + site_idx * per_site;
    if (D.states == 4) {
        switch (D.rate_cats) {
            case 1:
                return compute_root_loglikelihood_states4<1>(
                    D, clv_site, root_scaler, freqs, rate_weights,
                    site_idx, invar_indices, invar_proportion);
            case 4:
                return compute_root_loglikelihood_states4<4>(
                    D, clv_site, root_scaler, freqs, rate_weights,
                    site_idx, invar_indices, invar_proportion);
            case 8:
                return compute_root_loglikelihood_states4<8>(
                    D, clv_site, root_scaler, freqs, rate_weights,
                    site_idx, invar_indices, invar_proportion);
            default:
                return compute_root_loglikelihood_generic(
                    D, clv_site, root_scaler, freqs, rate_weights,
                    site_idx, invar_indices, invar_proportion);
        }
    }
    return compute_root_loglikelihood_generic(
        D, clv_site, root_scaler, freqs, rate_weights, site_idx,
        invar_indices, invar_proportion);
}

__global__ void ComputeRootLoglikelihoodKernel(
    DeviceTree D,
    const fp_t* root_clv,
    const unsigned* root_scaler,
    const int* __restrict__ d_invar_indices,
    double invar_proportion,
    size_t site_begin,
    size_t site_count,
    double* __restrict__ block_totals)
{
    double local_sum = 0.0;
    const size_t block_start =
        static_cast<size_t>(blockIdx.x) * static_cast<size_t>(blockDim.x) +
        static_cast<size_t>(threadIdx.x);
    const size_t stride =
        static_cast<size_t>(gridDim.x) * static_cast<size_t>(blockDim.x);
    const size_t site_end = site_begin + site_count;
    for (size_t site = site_begin + block_start; site < site_end; site += stride) {
        local_sum += compute_root_loglikelihood_at_site(
            D,
            root_clv,
            root_scaler,
            D.d_frequencies,
            D.d_rate_weights,
            d_invar_indices,
            invar_proportion,
            site);
    }

    const double block_sum = block_reduce_sum_double(local_sum);
    if (threadIdx.x == 0) block_totals[blockIdx.x] = block_sum;
}

__global__ void FinalizeRootLoglikelihoodKernel(
    const double* __restrict__ block_totals,
    size_t block_count,
    double* __restrict__ total)
{
    double local_sum = 0.0;
    for (size_t block = static_cast<size_t>(threadIdx.x);
         block < block_count;
         block += static_cast<size_t>(blockDim.x)) {
        local_sum += block_totals[block];
    }
    const double sum = block_reduce_sum_double(local_sum);
    if (threadIdx.x == 0) total[0] = sum;
}

double compute_root_loglikelihood(
    const DeviceTree& D,
    int root_id,
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
    if (!D.d_clv_up) {
        throw std::runtime_error("Device CLV buffers are not initialized.");
    }
    if (D.states <= 0 || D.rate_cats <= 0) {
        throw std::invalid_argument("Root likelihood requires positive state and rate-category counts.");
    }
    if (!std::isfinite(invar_proportion) || invar_proportion < 0.0 ||
        invar_proportion >= 1.0) {
        throw std::invalid_argument("Invariant-site proportion must be finite and in [0, 1).");
    }

    const size_t root_slot = static_cast<size_t>(root_id);
    const size_t per_node = D.per_node_elems();
    const fp_t* root_clv = D.d_clv_up + root_slot * per_node;
    const unsigned* root_scaler = D.d_site_scaler_up
        ? D.d_site_scaler_up + root_slot * D.scaler_elems()
        : nullptr;
    constexpr size_t site_begin = 0;
    const size_t site_count = D.sites;
    if (site_count == 0) {
        return 0.0;
    }

    dim3 block(256);
    const unsigned int grid_x = static_cast<unsigned int>(
        (site_count + block.x - 1) / block.x);
    dim3 grid(grid_x);
    struct RootReductionScratch {
        int device = -1;
        mlipper::gpu::DeviceBuffer<double> block_totals;
        mlipper::gpu::DeviceBuffer<double> total;
    };
    // Reuse reductions within one host thread, but never across CUDA devices:
    // DeviceBuffer does not carry device identity in its type.
    static thread_local RootReductionScratch reduction_scratch;
    int current_device = -1;
    CUDA_CHECK(cudaGetDevice(&current_device));
    if (reduction_scratch.device != current_device) {
        reduction_scratch.block_totals.reset();
        reduction_scratch.total.reset();
        reduction_scratch.device = current_device;
    }
    reduction_scratch.block_totals.ensureCapacity(grid_x);
    reduction_scratch.total.ensureCapacity(1);
    double total = 0.0;
    ComputeRootLoglikelihoodKernel<<<grid, block, 0, stream>>>(
        D,
        root_clv,
        root_scaler,
        d_invar_indices,
        invar_proportion,
        site_begin,
        site_count,
        reduction_scratch.block_totals.get());
    CUDA_CHECK(cudaGetLastError());
    FinalizeRootLoglikelihoodKernel<<<1, block, 0, stream>>>(
        reduction_scratch.block_totals.get(), grid_x,
        reduction_scratch.total.get());
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(
        &total,
        reduction_scratch.total.get(),
        sizeof(double),
        cudaMemcpyDeviceToHost,
        stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    if (!std::isfinite(total)) {
        throw std::runtime_error("Root likelihood evaluation produced a non-finite value.");
    }
    return total;
}

} // namespace mlipper::likelihood::root
