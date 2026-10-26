#include "placement_likelihood.cuh"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "util/mlipper_util.h"

namespace mlipper::likelihood::placement {
namespace {

constexpr double kLn2 = 0.69314718055994530942;
constexpr int kMaxRateCats = 8;

__device__ __forceinline__ unsigned int placement_scaler_shift_at(
    const unsigned* scaler_pool,
    const DeviceTree& D,
    int node_id,
    size_t site_idx,
    size_t rate_idx)
{
    if (!scaler_pool || node_id < 0) return 0u;
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const int slot = clv_node_slot(D, node_id);
    if (slot < 0) return 0u;
    const size_t base = static_cast<size_t>(slot) * scaler_span(D);
    if (D.per_rate_scaling) {
        return scaler_pool[base + site_idx * rate_count + rate_idx];
    }
    return scaler_pool[base + site_idx];
}

} // namespace

template<int RATE_CATS>
__device__ __forceinline__ void load_midpoint_pmat_pair_ratecat(
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

__device__ __forceinline__ void load_midpoint_pmat_pair_states4_generic(
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

// A candidate edge is represented by the target-side upward CLV and the
// parent/sibling-side outside CLV. Transform both to the proposed attachment
// midpoint before combining them with the query contribution.
template<int RATE_CATS>
__device__ void compute_midpoint_inner_inner_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site,
    bool active_thread,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat)
{
    if (!D.d_edge_midpoint_clv) return;
    if (op.clv_pool != static_cast<uint8_t>(CLV_POOL_DOWN)) return;

    const bool target_is_left =
        op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT);
    const int target_id = target_is_left ? op.left_id : op.right_id;
    if (op.parent_id < 0 || target_id < 0) return;

    const size_t site_off =
        static_cast<size_t>(site) * static_cast<size_t>(RATE_CATS) * 4;

    fp_t* edge_midpoint =
        edge_midpoint_clv_ptr<fp_t>(D, target_id, site_off);
    const fp_t* edge_outside =
        edge_outside_clv_ptr<const fp_t>(D, target_id, site_off);
    const fp_t* target_up =
        up_clv_ptr<const fp_t>(D, target_id, site_off);
    if (!target_up) return;
    const size_t matrix_offset = (size_t)target_id * (size_t)RATE_CATS * 16;
    const fp_t* target_mat = D.d_pmat_mid_prox
        ? D.d_pmat_mid_prox + matrix_offset
        : (D.d_pmat_mid ? D.d_pmat_mid + matrix_offset : nullptr);
    const fp_t* parent_mat = D.d_pmat_mid_dist
        ? D.d_pmat_mid_dist + matrix_offset
        : (D.d_pmat_mid ? D.d_pmat_mid + matrix_offset : nullptr);
    if (!target_mat || !parent_mat || !edge_outside) return;

    load_midpoint_pmat_pair_ratecat<RATE_CATS>(
        shared_target_mat,
        shared_parent_mat,
        target_mat,
        parent_mat);
    if (!active_thread) return;

    unsigned int* edge_midpoint_scaler = edge_midpoint_scaler_ptr(D, target_id, site);
    unsigned int* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site);
    unsigned int* target_up_scaler = up_scaler_ptr(D, target_id, site);

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        unsigned int inherited_shift =
            read_clv_scaler_shift(D, edge_outside_scaler, r);
        if (target_up_scaler) {
            inherited_shift += read_clv_scaler_shift(D, target_up_scaler, r);
        }
        write_clv_scaler_shift(D, edge_midpoint_scaler, r, inherited_shift);

        const fp_t* Mtarget = shared_target_mat + (size_t)r * 16;
        const fp_t* Mparent = shared_parent_mat + (size_t)r * 16;
        const fp4_t Pup = reinterpret_cast<const fp4_t*>(
            target_up + static_cast<size_t>(r) * 4)[0];
        fp_t* Pmid = edge_midpoint + static_cast<size_t>(r) * 4;
        const fp4_t Pbase = reinterpret_cast<const fp4_t*>(
            edge_outside + static_cast<size_t>(r) * 4)[0];

        const fp_t p0 = fp_dot4(make_fp4(Mparent[0], Mparent[1], Mparent[2], Mparent[3]), Pbase) *
                        fp_dot4(make_fp4(Mtarget[0], Mtarget[1], Mtarget[2], Mtarget[3]), Pup);
        const fp_t p1 = fp_dot4(make_fp4(Mparent[4], Mparent[5], Mparent[6], Mparent[7]), Pbase) *
                        fp_dot4(make_fp4(Mtarget[4], Mtarget[5], Mtarget[6], Mtarget[7]), Pup);
        const fp_t p2 = fp_dot4(make_fp4(Mparent[8], Mparent[9], Mparent[10], Mparent[11]), Pbase) *
                        fp_dot4(make_fp4(Mtarget[8], Mtarget[9], Mtarget[10], Mtarget[11]), Pup);
        const fp_t p3 = fp_dot4(make_fp4(Mparent[12], Mparent[13], Mparent[14], Mparent[15]), Pbase) *
                        fp_dot4(make_fp4(Mtarget[12], Mtarget[13], Mtarget[14], Mtarget[15]), Pup);

        Pmid[0] = p0;
        Pmid[1] = p1;
        Pmid[2] = p2;
        Pmid[3] = p3;

        scale_clv_states4_if_needed(D, edge_midpoint_scaler, (unsigned int)r, Pmid);

    }
}

__device__ void compute_midpoint_inner_inner_states4_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site,
    bool active_thread,
    fp_t* shared_target_mat,
    fp_t* shared_parent_mat)
{
    if (!D.d_edge_midpoint_clv || D.states != 4) return;
    if (D.rate_cats <= 0 ||
        D.rate_cats > kMaxRateCats) return;
    if (op.clv_pool != static_cast<uint8_t>(CLV_POOL_DOWN)) return;

    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_id = target_is_left ? op.left_id : op.right_id;
    if (op.parent_id < 0 || target_id < 0) return;

    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t matrix_offset = static_cast<size_t>(target_id) * rate_count * 16;
    const fp_t* target_mat = D.d_pmat_mid_prox
        ? D.d_pmat_mid_prox + matrix_offset
        : (D.d_pmat_mid ? D.d_pmat_mid + matrix_offset : nullptr);
    const fp_t* parent_mat = D.d_pmat_mid_dist
        ? D.d_pmat_mid_dist + matrix_offset
        : (D.d_pmat_mid ? D.d_pmat_mid + matrix_offset : nullptr);
    if (!target_mat || !parent_mat) return;

    load_midpoint_pmat_pair_states4_generic(
        shared_target_mat,
        shared_parent_mat,
        target_mat,
        parent_mat,
        D.rate_cats);
    if (!active_thread) return;

    const size_t site_off = static_cast<size_t>(site) * rate_count * 4;
    fp_t* edge_midpoint = edge_midpoint_clv_ptr<fp_t>(D, target_id, site_off);
    const fp_t* edge_outside = edge_outside_clv_ptr<const fp_t>(D, target_id, site_off);
    const fp_t* target_up = up_clv_ptr<const fp_t>(D, target_id, site_off);
    if (!edge_midpoint || !edge_outside || !target_up) return;

    unsigned int* edge_midpoint_scaler = edge_midpoint_scaler_ptr(D, target_id, site);
    unsigned int* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site);
    unsigned int* target_up_scaler = up_scaler_ptr(D, target_id, site);

    for (int r = 0; r < D.rate_cats; ++r) {
        unsigned int inherited_shift = read_clv_scaler_shift(D, edge_outside_scaler, r);
        if (target_up_scaler) {
            inherited_shift += read_clv_scaler_shift(D, target_up_scaler, r);
        }
        write_clv_scaler_shift(D, edge_midpoint_scaler, r, inherited_shift);

        const fp_t* Mtarget = shared_target_mat + static_cast<size_t>(r) * 16;
        const fp_t* Mparent = shared_parent_mat + static_cast<size_t>(r) * 16;
        const fp4_t Pup = reinterpret_cast<const fp4_t*>(target_up + static_cast<size_t>(r) * 4)[0];
        const fp4_t Pbase = reinterpret_cast<const fp4_t*>(edge_outside + static_cast<size_t>(r) * 4)[0];
        fp_t* Pmid = edge_midpoint + static_cast<size_t>(r) * 4;

        Pmid[0] = fp_dot4(make_fp4(Mparent[0], Mparent[1], Mparent[2], Mparent[3]), Pbase) *
                  fp_dot4(make_fp4(Mtarget[0], Mtarget[1], Mtarget[2], Mtarget[3]), Pup);
        Pmid[1] = fp_dot4(make_fp4(Mparent[4], Mparent[5], Mparent[6], Mparent[7]), Pbase) *
                  fp_dot4(make_fp4(Mtarget[4], Mtarget[5], Mtarget[6], Mtarget[7]), Pup);
        Pmid[2] = fp_dot4(make_fp4(Mparent[8], Mparent[9], Mparent[10], Mparent[11]), Pbase) *
                  fp_dot4(make_fp4(Mtarget[8], Mtarget[9], Mtarget[10], Mtarget[11]), Pup);
        Pmid[3] = fp_dot4(make_fp4(Mparent[12], Mparent[13], Mparent[14], Mparent[15]), Pbase) *
                  fp_dot4(make_fp4(Mtarget[12], Mtarget[13], Mtarget[14], Mtarget[15]), Pup);

        scale_clv_states4_if_needed(D, edge_midpoint_scaler, static_cast<unsigned int>(r), Pmid);
    }
}

__global__ void UpdatePlacementMidpointPartialsKernel(
    DeviceTree D,
    const NodeOpInfo* d_ops,
    int op_offset,
    int num_ops)
{
    const int op_local = op_offset + static_cast<int>(blockIdx.y);
    const unsigned int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (!d_ops || op_local >= num_ops) return;
    const NodeOpInfo op = d_ops[op_local];
    const bool active_thread = site < D.sites;
    __shared__ fp_t shared_target_mat[kMaxRateCats * 16];
    __shared__ fp_t shared_parent_mat[kMaxRateCats * 16];
    if (D.states != 4) return;
    switch (D.rate_cats) {
        case 1:
            compute_midpoint_inner_inner_ratecat<1>(
                D, op, site, active_thread,
                shared_target_mat, shared_parent_mat);
            break;
        case 4:
            compute_midpoint_inner_inner_ratecat<4>(
                D, op, site, active_thread,
                shared_target_mat, shared_parent_mat);
            break;
        case 8:
            compute_midpoint_inner_inner_ratecat<8>(
                D, op, site, active_thread,
                shared_target_mat, shared_parent_mat);
            break;
        default:
            compute_midpoint_inner_inner_states4_generic(
                D, op, site, active_thread,
                shared_target_mat, shared_parent_mat);
            break;
    }
}

void update_midpoint_partials(
    const DeviceTree& D,
    const NodeOpInfo* d_ops,
    int num_ops,
    cudaStream_t stream)
{
    if (num_ops <= 0) return;
    if (!d_ops) {
        throw std::runtime_error("Placement midpoint update requires operations.");
    }
    if (D.states != 4 || D.rate_cats <= 0 || D.rate_cats > kMaxRateCats) {
        throw std::invalid_argument(
            "Placement midpoint update requires four states and 1-8 rate categories.");
    }
    constexpr int kMaximumGridY = 65535;
    dim3 block(256);
    const unsigned int grid_x = static_cast<unsigned int>(
        (D.sites + block.x - 1) / block.x);
    for (int op_offset = 0; op_offset < num_ops; op_offset += kMaximumGridY) {
        const int chunk_size = std::min(kMaximumGridY, num_ops - op_offset);
        dim3 grid(grid_x, static_cast<unsigned int>(chunk_size));
        UpdatePlacementMidpointPartialsKernel<<<grid, block, 0, stream>>>(
            D, d_ops, op_offset, num_ops);
        const cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            throw std::runtime_error(
                std::string("UpdatePlacementMidpointPartialsKernel: ") +
                cudaGetErrorString(error));
        }
    }
}


// Convert per-rate likelihoods and scaler exponents into one weighted site
// log likelihood. Subtracting the minimum shift keeps all rate terms on a
// common numerical scale before summation.
template<int RATE_CATS>
__device__ __forceinline__ double edge_site_loglikelihood_sum_ratecat(
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
    if (!isfinite(static_cast<double>(site_lk)) || site_lk < fp_t(0)) return NAN;
    const fp_t positive_likelihood = site_lk > fp_t(0) ? site_lk : FP_EPS;
    return static_cast<double>(site_weight) *
        (static_cast<double>(fp_log(positive_likelihood))
         - static_cast<double>(site_min_shift) * kLn2);
}

__device__ __forceinline__ double edge_site_loglikelihood_sum_generic(
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
    if (!isfinite(static_cast<double>(site_lk)) || site_lk < fp_t(0)) return NAN;
    const fp_t positive_likelihood = site_lk > fp_t(0) ? site_lk : FP_EPS;
    return static_cast<double>(site_weight) *
        (static_cast<double>(fp_log(positive_likelihood))
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
__global__ void ComputeEdgeLoglikelihoodKernel(
    DeviceTree D,
    const NodeOpInfo* __restrict__ d_ops,
    const fp_t* __restrict__ d_pendant_pmats,
    const fp_t* __restrict__ d_distal_pmats,
    const fp_t* __restrict__ d_proximal_pmats,
    size_t per_query,
    size_t per_node_pmat,
    int op_batch_offset,
    fp_t* __restrict__ d_out)
{
    const int op_local = op_batch_offset + static_cast<int>(blockIdx.y);
    if (!d_ops || op_local >= D.N) {
        if (threadIdx.x == 0) d_out[op_local] = static_cast<fp_t>(NAN);
        return;
    }
    const NodeOpInfo op = d_ops[op_local];
    const int target_id = node_op_target_id(op);
    if (target_id < 0 || target_id >= D.N) {
        if (threadIdx.x == 0) d_out[op_local] = static_cast<fp_t>(NAN);
        return;
    }

    const size_t op_offset = static_cast<size_t>(op_local);
    const size_t target_offset = static_cast<size_t>(target_id);
    const int query_idx = node_op_query_idx(op);
    if (query_idx < 0 || query_idx >= D.query_capacity) {
        if (threadIdx.x == 0) d_out[op_local] = static_cast<fp_t>(NAN);
        return;
    }
    const size_t query_offset = static_cast<size_t>(query_idx);
    const fp_t* pendant_pmat = d_pendant_pmats ? d_pendant_pmats + op_offset * per_query : nullptr;
    const fp_t* distal_pmat  = d_distal_pmats ? d_distal_pmats + target_offset * per_node_pmat : nullptr;
    const fp_t* prox_pmat    = d_proximal_pmats ? d_proximal_pmats + target_offset * per_node_pmat : nullptr;
    if (!pendant_pmat || !distal_pmat || !prox_pmat) return;

    if (!D.d_query_clv || !D.d_edge_outside_clv || !D.d_clv_up || !D.d_rate_weights || !D.d_frequencies) return;
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
        const fp_t* query_clv = D.d_query_clv + query_offset * per_node + site_offset;
        const fp_t* distal_clv = edge_outside_clv_ptr<const fp_t>(D, target_id, site_offset);
        const fp_t* prox_clv = up_clv_ptr<const fp_t>(D, target_id, site_offset);
        if (!distal_clv || !prox_clv) {
            local_sum = NAN;
            break;
        }

        fp_t rate_vals[kRateCatsStorage];
        unsigned int rate_shifts[kRateCatsStorage];
        unsigned int site_min_shift = 0u;
        bool have_positive = false;
        for (int rc = 0; rc < rate_cat_count; ++rc) {
            const size_t rate_idx = static_cast<size_t>(rc);
            const unsigned int distal_shift =
                placement_scaler_shift_at(D.d_edge_outside_scaler, D, target_id, site, rate_idx);
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
            local_sum += edge_site_loglikelihood_sum_generic(
                rate_vals,
                rate_shifts,
                rate_cat_count,
                site_min_shift,
                D.d_rate_weights,
                site_weight);
        } else {
            local_sum += edge_site_loglikelihood_sum_ratecat<kRateCatsStorage>(
                rate_vals,
                rate_shifts,
                site_min_shift,
                D.d_rate_weights,
                site_weight);
        }
    }

    const double block_sum = block_reduce_sum_double(local_sum);
    if (threadIdx.x == 0) d_out[op_local] = static_cast<fp_t>(block_sum);
}

void compute_edge_loglikelihoods(
    const DeviceTree& D,
    const NodeOpInfo* d_ops,
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
    if (!D.d_query_clv || !D.d_edge_outside_clv || !D.d_clv_up ||
        !D.d_rate_weights || !D.d_frequencies) {
        throw std::runtime_error("Placement buffers not initialized.");
    }
    if (D.states <= 0 || D.rate_cats <= 0 || D.rate_cats > kMaxRateCats) {
        throw std::invalid_argument(
            "Placement likelihood requires positive states and between 1 and 8 rate categories.");
    }

    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t per_query = rate_count * state_count * state_count;
    const size_t per_node_pmat = per_query;
    auto check_launch = [&](const char* stage) {
        const cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            throw std::runtime_error(
                std::string(stage) + ": " + cudaGetErrorString(err));
        }
    };

    constexpr int kMaximumGridY = 65535;
    dim3 block(256);
    for (int op_offset = 0; op_offset < num_ops; op_offset += kMaximumGridY) {
        const int chunk_size = std::min(kMaximumGridY, num_ops - op_offset);
        dim3 grid(1, static_cast<unsigned int>(chunk_size));
        if (D.states == 4) {
            switch (D.rate_cats) {
                case 1:
                    ComputeEdgeLoglikelihoodKernel<PlacementRateCatStates4<1>><<<grid, block, 0, stream>>>(
                        D, d_ops, d_pendant_pmats, d_distal_pmats, d_proximal_pmats,
                        per_query, per_node_pmat, op_offset, d_out);
                    check_launch("ComputeEdgeLoglikelihoodKernel<PlacementRateCatStates4<1>>");
                    continue;
                case 4:
                    ComputeEdgeLoglikelihoodKernel<PlacementRateCatStates4<4>><<<grid, block, 0, stream>>>(
                        D, d_ops, d_pendant_pmats, d_distal_pmats, d_proximal_pmats,
                        per_query, per_node_pmat, op_offset, d_out);
                    check_launch("ComputeEdgeLoglikelihoodKernel<PlacementRateCatStates4<4>>");
                    continue;
                case 8:
                    ComputeEdgeLoglikelihoodKernel<PlacementRateCatStates4<8>><<<grid, block, 0, stream>>>(
                        D, d_ops, d_pendant_pmats, d_distal_pmats, d_proximal_pmats,
                        per_query, per_node_pmat, op_offset, d_out);
                    check_launch("ComputeEdgeLoglikelihoodKernel<PlacementRateCatStates4<8>>");
                    continue;
                default:
                    break;
            }
        }
        ComputeEdgeLoglikelihoodKernel<PlacementRateCatGeneric><<<grid, block, 0, stream>>>(
            D,
            d_ops,
            d_pendant_pmats,
            d_distal_pmats,
            d_proximal_pmats,
            per_query,
            per_node_pmat,
            op_offset,
            d_out);
        check_launch("ComputeEdgeLoglikelihoodKernel<PlacementRateCatGeneric>");
    }
}

} // namespace mlipper::likelihood::placement
