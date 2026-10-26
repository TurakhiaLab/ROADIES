#include <cuda_runtime.h>
#include <cmath>
#include "util/mlipper_util.h"
#include "partials.cuh"

namespace mlipper::likelihood::partials {

// CLVs are laid out as [node][site][rate][state]. Upward values summarize a
// node's descendants; downward values summarize the complementary side of the
// tree. Scaling exponents follow the same site/rate policy as their CLV and
// must be inherited whenever two messages are multiplied. The specialized DNA
// paths below preserve this contract while using fp4 vector operations.

__device__ __forceinline__ unsigned int* scaler_ptr_for_pool(
    const DeviceTree& D,
    uint8_t clv_pool,
    int node_id,
    unsigned int site)
{
    if (clv_pool == static_cast<uint8_t>(CLV_POOL_DOWN)) {
        return down_scaler_ptr(D, node_id, site);
    }
    return up_scaler_ptr(D, node_id, site);
}

__device__ __forceinline__ fp4_t tip_clv_states4(unsigned int mask)
{
    return make_fp4(
        (mask & 1u) ? fp_t(1) : fp_t(0),
        (mask & 2u) ? fp_t(1) : fp_t(0),
        (mask & 4u) ? fp_t(1) : fp_t(0),
        (mask & 8u) ? fp_t(1) : fp_t(0));
}

__device__ __forceinline__ void store_fp4(fp_t* dst, const fp4_t& value)
{
    dst[0] = value.x;
    dst[1] = value.y;
    dst[2] = value.z;
    dst[3] = value.w;
}

__device__ __forceinline__ fp_t masked_tip_sum_states4(
    const fp_t* row,
    unsigned int mask)
{
    fp_t sum = fp_t(0);
    if (mask & 1u) sum += row[0];
    if (mask & 2u) sum += row[1];
    if (mask & 4u) sum += row[2];
    if (mask & 8u) sum += row[3];
    return sum;
}

__device__ __forceinline__ void write_downward_inherited_scalers_states4(
    const DeviceTree& D,
    unsigned int* parent_scaler,
    unsigned int* sibling_scaler,
    unsigned int* target_up_scaler,
    unsigned int* down_scaler,
    unsigned int* edge_midpoint_scaler,
    unsigned int* edge_outside_scaler,
    unsigned int rate_idx,
    bool write_edge_outside)
{
    const unsigned int down_inherited =
        read_clv_scaler_shift(D, parent_scaler, rate_idx) +
        read_clv_scaler_shift(D, sibling_scaler, rate_idx);
    const unsigned int mid_inherited =
        down_inherited + read_clv_scaler_shift(D, target_up_scaler, rate_idx);
    write_clv_scaler_shift(D, down_scaler, rate_idx, down_inherited);
    write_clv_scaler_shift(D, edge_midpoint_scaler, rate_idx, mid_inherited);
    if (write_edge_outside) {
        write_clv_scaler_shift(D, edge_outside_scaler, rate_idx, down_inherited);
    }
}

__device__ __forceinline__ void build_midpoint_states4(
    const fp_t* half_mat,
    fp_t p0, fp_t p1, fp_t p2, fp_t p3,
    const fp4_t& target_up,
    fp_t* out_mid)
{
    const fp4_t parent_vec = make_fp4(p0, p1, p2, p3);
    out_mid[0] = fp_dot4(make_fp4(half_mat[0], half_mat[4], half_mat[8],  half_mat[12]), parent_vec) *
                 fp_dot4(make_fp4(half_mat[0], half_mat[4], half_mat[8],  half_mat[12]), target_up);
    out_mid[1] = fp_dot4(make_fp4(half_mat[1], half_mat[5], half_mat[9],  half_mat[13]), parent_vec) *
                 fp_dot4(make_fp4(half_mat[1], half_mat[5], half_mat[9],  half_mat[13]), target_up);
    out_mid[2] = fp_dot4(make_fp4(half_mat[2], half_mat[6], half_mat[10], half_mat[14]), parent_vec) *
                 fp_dot4(make_fp4(half_mat[2], half_mat[6], half_mat[10], half_mat[14]), target_up);
    out_mid[3] = fp_dot4(make_fp4(half_mat[3], half_mat[7], half_mat[11], half_mat[15]), parent_vec) *
                 fp_dot4(make_fp4(half_mat[3], half_mat[7], half_mat[11], half_mat[15]), target_up);
}

__device__ __forceinline__ fp_t compute_tip_tip_states4_rate(
    const fp_t* left_mat,
    const fp_t* right_mat,
    unsigned int left_mask,
    unsigned int right_mask,
    fp_t* out)
{
    fp_t max_val = fp_t(0);
    const fp_t* left_row = left_mat;
    const fp_t* right_row = right_mat;
#pragma unroll
    for (int i = 0; i < 4; ++i) {
        const fp_t left_term = masked_tip_sum_states4(left_row, left_mask);
        const fp_t right_term = masked_tip_sum_states4(right_row, right_mask);
        const fp_t value = left_term * right_term;
        out[i] = value;
        if (value > max_val) max_val = value;
        left_row += 4;
        right_row += 4;
    }
    return max_val;
}

__device__ __forceinline__ fp_t compute_tip_inner_states4_rate(
    const fp_t* tip_mat,
    const fp_t* inner_mat,
    const fp_t* inner_clv,
    unsigned int tip_mask,
    fp_t* out)
{
    const fp4_t inner = make_fp4(inner_clv[0], inner_clv[1], inner_clv[2], inner_clv[3]);
    fp_t max_val = fp_t(0);
    const fp_t* tip_row = tip_mat;
    const fp_t* inner_row = inner_mat;
#pragma unroll
    for (int i = 0; i < 4; ++i) {
        const fp_t left_term = masked_tip_sum_states4(tip_row, tip_mask);
        const fp_t right_term = fp_dot4(
            make_fp4(inner_row[0], inner_row[1], inner_row[2], inner_row[3]),
            inner);
        const fp_t value = left_term * right_term;
        out[i] = value;
        if (value > max_val) max_val = value;
        tip_row += 4;
        inner_row += 4;
    }
    return max_val;
}

__device__ __forceinline__ fp_t compute_inner_inner_states4_rate(
    const fp_t* left_mat,
    const fp_t* right_mat,
    const fp_t* left_clv,
    const fp_t* right_clv,
    fp_t* out)
{
    const fp4_t left = make_fp4(left_clv[0], left_clv[1], left_clv[2], left_clv[3]);
    const fp4_t right = make_fp4(right_clv[0], right_clv[1], right_clv[2], right_clv[3]);
    fp_t max_val = fp_t(0);
    const fp_t* left_row = left_mat;
    const fp_t* right_row = right_mat;
#pragma unroll
    for (int i = 0; i < 4; ++i) {
        const fp_t left_term = fp_dot4(
            make_fp4(left_row[0], left_row[1], left_row[2], left_row[3]),
            left);
        const fp_t right_term = fp_dot4(
            make_fp4(right_row[0], right_row[1], right_row[2], right_row[3]),
            right);
        const fp_t value = left_term * right_term;
        out[i] = value;
        if (value > max_val) max_val = value;
        left_row += 4;
        right_row += 4;
    }
    return max_val;
}

__global__ void InitializeTipPartialsKernel(const DeviceTree D)
{
    if (!D.d_tipchars || !D.d_tip_node_ids || !D.d_clv_up) return;

    const size_t tip_site_idx =
        static_cast<size_t>(blockIdx.x) * static_cast<size_t>(blockDim.x) +
        static_cast<size_t>(threadIdx.x);
    const size_t total_tip_sites = static_cast<size_t>(D.tips) * D.sites;
    if (tip_site_idx >= total_tip_sites) return;

    const size_t tip_idx = tip_site_idx / D.sites;
    const size_t site = tip_site_idx % D.sites;
    const int node_id = D.d_tip_node_ids[tip_idx];
    if (node_id < 0 || node_id >= D.capacity_N) return;

    const unsigned int mask = D.d_tipmap[D.d_tipchars[tip_site_idx]];
    const size_t per_node = per_node_span(D);
    const size_t site_off = site * static_cast<size_t>(D.rate_cats) * static_cast<size_t>(D.states);
    fp_t* tip_up = up_clv_ptr<fp_t>(D, node_id, site_off);
    unsigned int* tip_scaler = up_scaler_ptr(D, node_id, site);
    if (!tip_up) return;

    if (D.states == 4) {
        const fp4_t tip = tip_clv_states4(mask);
        for (int r = 0; r < D.rate_cats; ++r) {
            if (tip_scaler) {
                write_clv_scaler_shift(D, tip_scaler, r, 0u);
            }
            store_fp4(tip_up + static_cast<size_t>(r) * 4, tip);
        }
        return;
    }

    for (int r = 0; r < D.rate_cats; ++r) {
        if (tip_scaler) {
            write_clv_scaler_shift(D, tip_scaler, r, 0u);
        }
        fp_t* out = tip_up + static_cast<size_t>(r) * static_cast<size_t>(D.states);
        for (int s = 0; s < D.states; ++s) {
            out[s] = (mask & (1u << s)) ? fp_t(1) : fp_t(0);
        }
    }
}

// Downward operations combine the parent's outside contribution with the
// upward contribution from the target's sibling. Separate cases avoid
// materializing tip CLVs when either side is a directly encoded tip.
__device__ __forceinline__ void compute_downward_inner_inner_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    if (op.clv_pool != static_cast<uint8_t>(CLV_POOL_DOWN)) return;
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_id  = target_is_left ? op.left_id  : op.right_id;
    const int sibling_id = target_is_left ? op.right_id : op.left_id;
    if (target_id < 0 || sibling_id < 0) return;

    const unsigned int states    = (unsigned int)D.states;
    const unsigned int rate_cats = (unsigned int)D.rate_cats;
    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * (size_t)states * (size_t)rate_cats;

    const fp_t* parent_down = down_clv_ptr<const fp_t>(D, op.parent_id, site_off);
    const fp_t* sibling_up  = up_clv_ptr<const fp_t>(D, sibling_id, site_off);
    fp_t*       target_down = down_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_midpoint  = edge_midpoint_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_outside    = edge_outside_clv_ptr<fp_t>(D, target_id, site_off);
    const fp_t* target_up   = up_clv_ptr<const fp_t>(D, target_id, site_off);
    if (!parent_down || !sibling_up || !target_down) return;

    const fp_t* target_mat  = D.d_pmat + (size_t)target_id  * rate_cats * states * states;
    const fp_t* target_mat_half = D.d_pmat_mid
        ? (D.d_pmat_mid + (size_t)target_id * rate_cats * states * states)
        : target_mat;
    const fp_t* sibling_mat = D.d_pmat + (size_t)sibling_id * rate_cats * states * states;
    unsigned int* parent_scaler = down_scaler_ptr(D, op.parent_id, site);
    unsigned int* sibling_scaler = up_scaler_ptr(D, sibling_id, site);
    unsigned int* target_up_scaler = up_scaler_ptr(D, target_id, site);
    unsigned int* down_scaler = down_scaler_ptr(D, target_id, site);
    unsigned int* edge_midpoint_scaler = edge_midpoint_scaler_ptr(D, target_id, site);
    unsigned int* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site);

    for (unsigned int r = 0; r < rate_cats; ++r) {
        const unsigned int parent_shift = read_clv_scaler_shift(D, parent_scaler, r);
        const unsigned int sibling_shift = read_clv_scaler_shift(D, sibling_scaler, r);
        const unsigned int target_up_shift = read_clv_scaler_shift(D, target_up_scaler, r);
        const unsigned int down_inherited = parent_shift + sibling_shift;
        const unsigned int mid_inherited = down_inherited + target_up_shift;
        write_clv_scaler_shift(D, down_scaler, r, down_inherited);
        write_clv_scaler_shift(D, edge_midpoint_scaler, r, mid_inherited);
        if (edge_outside) write_clv_scaler_shift(D, edge_outside_scaler, r, down_inherited);

        const fp_t* Tmat = target_mat  + (size_t)r * states * states;
        const fp_t* Thalf= target_mat_half + (size_t)r * states * states;
        const fp_t* Smat = sibling_mat + (size_t)r * states * states;
        const fp_t* Ppar = parent_down + (size_t)r * states;
        const fp_t* Psib = sibling_up  + (size_t)r * states;
        fp_t*       Pout = target_down + (size_t)r * states;
        fp_t*       Pmid  = (edge_midpoint && target_up) ? (edge_midpoint + (size_t)r * states) : nullptr;
        const fp_t* Pup   = target_up ? (target_up + (size_t)r * states) : nullptr;
        fp_t*       Pbase = edge_outside ? (edge_outside + (size_t)r * states) : nullptr;

        double sib_to_parent[64];
        for (unsigned int j = 0; j < states; ++j) {
            const fp_t* row = Smat + j * states;
            double acc = 0.0;
            for (unsigned int k = 0; k < states; ++k) acc += row[k] * Psib[k];
            sib_to_parent[j] = acc;
        }

        double col_scale_max_val = 0.0;
        for (unsigned int i = 0; i < states; ++i) {
            double acc = 0.0;
            for (unsigned int j = 0; j < states; ++j) {
                const size_t matrix_index =
                    D.downward_pmat_indexing == DownwardPmatIndexing::Rows
                    ? static_cast<size_t>(i) * states + j
                    : static_cast<size_t>(j) * states + i;
                acc += Tmat[matrix_index] *
                    (Ppar[j] * sib_to_parent[j]);
            }
            Pout[i] = acc;
            if (acc > col_scale_max_val) col_scale_max_val = acc;
        }

        double pbase_max_val = 0.0;
        double pmid_max_val = 0.0;
        if (Pmid) {
            // Cache parent_down * sibling_up (after sibling branch matrix) per state.
            if (Pbase) {
                for (unsigned int j = 0; j < states; ++j) {
                    Pbase[j] = Ppar[j] * sib_to_parent[j];
                    if (Pbase[j] > pbase_max_val) pbase_max_val = Pbase[j];
                }
            }
            for (unsigned int i = 0; i < states; ++i) {
                const fp_t* Throw = Thalf + i * states;
                double par_acc = 0.0, tgt_acc = 0.0;
                for (unsigned int j = 0; j < states; ++j) {
                    const double pj = Pbase ? Pbase[j] : (Ppar[j] * sib_to_parent[j]);
                    par_acc += Throw[j] * pj;
                    tgt_acc += Throw[j] * Pup[j];
                }
                const double val = par_acc * tgt_acc;
                Pmid[i] = val;
                if (val > pmid_max_val) pmid_max_val = val;
            }
        }

        {
            const unsigned int down_shift_local = clv_scale_shift(col_scale_max_val);
            if (down_shift_local) {
            add_clv_scaler_shift(D, down_scaler, r, down_shift_local);
            for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pout[j], down_shift_local);
            }
            const unsigned int base_shift_local = Pbase ? clv_scale_shift(pbase_max_val) : 0u;
            if (base_shift_local) {
            add_clv_scaler_shift(D, edge_outside_scaler, r, base_shift_local);
            if (Pbase) {
                for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pbase[j], base_shift_local);
            }
            }
            const unsigned int mid_shift_local = Pmid ? clv_scale_shift(pmid_max_val) : 0u;
            if (mid_shift_local) {
            add_clv_scaler_shift(D, edge_midpoint_scaler, r, mid_shift_local);
            if (Pmid) {
                for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pmid[j], mid_shift_local);
            }
            }
        }
    }
}

__device__ __forceinline__ void compute_downward_inner_tip_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    if (op.clv_pool != static_cast<uint8_t>(CLV_POOL_DOWN)) return;
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_id       = target_is_left ? op.left_id  : op.right_id;
    const int sibling_tip_idx = target_is_left ? op.right_tip_index : op.left_tip_index;
    if (target_id < 0 || sibling_tip_idx < 0) return;

    const unsigned int states    = (unsigned int)D.states;
    const unsigned int rate_cats = (unsigned int)D.rate_cats;
    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * (size_t)states * (size_t)rate_cats;

    const fp_t* parent_down = down_clv_ptr<const fp_t>(D, op.parent_id, site_off);
    fp_t*       target_down = down_clv_ptr<fp_t>(D, target_id, site_off);
    const fp_t* target_up   = up_clv_ptr<const fp_t>(D, target_id, site_off);
    if (!parent_down || !target_down) return;

    const fp_t* target_mat  = D.d_pmat + (size_t)target_id * rate_cats * states * states;
    const fp_t* target_mat_half = D.d_pmat_mid
        ? (D.d_pmat_mid + (size_t)target_id * rate_cats * states * states)
        : target_mat;
    const size_t sibling_node_id = static_cast<size_t>(
        target_is_left ? op.right_id : op.left_id);
    const fp_t* sibling_mat =
        D.d_pmat + sibling_node_id * rate_cats * states * states;
    fp_t*         edge_midpoint  = edge_midpoint_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*         edge_outside    = edge_outside_clv_ptr<fp_t>(D, target_id, site_off);
    unsigned int* parent_scaler = down_scaler_ptr(D, op.parent_id, site);
    unsigned int* sibling_scaler = up_scaler_ptr(D, target_is_left ? op.right_id : op.left_id, site);
    unsigned int* target_up_scaler = up_scaler_ptr(D, target_id, site);
    unsigned int* down_scaler = down_scaler_ptr(D, target_id, site);
    unsigned int* edge_midpoint_scaler = edge_midpoint_scaler_ptr(D, target_id, site);
    unsigned int* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site);

    const unsigned char* tipchars = D.d_tipchars + (size_t)sibling_tip_idx * D.sites;

    for (unsigned int r = 0; r < rate_cats; ++r) {
        const unsigned int down_inherited =
            read_clv_scaler_shift(D, parent_scaler, r) +
            read_clv_scaler_shift(D, sibling_scaler, r);
        const unsigned int mid_inherited =
            down_inherited + read_clv_scaler_shift(D, target_up_scaler, r);
        write_clv_scaler_shift(D, down_scaler, r, down_inherited);
        write_clv_scaler_shift(D, edge_midpoint_scaler, r, mid_inherited);
        if (edge_outside) write_clv_scaler_shift(D, edge_outside_scaler, r, down_inherited);

        const unsigned int mask = D.d_tipmap[tipchars[site]];
        const fp_t* Tmat = target_mat  + (size_t)r * states * states;
        const fp_t* Thalf= target_mat_half + (size_t)r * states * states;
        const fp_t* Smat = sibling_mat + (size_t)r * states * states;
        const fp_t* Ppar = parent_down + (size_t)r * states;
        const fp_t* Pup  = target_up ? (target_up + (size_t)r * states) : nullptr;
        fp_t*       Pout = target_down + (size_t)r * states;
        fp_t*       Pmid = (target_up && D.d_edge_midpoint_clv)
            ? (edge_midpoint ? (edge_midpoint + (size_t)r * states) : nullptr)
            : nullptr;
        fp_t*       Pbase = edge_outside ? (edge_outside + (size_t)r * states) : nullptr;

        double sib_to_parent[64];
        for (unsigned int j = 0; j < states; ++j) {
            const fp_t* row = Smat + j * states;
            double acc = 0.0;
            for (unsigned int k = 0; k < states; ++k)
                if (mask & (1u << k)) acc += row[k];
            sib_to_parent[j] = acc;
        }

        double col_scale_max_val = 0.0;
        for (unsigned int i = 0; i < states; ++i) {
            double acc = 0.0;
            for (unsigned int j = 0; j < states; ++j) {
                const size_t matrix_index =
                    D.downward_pmat_indexing == DownwardPmatIndexing::Rows
                    ? static_cast<size_t>(i) * states + j
                    : static_cast<size_t>(j) * states + i;
                acc += Tmat[matrix_index] *
                    (Ppar[j] * sib_to_parent[j]);
            }
            Pout[i] = acc;
            if (acc > col_scale_max_val) col_scale_max_val = acc;
        }

        double pbase_max_val = 0.0;
        double pmid_max_val = 0.0;
        if (Pmid) {
            if (Pbase) {
                for (unsigned int j = 0; j < states; ++j) {
                    Pbase[j] = Ppar[j] * sib_to_parent[j];
                    if (Pbase[j] > pbase_max_val) pbase_max_val = Pbase[j];
                }
            }
            for (unsigned int i = 0; i < states; ++i) {
                const fp_t* Throw = Thalf + i * states;
                double par_acc = 0.0, tgt_acc = 0.0;
                for (unsigned int j = 0; j < states; ++j) {
                    const double pj = Pbase ? Pbase[j] : (Ppar[j] * sib_to_parent[j]);
                    par_acc += Throw[j] * pj;
                    tgt_acc += Throw[j] * Pup[j];
                }
                Pmid[i] = par_acc * tgt_acc;
                if (Pmid[i] > pmid_max_val) pmid_max_val = Pmid[i];
            }
        }

        {
            const unsigned int down_shift_local = clv_scale_shift(col_scale_max_val);
            if (down_shift_local) {
            add_clv_scaler_shift(D, down_scaler, r, down_shift_local);
            for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pout[j], down_shift_local);
            }
            const unsigned int base_shift_local = Pbase ? clv_scale_shift(pbase_max_val) : 0u;
            if (base_shift_local) {
            add_clv_scaler_shift(D, edge_outside_scaler, r, base_shift_local);
            if (Pbase) {
                for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pbase[j], base_shift_local);
            }
            }
            const unsigned int mid_shift_local = Pmid ? clv_scale_shift(pmid_max_val) : 0u;
            if (mid_shift_local) {
            add_clv_scaler_shift(D, edge_midpoint_scaler, r, mid_shift_local);
            if (Pmid) {
                for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pmid[j], mid_shift_local);
            }
            }
        }
    }
}

__device__ __forceinline__ void compute_downward_tip_tip_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    if (op.clv_pool != static_cast<uint8_t>(CLV_POOL_DOWN)) return;
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_tip_idx  = target_is_left ? op.left_tip_index : op.right_tip_index;
    const int sibling_tip_idx = target_is_left ? op.right_tip_index : op.left_tip_index;
    const int target_id       = target_is_left ? op.left_id : op.right_id;
    if (target_tip_idx < 0 || sibling_tip_idx < 0 || target_id < 0) return;

    const unsigned int states    = (unsigned int)D.states;
    const unsigned int rate_cats = (unsigned int)D.rate_cats;
    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * (size_t)states * (size_t)rate_cats;

    const fp_t* parent_down = down_clv_ptr<const fp_t>(D, op.parent_id, site_off);
    fp_t*       target_down = down_clv_ptr<fp_t>(D, target_id, site_off);
    const fp_t* target_up   = up_clv_ptr<const fp_t>(D, target_id, site_off);
    if (!parent_down || !target_down) return;

    const fp_t* target_mat  = D.d_pmat + (size_t)target_id * rate_cats * states * states;
    const fp_t* target_mat_half = D.d_pmat_mid
        ? (D.d_pmat_mid + (size_t)target_id * rate_cats * states * states)
        : target_mat;
    const size_t sibling_node_id = static_cast<size_t>(
        target_is_left ? op.right_id : op.left_id);
    const fp_t* sibling_mat =
        D.d_pmat + sibling_node_id * rate_cats * states * states;
    fp_t*       edge_midpoint  = edge_midpoint_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_outside    = edge_outside_clv_ptr<fp_t>(D, target_id, site_off);
    unsigned int* parent_scaler = down_scaler_ptr(D, op.parent_id, site);
    unsigned int* sibling_scaler = up_scaler_ptr(D, target_is_left ? op.right_id : op.left_id, site);
    unsigned int* target_up_scaler = up_scaler_ptr(D, target_id, site);
    unsigned int* down_scaler = down_scaler_ptr(D, target_id, site);
    unsigned int* edge_midpoint_scaler = edge_midpoint_scaler_ptr(D, target_id, site);
    unsigned int* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site);
    const unsigned char* tipchars = D.d_tipchars + (size_t)sibling_tip_idx * D.sites;

    for (unsigned int r = 0; r < rate_cats; ++r) {
        const unsigned int down_inherited =
            read_clv_scaler_shift(D, parent_scaler, r) +
            read_clv_scaler_shift(D, sibling_scaler, r);
        const unsigned int mid_inherited =
            down_inherited + read_clv_scaler_shift(D, target_up_scaler, r);
        write_clv_scaler_shift(D, down_scaler, r, down_inherited);
        write_clv_scaler_shift(D, edge_midpoint_scaler, r, mid_inherited);
        if (edge_outside) write_clv_scaler_shift(D, edge_outside_scaler, r, down_inherited);

        const unsigned int mask = D.d_tipmap[tipchars[site]];
        const fp_t* Tmat = target_mat  + (size_t)r * states * states;
        const fp_t* Thalf= target_mat_half + (size_t)r * states * states;
        const fp_t* Smat = sibling_mat + (size_t)r * states * states;
        const fp_t* Ppar = parent_down + (size_t)r * states;
        const fp_t* Pup  = target_up ? (target_up + (size_t)r * states) : nullptr;
        fp_t*       Pout = target_down + (size_t)r * states;
        fp_t*       Pmid = (target_up && edge_midpoint) ? (edge_midpoint + (size_t)r * states) : nullptr;
        fp_t*       Pbase = edge_outside ? (edge_outside + (size_t)r * states) : nullptr;

        double sib_to_parent[64];
        for (unsigned int j = 0; j < states; ++j) {
            const fp_t* row = Smat + j * states;
            double acc = 0.0;
            for (unsigned int k = 0; k < states; ++k)
                if (mask & (1u << k)) acc += row[k];
            sib_to_parent[j] = acc;
        }

        double col_scale_max_val = 0.0;
        for (unsigned int i = 0; i < states; ++i) {
            double acc = 0.0;
            for (unsigned int j = 0; j < states; ++j) {
                const size_t matrix_index =
                    D.downward_pmat_indexing == DownwardPmatIndexing::Rows
                    ? static_cast<size_t>(i) * states + j
                    : static_cast<size_t>(j) * states + i;
                acc += Tmat[matrix_index] *
                    (Ppar[j] * sib_to_parent[j]);
            }
            Pout[i] = acc;
            if (acc > col_scale_max_val) col_scale_max_val = acc;
        }

        double pbase_max_val = 0.0;
        double pmid_max_val = 0.0;
        if (Pmid) {
            if (Pbase) {
                for (unsigned int j = 0; j < states; ++j) {
                    Pbase[j] = Ppar[j] * sib_to_parent[j];
                    if (Pbase[j] > pbase_max_val) pbase_max_val = Pbase[j];
                }
            }
            for (unsigned int i = 0; i < states; ++i) {
                const fp_t* Throw = Thalf + i * states;
                double par_acc = 0.0, tgt_acc = 0.0;
                for (unsigned int j = 0; j < states; ++j) {
                    const double pj = Pbase ? Pbase[j] : (Ppar[j] * sib_to_parent[j]);
                    par_acc += Throw[j] * pj;
                    tgt_acc += Throw[j] * Pup[j];
                }
                Pmid[i] = par_acc * tgt_acc;
                if (Pmid[i] > pmid_max_val) pmid_max_val = Pmid[i];
            }
        }

        {
            const unsigned int down_shift_local = clv_scale_shift(col_scale_max_val);
            if (down_shift_local) {
            add_clv_scaler_shift(D, down_scaler, r, down_shift_local);
            for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pout[j], down_shift_local);
            }
            const unsigned int base_shift_local = Pbase ? clv_scale_shift(pbase_max_val) : 0u;
            if (base_shift_local) {
            add_clv_scaler_shift(D, edge_outside_scaler, r, base_shift_local);
            if (Pbase) {
                for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pbase[j], base_shift_local);
            }
            }
            const unsigned int mid_shift_local = Pmid ? clv_scale_shift(pmid_max_val) : 0u;
            if (mid_shift_local) {
            add_clv_scaler_shift(D, edge_midpoint_scaler, r, mid_shift_local);
            if (Pmid) {
                for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pmid[j], mid_shift_local);
            }
            }
        }
    }
}

__device__ __forceinline__ void compute_downward_tip_inner_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    if (op.clv_pool != static_cast<uint8_t>(CLV_POOL_DOWN)) return;
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_tip_idx  = target_is_left ? op.left_tip_index : op.right_tip_index;
    const int sibling_id      = target_is_left ? op.right_id : op.left_id;
    if (target_tip_idx < 0 || sibling_id < 0) return;

    const unsigned int states    = (unsigned int)D.states;
    const unsigned int rate_cats = (unsigned int)D.rate_cats;
    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * (size_t)states * (size_t)rate_cats;

    const fp_t* parent_down = down_clv_ptr<const fp_t>(D, op.parent_id, site_off);
    const fp_t* sibling_up  = up_clv_ptr<const fp_t>(D, sibling_id, site_off);
    fp_t*       target_down = down_clv_ptr<fp_t>(D, target_is_left ? op.left_id : op.right_id, site_off);
    const fp_t* target_up   = up_clv_ptr<const fp_t>(D, target_is_left ? op.left_id : op.right_id, site_off);
    if (!parent_down || !sibling_up || !target_down) return;

    const size_t target_node_id = static_cast<size_t>(
        target_is_left ? op.left_id : op.right_id);
    const fp_t* target_mat =
        D.d_pmat + target_node_id * rate_cats * states * states;
    const fp_t* target_mat_half = D.d_pmat_mid
        ? (D.d_pmat_mid + (size_t)(target_is_left ? op.left_id : op.right_id) * rate_cats * states * states)
        : target_mat;
    const fp_t* sibling_mat = D.d_pmat + (size_t)sibling_id * rate_cats * states * states;
    const int target_id       = target_is_left ? op.left_id : op.right_id;
    fp_t*       edge_midpoint  = edge_midpoint_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_outside    = edge_outside_clv_ptr<fp_t>(D, target_id, site_off);
    unsigned int* parent_scaler = down_scaler_ptr(D, op.parent_id, site);
    unsigned int* sibling_scaler = up_scaler_ptr(D, sibling_id, site);
    unsigned int* target_up_scaler = up_scaler_ptr(D, target_id, site);
    unsigned int* down_scaler = down_scaler_ptr(D, target_id, site);
    unsigned int* edge_midpoint_scaler = edge_midpoint_scaler_ptr(D, target_id, site);
    unsigned int* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site);

    const unsigned char* tipchars = D.d_tipchars + (size_t)target_tip_idx * D.sites;
    const unsigned int tmask = D.d_tipmap[tipchars[site]];

    for (unsigned int r = 0; r < rate_cats; ++r) {
        const unsigned int down_inherited =
            read_clv_scaler_shift(D, parent_scaler, r) +
            read_clv_scaler_shift(D, sibling_scaler, r);
        const unsigned int mid_inherited =
            down_inherited + read_clv_scaler_shift(D, target_up_scaler, r);
        write_clv_scaler_shift(D, down_scaler, r, down_inherited);
        write_clv_scaler_shift(D, edge_midpoint_scaler, r, mid_inherited);
        if (edge_outside) write_clv_scaler_shift(D, edge_outside_scaler, r, down_inherited);

        const fp_t* Tmat = target_mat  + (size_t)r * states * states;
        const fp_t* Thalf= target_mat_half + (size_t)r * states * states;
        const fp_t* Smat = sibling_mat + (size_t)r * states * states;
        const fp_t* Ppar = parent_down + (size_t)r * states;
        const fp_t* Psib = sibling_up  + (size_t)r * states;
        const fp_t* Pup  = target_up ? (target_up + (size_t)r * states) : nullptr;
        fp_t*       Pout = target_down + (size_t)r * states;
        fp_t*       Pmid = (target_up && edge_midpoint) ? (edge_midpoint + (size_t)r * states) : nullptr;
        fp_t*       Pbase = edge_outside ? (edge_outside + (size_t)r * states) : nullptr;

        double sib_to_parent[64];
        for (unsigned int j = 0; j < states; ++j) {
            const fp_t* row = Smat + j * states;
            double acc = 0.0;
            for (unsigned int k = 0; k < states; ++k) acc += row[k] * Psib[k];
            sib_to_parent[j] = acc;
        }

        double col_scale_max_val = 0.0;
        for (unsigned int i = 0; i < states; ++i) {
            double acc = 0.0;
            for (unsigned int j = 0; j < states; ++j) {
                const size_t matrix_index =
                    D.downward_pmat_indexing == DownwardPmatIndexing::Rows
                    ? static_cast<size_t>(i) * states + j
                    : static_cast<size_t>(j) * states + i;
                acc += Tmat[matrix_index] *
                    (Ppar[j] * sib_to_parent[j]);
            }
            Pout[i] = (tmask & (1u << i)) ? acc : 0.0;
            if (Pout[i] > col_scale_max_val) col_scale_max_val = Pout[i];
        }

        double pbase_max_val = 0.0;
        double pmid_max_val = 0.0;
        if (Pmid) {
            if (Pbase) {
                for (unsigned int j = 0; j < states; ++j) {
                    Pbase[j] = Ppar[j] * sib_to_parent[j];
                    if (Pbase[j] > pbase_max_val) pbase_max_val = Pbase[j];
                }
            }
            for (unsigned int i = 0; i < states; ++i) {
                const fp_t* Throw = Thalf + i * states;
                double par_acc = 0.0, tgt_acc = 0.0;
                for (unsigned int j = 0; j < states; ++j) {
                    const double pj = Pbase ? Pbase[j] : (Ppar[j] * sib_to_parent[j]);
                    par_acc += Throw[j] * pj;
                    tgt_acc += Throw[j] * Pup[j];
                }
                Pmid[i] = (tmask & (1u << i)) ? (par_acc * tgt_acc) : 0.0;
                if (Pmid[i] > pmid_max_val) pmid_max_val = Pmid[i];
            }
        }

        {
            const unsigned int down_shift_local = clv_scale_shift(col_scale_max_val);
            if (down_shift_local) {
            add_clv_scaler_shift(D, down_scaler, r, down_shift_local);
            for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pout[j], down_shift_local);
            }
            const unsigned int base_shift_local = Pbase ? clv_scale_shift(pbase_max_val) : 0u;
            if (base_shift_local) {
            add_clv_scaler_shift(D, edge_outside_scaler, r, base_shift_local);
            if (Pbase) {
                for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pbase[j], base_shift_local);
            }
            }
            const unsigned int mid_shift_local = Pmid ? clv_scale_shift(pmid_max_val) : 0u;
            if (mid_shift_local) {
            add_clv_scaler_shift(D, edge_midpoint_scaler, r, mid_shift_local);
            if (Pmid) {
                for (unsigned int j = 0; j < states; ++j)
                    fp_scale_pow2(Pmid[j], mid_shift_local);
            }
            }
        }
    }
}

// Four-state specializations keep rate-category counts compile-time constants
// so CUDA can unroll the small matrix operations.
template<int RATE_CATS>
__device__ __forceinline__ void compute_downward_inner_inner_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    if (op.clv_pool != static_cast<uint8_t>(CLV_POOL_DOWN)) return;
    if (op.left_tip_index >= 0 || op.right_tip_index >= 0) return;

    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_id  = target_is_left ? op.left_id  : op.right_id;
    const int sibling_id = target_is_left ? op.right_id : op.left_id;
    if (target_id < 0 || sibling_id < 0) return;

    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * (size_t)RATE_CATS * 4;

    const fp_t* parent_down = down_clv_ptr<const fp_t>(D, op.parent_id, site_off);
    const fp_t* sibling_up  = up_clv_ptr<const fp_t>(D, sibling_id, site_off);
    fp_t*       target_down = down_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_midpoint  = edge_midpoint_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_outside    = edge_outside_clv_ptr<fp_t>(D, target_id, site_off);
    const fp_t* target_up   = up_clv_ptr<const fp_t>(D, target_id, site_off);
    if (!parent_down || !target_down || !sibling_up || !target_up) return;

    const fp_t* target_mat  = D.d_pmat + (size_t)target_id  * (size_t)RATE_CATS * 16;
    const fp_t* target_mat_half = D.d_pmat_mid
        ? (D.d_pmat_mid + (size_t)target_id * (size_t)RATE_CATS * 16)
        : target_mat;
    const fp_t* sibling_mat = D.d_pmat + (size_t)sibling_id * (size_t)RATE_CATS * 16;
    unsigned int* parent_scaler = down_scaler_ptr(D, op.parent_id, site);
    unsigned int* sibling_scaler = up_scaler_ptr(D, sibling_id, site);
    unsigned int* target_up_scaler = up_scaler_ptr(D, target_id, site);
    unsigned int* down_scaler = down_scaler_ptr(D, target_id, site);
    unsigned int* edge_midpoint_scaler = edge_midpoint_scaler_ptr(D, target_id, site);
    unsigned int* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site);

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        write_downward_inherited_scalers_states4(
            D,
            parent_scaler,
            sibling_scaler,
            target_up_scaler,
            down_scaler,
            edge_midpoint_scaler,
            edge_outside_scaler,
            (unsigned int)r,
            edge_outside != nullptr);

        const fp_t* Tmat = target_mat  + (size_t)r * 16;
        const fp_t* Thalf= target_mat_half + (size_t)r * 16;
        const fp_t* Smat = sibling_mat + (size_t)r * 16;
        const fp4_t Ppar = reinterpret_cast<const fp4_t*>(parent_down)[r];
        const fp4_t Psib = reinterpret_cast<const fp4_t*>(sibling_up)[r];
        const fp4_t Pup  = reinterpret_cast<const fp4_t*>(target_up)[r];
        fp_t*       Pout = target_down + (size_t)r * 4;

        const fp_t sib0 = fp_dot4(make_fp4(Smat[0], Smat[1], Smat[2], Smat[3]), Psib);
        const fp_t sib1 = fp_dot4(make_fp4(Smat[4], Smat[5], Smat[6], Smat[7]), Psib);
        const fp_t sib2 = fp_dot4(make_fp4(Smat[8], Smat[9], Smat[10], Smat[11]), Psib);
        const fp_t sib3 = fp_dot4(make_fp4(Smat[12], Smat[13], Smat[14], Smat[15]), Psib);

        const fp_t p0 = Ppar.x * sib0;
        const fp_t p1 = Ppar.y * sib1;
        const fp_t p2 = Ppar.z * sib2;
        const fp_t p3 = Ppar.w * sib3;


        if (D.downward_pmat_indexing == DownwardPmatIndexing::Rows) {
            Pout[0] = Tmat[0] * p0 + Tmat[1] * p1 + Tmat[2] * p2 + Tmat[3] * p3;
            Pout[1] = Tmat[4] * p0 + Tmat[5] * p1 + Tmat[6] * p2 + Tmat[7] * p3;
            Pout[2] = Tmat[8] * p0 + Tmat[9] * p1 + Tmat[10] * p2 + Tmat[11] * p3;
            Pout[3] = Tmat[12] * p0 + Tmat[13] * p1 + Tmat[14] * p2 + Tmat[15] * p3;
        } else {
            Pout[0] = Tmat[0] * p0 + Tmat[4] * p1 + Tmat[8] * p2 + Tmat[12] * p3;
            Pout[1] = Tmat[1] * p0 + Tmat[5] * p1 + Tmat[9] * p2 + Tmat[13] * p3;
            Pout[2] = Tmat[2] * p0 + Tmat[6] * p1 + Tmat[10] * p2 + Tmat[14] * p3;
            Pout[3] = Tmat[3] * p0 + Tmat[7] * p1 + Tmat[11] * p2 + Tmat[15] * p3;
        }
        if (edge_outside) {
            fp_t* Pbase = edge_outside + (size_t)r * 4;
            Pbase[0] = p0;
            Pbase[1] = p1;
            Pbase[2] = p2;
            Pbase[3] = p3;
        }

        if (edge_midpoint) {
            fp_t* Pmid = edge_midpoint + (size_t)r * 4;
            build_midpoint_states4(Thalf, p0, p1, p2, p3, Pup, Pmid);
        }

        scale_clv_states4_if_needed(D, down_scaler, (unsigned int)r, Pout);
        if (edge_outside) {
            scale_clv_states4_if_needed(
                D,
                edge_outside_scaler,
                (unsigned int)r,
                edge_outside + (size_t)r * 4);
        }
        if (edge_midpoint) {
            scale_clv_states4_if_needed(
                D,
                edge_midpoint_scaler,
                (unsigned int)r,
                edge_midpoint + (size_t)r * 4);
        }
    }
}

template<int RATE_CATS>
__device__ __forceinline__ void compute_downward_inner_tip_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    if (op.clv_pool != static_cast<uint8_t>(CLV_POOL_DOWN)) return;
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_id       = target_is_left ? op.left_id  : op.right_id;
    const int sibling_tip_idx = target_is_left ? op.right_tip_index : op.left_tip_index;
    if (target_id < 0 || sibling_tip_idx < 0) return;

    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * (size_t)RATE_CATS * 4;

    const fp_t* parent_down = down_clv_ptr<const fp_t>(D, op.parent_id, site_off);
    fp_t*       target_down = down_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_midpoint  = edge_midpoint_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_outside    = edge_outside_clv_ptr<fp_t>(D, target_id, site_off);
    const fp_t* target_up   = up_clv_ptr<const fp_t>(D, target_id, site_off);
    if (!parent_down || !target_down || !target_up) return;

    const fp_t* target_mat  = D.d_pmat + (size_t)target_id * (size_t)RATE_CATS * 16;
    const fp_t* target_mat_half = D.d_pmat_mid
        ? (D.d_pmat_mid + (size_t)target_id * (size_t)RATE_CATS * 16)
        : target_mat;
    const fp_t* sibling_mat = D.d_pmat + (size_t)(target_is_left ? op.right_id : op.left_id) * (size_t)RATE_CATS * 16;
    unsigned int* parent_scaler = down_scaler_ptr(D, op.parent_id, site);
    unsigned int* sibling_scaler = up_scaler_ptr(D, target_is_left ? op.right_id : op.left_id, site);
    unsigned int* target_up_scaler = up_scaler_ptr(D, target_id, site);
    unsigned int* down_scaler = down_scaler_ptr(D, target_id, site);
    unsigned int* edge_midpoint_scaler = edge_midpoint_scaler_ptr(D, target_id, site);
    unsigned int* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site);

    const unsigned char* tipchars = D.d_tipchars + (size_t)sibling_tip_idx * D.sites;

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        write_downward_inherited_scalers_states4(
            D,
            parent_scaler,
            sibling_scaler,
            target_up_scaler,
            down_scaler,
            edge_midpoint_scaler,
            edge_outside_scaler,
            (unsigned int)r,
            edge_outside != nullptr);

        const unsigned int mask = D.d_tipmap[tipchars[site]];
        const fp_t* Tmat = target_mat  + (size_t)r * 16;
        const fp_t* Thalf= target_mat_half + (size_t)r * 16;
        const fp_t* Smat = sibling_mat + (size_t)r * 16;
        const fp4_t Ppar = reinterpret_cast<const fp4_t*>(parent_down)[r];
        const fp4_t Pup  = reinterpret_cast<const fp4_t*>(target_up)[r];
        fp_t*       Pout = target_down + (size_t)r * 4;

        const fp_t sib0 =
            ((mask & 1u) ? Smat[0] : fp_t(0)) +
            ((mask & 2u) ? Smat[1] : fp_t(0)) +
            ((mask & 4u) ? Smat[2] : fp_t(0)) +
            ((mask & 8u) ? Smat[3] : fp_t(0));
        const fp_t sib1 =
            ((mask & 1u) ? Smat[4] : fp_t(0)) +
            ((mask & 2u) ? Smat[5] : fp_t(0)) +
            ((mask & 4u) ? Smat[6] : fp_t(0)) +
            ((mask & 8u) ? Smat[7] : fp_t(0));
        const fp_t sib2 =
            ((mask & 1u) ? Smat[8] : fp_t(0)) +
            ((mask & 2u) ? Smat[9] : fp_t(0)) +
            ((mask & 4u) ? Smat[10] : fp_t(0)) +
            ((mask & 8u) ? Smat[11] : fp_t(0));
        const fp_t sib3 =
            ((mask & 1u) ? Smat[12] : fp_t(0)) +
            ((mask & 2u) ? Smat[13] : fp_t(0)) +
            ((mask & 4u) ? Smat[14] : fp_t(0)) +
            ((mask & 8u) ? Smat[15] : fp_t(0));

        const fp_t p0 = Ppar.x * sib0;
        const fp_t p1 = Ppar.y * sib1;
        const fp_t p2 = Ppar.z * sib2;
        const fp_t p3 = Ppar.w * sib3;

        if (D.downward_pmat_indexing == DownwardPmatIndexing::Rows) {
            Pout[0] = Tmat[0] * p0 + Tmat[1] * p1 + Tmat[2] * p2 + Tmat[3] * p3;
            Pout[1] = Tmat[4] * p0 + Tmat[5] * p1 + Tmat[6] * p2 + Tmat[7] * p3;
            Pout[2] = Tmat[8] * p0 + Tmat[9] * p1 + Tmat[10] * p2 + Tmat[11] * p3;
            Pout[3] = Tmat[12] * p0 + Tmat[13] * p1 + Tmat[14] * p2 + Tmat[15] * p3;
        } else {
            Pout[0] = Tmat[0] * p0 + Tmat[4] * p1 + Tmat[8] * p2 + Tmat[12] * p3;
            Pout[1] = Tmat[1] * p0 + Tmat[5] * p1 + Tmat[9] * p2 + Tmat[13] * p3;
            Pout[2] = Tmat[2] * p0 + Tmat[6] * p1 + Tmat[10] * p2 + Tmat[14] * p3;
            Pout[3] = Tmat[3] * p0 + Tmat[7] * p1 + Tmat[11] * p2 + Tmat[15] * p3;
        }
        if (edge_outside) {
            fp_t* Pbase = edge_outside + (size_t)r * 4;
            Pbase[0] = p0;
            Pbase[1] = p1;
            Pbase[2] = p2;
            Pbase[3] = p3;
        }

        if (edge_midpoint) {
            fp_t* Pmid = edge_midpoint + (size_t)r * 4;
            build_midpoint_states4(Thalf, p0, p1, p2, p3, Pup, Pmid);
        }

        scale_clv_states4_if_needed(D, down_scaler, (unsigned int)r, Pout);
        if (edge_outside) {
            scale_clv_states4_if_needed(
                D,
                edge_outside_scaler,
                (unsigned int)r,
                edge_outside + (size_t)r * 4);
        }
        if (edge_midpoint) {
            scale_clv_states4_if_needed(
                D,
                edge_midpoint_scaler,
                (unsigned int)r,
                edge_midpoint + (size_t)r * 4);
        }
    }
}

template<int RATE_CATS>
__device__ __forceinline__ void compute_downward_tip_inner_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    if (op.clv_pool != static_cast<uint8_t>(CLV_POOL_DOWN)) return;
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_tip_idx  = target_is_left ? op.left_tip_index : op.right_tip_index;
    const int target_id       = target_is_left ? op.left_id : op.right_id;
    const int sibling_id      = target_is_left ? op.right_id : op.left_id;
    if (target_tip_idx < 0 || sibling_id < 0) return;

    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * (size_t)RATE_CATS * 4;

    const fp_t* parent_down = down_clv_ptr<const fp_t>(D, op.parent_id, site_off);
    const fp_t* sibling_up  = up_clv_ptr<const fp_t>(D, sibling_id, site_off);
    fp_t*       target_down = down_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_midpoint  = edge_midpoint_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_outside    = edge_outside_clv_ptr<fp_t>(D, target_id, site_off);
    const fp_t* target_up   = up_clv_ptr<const fp_t>(D, target_id, site_off);
    if (!parent_down || !target_down || !sibling_up || !target_up) return;

    const fp_t* target_mat  = D.d_pmat + (size_t)target_id * (size_t)RATE_CATS * 16;
    const fp_t* target_mat_half = D.d_pmat_mid
        ? (D.d_pmat_mid + (size_t)target_id * (size_t)RATE_CATS * 16)
        : target_mat;
    const fp_t* sibling_mat = D.d_pmat + (size_t)sibling_id * (size_t)RATE_CATS * 16;
    unsigned int* parent_scaler = down_scaler_ptr(D, op.parent_id, site);
    unsigned int* sibling_scaler = up_scaler_ptr(D, sibling_id, site);
    unsigned int* target_up_scaler = up_scaler_ptr(D, target_id, site);
    unsigned int* down_scaler = down_scaler_ptr(D, target_id, site);
    unsigned int* edge_midpoint_scaler = edge_midpoint_scaler_ptr(D, target_id, site);
    unsigned int* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site);

    const unsigned char* tipchars = D.d_tipchars + (size_t)target_tip_idx * D.sites;
    const unsigned int tmask = D.d_tipmap[tipchars[site]];

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        write_downward_inherited_scalers_states4(
            D,
            parent_scaler,
            sibling_scaler,
            target_up_scaler,
            down_scaler,
            edge_midpoint_scaler,
            edge_outside_scaler,
            (unsigned int)r,
            edge_outside != nullptr);

        const fp_t* Tmat  = target_mat  + (size_t)r * 16;
        const fp_t* Thalf = target_mat_half + (size_t)r * 16;
        const fp_t* Smat  = sibling_mat + (size_t)r * 16;
        const fp4_t Ppar = reinterpret_cast<const fp4_t*>(parent_down)[r];
        const fp4_t Psib = reinterpret_cast<const fp4_t*>(sibling_up)[r];
        const fp4_t Pup  = reinterpret_cast<const fp4_t*>(target_up)[r];
        fp_t*       Pout = target_down + (size_t)r * 4;

        const fp_t sib0 = fp_dot4(make_fp4(Smat[0], Smat[1], Smat[2], Smat[3]), Psib);
        const fp_t sib1 = fp_dot4(make_fp4(Smat[4], Smat[5], Smat[6], Smat[7]), Psib);
        const fp_t sib2 = fp_dot4(make_fp4(Smat[8], Smat[9], Smat[10], Smat[11]), Psib);
        const fp_t sib3 = fp_dot4(make_fp4(Smat[12], Smat[13], Smat[14], Smat[15]), Psib);

        const fp_t p0 = Ppar.x * sib0;
        const fp_t p1 = Ppar.y * sib1;
        const fp_t p2 = Ppar.z * sib2;
        const fp_t p3 = Ppar.w * sib3;


        if (D.downward_pmat_indexing == DownwardPmatIndexing::Rows) {
            Pout[0] = Tmat[0] * p0 + Tmat[1] * p1 + Tmat[2] * p2 + Tmat[3] * p3;
            Pout[1] = Tmat[4] * p0 + Tmat[5] * p1 + Tmat[6] * p2 + Tmat[7] * p3;
            Pout[2] = Tmat[8] * p0 + Tmat[9] * p1 + Tmat[10] * p2 + Tmat[11] * p3;
            Pout[3] = Tmat[12] * p0 + Tmat[13] * p1 + Tmat[14] * p2 + Tmat[15] * p3;
        } else {
            Pout[0] = Tmat[0] * p0 + Tmat[4] * p1 + Tmat[8] * p2 + Tmat[12] * p3;
            Pout[1] = Tmat[1] * p0 + Tmat[5] * p1 + Tmat[9] * p2 + Tmat[13] * p3;
            Pout[2] = Tmat[2] * p0 + Tmat[6] * p1 + Tmat[10] * p2 + Tmat[14] * p3;
            Pout[3] = Tmat[3] * p0 + Tmat[7] * p1 + Tmat[11] * p2 + Tmat[15] * p3;
        }
        if (!(tmask & 1u)) Pout[0] = 0.0;
        if (!(tmask & 2u)) Pout[1] = 0.0;
        if (!(tmask & 4u)) Pout[2] = 0.0;
        if (!(tmask & 8u)) Pout[3] = 0.0;

        if (edge_outside) {
            fp_t* Pbase = edge_outside + (size_t)r * 4;
            Pbase[0] = p0;
            Pbase[1] = p1;
            Pbase[2] = p2;
            Pbase[3] = p3;
        }

        if (edge_midpoint) {
            fp_t* Pmid = edge_midpoint + (size_t)r * 4;
            build_midpoint_states4(Thalf, p0, p1, p2, p3, Pup, Pmid);
            if (!(tmask & 1u)) Pmid[0] = 0.0;
            if (!(tmask & 2u)) Pmid[1] = 0.0;
            if (!(tmask & 4u)) Pmid[2] = 0.0;
            if (!(tmask & 8u)) Pmid[3] = 0.0;
        }

        scale_clv_states4_if_needed(D, down_scaler, (unsigned int)r, Pout);
        if (edge_outside) {
            scale_clv_states4_if_needed(
                D,
                edge_outside_scaler,
                (unsigned int)r,
                edge_outside + (size_t)r * 4);
        }
        if (edge_midpoint) {
            scale_clv_states4_if_needed(
                D,
                edge_midpoint_scaler,
                (unsigned int)r,
                edge_midpoint + (size_t)r * 4);
        }
    }
}

// target tip, sibling tip (states=4, rate-specific)
template<int RATE_CATS>
__device__ __forceinline__ void compute_downward_tip_tip_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    if (op.clv_pool != static_cast<uint8_t>(CLV_POOL_DOWN)) return;
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_tip_idx  = target_is_left ? op.left_tip_index : op.right_tip_index;
    const int sibling_tip_idx = target_is_left ? op.right_tip_index : op.left_tip_index;
    const int target_id       = target_is_left ? op.left_id : op.right_id;
    if (target_tip_idx < 0 || sibling_tip_idx < 0 || target_id < 0) return;

    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * (size_t)RATE_CATS * 4;

    const fp_t* parent_down = down_clv_ptr<const fp_t>(D, op.parent_id, site_off);
    fp_t*       target_down = down_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_midpoint  = edge_midpoint_clv_ptr<fp_t>(D, target_id, site_off);
    fp_t*       edge_outside    = edge_outside_clv_ptr<fp_t>(D, target_id, site_off);
    const fp_t* target_up   = up_clv_ptr<const fp_t>(D, target_id, site_off);
    if (!parent_down || !target_down || !target_up) return;

    const fp_t* target_mat  = D.d_pmat + (size_t)target_id * (size_t)RATE_CATS * 16;
    const fp_t* target_mat_half = D.d_pmat_mid
        ? (D.d_pmat_mid + (size_t)target_id * (size_t)RATE_CATS * 16)
        : target_mat;
    const fp_t* sibling_mat = D.d_pmat + (size_t)(target_is_left ? op.right_id : op.left_id) * (size_t)RATE_CATS * 16;
    unsigned int* parent_scaler = down_scaler_ptr(D, op.parent_id, site);
    unsigned int* sibling_scaler = up_scaler_ptr(D, target_is_left ? op.right_id : op.left_id, site);
    unsigned int* target_up_scaler = up_scaler_ptr(D, target_id, site);
    unsigned int* down_scaler = down_scaler_ptr(D, target_id, site);
    unsigned int* edge_midpoint_scaler = edge_midpoint_scaler_ptr(D, target_id, site);
    unsigned int* edge_outside_scaler = edge_outside_scaler_ptr(D, target_id, site);

    const unsigned char* tipchars = D.d_tipchars + (size_t)sibling_tip_idx * D.sites;

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        write_downward_inherited_scalers_states4(
            D,
            parent_scaler,
            sibling_scaler,
            target_up_scaler,
            down_scaler,
            edge_midpoint_scaler,
            edge_outside_scaler,
            (unsigned int)r,
            edge_outside != nullptr);

        const unsigned int mask = D.d_tipmap[tipchars[site]];
        const fp_t* Tmat  = target_mat  + (size_t)r * 16;
        const fp_t* Thalf = target_mat_half + (size_t)r * 16;
        const fp_t* Smat  = sibling_mat + (size_t)r * 16;
        const fp4_t Ppar = reinterpret_cast<const fp4_t*>(parent_down)[r];
        const fp4_t Pup  = reinterpret_cast<const fp4_t*>(target_up)[r];
        fp_t*       Pout = target_down + (size_t)r * 4;

        const fp_t sib0 =
            ((mask & 1u) ? Smat[0] : fp_t(0)) +
            ((mask & 2u) ? Smat[1] : fp_t(0)) +
            ((mask & 4u) ? Smat[2] : fp_t(0)) +
            ((mask & 8u) ? Smat[3] : fp_t(0));
        const fp_t sib1 =
            ((mask & 1u) ? Smat[4] : fp_t(0)) +
            ((mask & 2u) ? Smat[5] : fp_t(0)) +
            ((mask & 4u) ? Smat[6] : fp_t(0)) +
            ((mask & 8u) ? Smat[7] : fp_t(0));
        const fp_t sib2 =
            ((mask & 1u) ? Smat[8] : fp_t(0)) +
            ((mask & 2u) ? Smat[9] : fp_t(0)) +
            ((mask & 4u) ? Smat[10] : fp_t(0)) +
            ((mask & 8u) ? Smat[11] : fp_t(0));
        const fp_t sib3 =
            ((mask & 1u) ? Smat[12] : fp_t(0)) +
            ((mask & 2u) ? Smat[13] : fp_t(0)) +
            ((mask & 4u) ? Smat[14] : fp_t(0)) +
            ((mask & 8u) ? Smat[15] : fp_t(0));

        const fp_t p0 = Ppar.x * sib0;
        const fp_t p1 = Ppar.y * sib1;
        const fp_t p2 = Ppar.z * sib2;
        const fp_t p3 = Ppar.w * sib3;

        if (D.downward_pmat_indexing == DownwardPmatIndexing::Rows) {
            Pout[0] = Tmat[0] * p0 + Tmat[1] * p1 + Tmat[2] * p2 + Tmat[3] * p3;
            Pout[1] = Tmat[4] * p0 + Tmat[5] * p1 + Tmat[6] * p2 + Tmat[7] * p3;
            Pout[2] = Tmat[8] * p0 + Tmat[9] * p1 + Tmat[10] * p2 + Tmat[11] * p3;
            Pout[3] = Tmat[12] * p0 + Tmat[13] * p1 + Tmat[14] * p2 + Tmat[15] * p3;
        } else {
            Pout[0] = Tmat[0] * p0 + Tmat[4] * p1 + Tmat[8] * p2 + Tmat[12] * p3;
            Pout[1] = Tmat[1] * p0 + Tmat[5] * p1 + Tmat[9] * p2 + Tmat[13] * p3;
            Pout[2] = Tmat[2] * p0 + Tmat[6] * p1 + Tmat[10] * p2 + Tmat[14] * p3;
            Pout[3] = Tmat[3] * p0 + Tmat[7] * p1 + Tmat[11] * p2 + Tmat[15] * p3;
        }

        if (edge_outside) {
            fp_t* Pbase = edge_outside + (size_t)r * 4;
            Pbase[0] = p0;
            Pbase[1] = p1;
            Pbase[2] = p2;
            Pbase[3] = p3;
        }

        if (edge_midpoint) {
            fp_t* Pmid = edge_midpoint + (size_t)r * 4;
            build_midpoint_states4(Thalf, p0, p1, p2, p3, Pup, Pmid);
        }

        scale_clv_states4_if_needed(D, down_scaler, (unsigned int)r, Pout);
        if (edge_outside) {
            scale_clv_states4_if_needed(
                D,
                edge_outside_scaler,
                (unsigned int)r,
                edge_outside + (size_t)r * 4);
        }
        if (edge_midpoint) {
            scale_clv_states4_if_needed(
                D,
                edge_midpoint_scaler,
                (unsigned int)r,
                edge_midpoint + (size_t)r * 4);
        }
    }
}

// Upward dispatch preserves the [site][rate][state] CLV layout and one scaler
// shift per site/rate. DNA4 uses compile-time rate counts when available; other
// state/rate combinations retain the same contract through generic helpers.
template<int RATE_CATS>
__device__ __forceinline__ void compute_tip_tip_site_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const size_t span = (size_t)4 * RATE_CATS;

    const unsigned char* left_tip  = D.d_tipchars + (size_t)op.left_tip_index  * D.sites;
    const unsigned char* right_tip = D.d_tipchars + (size_t)op.right_tip_index * D.sites;

    const unsigned int j = (unsigned int)left_tip[site];
    const unsigned int k = (unsigned int)right_tip[site];

    const unsigned int jmask_base = D.d_tipmap[j];
    const unsigned int kmask_base = D.d_tipmap[k];

    const fp_t* __restrict__ jmat_base =
        D.d_pmat + (size_t)op.left_id  * RATE_CATS * 4 * 4;
    const fp_t* __restrict__ kmat_base =
        D.d_pmat + (size_t)op.right_id * RATE_CATS * 4 * 4;

    fp_t* parent_clv = clv_write_ptr_for_node<fp_t>(D, op, op.parent_id);
    if (!parent_clv) return;
    fp_t* __restrict__ dst = parent_clv + (size_t)site * span;

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, RATE_CATS);

#pragma unroll
    for (int r = 0; r < RATE_CATS; ++r) {
        write_clv_scaler_shift(D, site_scaler_ptr, r, 0u);
        const fp_t* __restrict__ jmat = jmat_base + (size_t)r * 4 * 4;
        const fp_t* __restrict__ kmat = kmat_base + (size_t)r * 4 * 4;
        fp_t* __restrict__ Pout = dst + (size_t)r * 4;
        const fp_t max_val =
            compute_tip_tip_states4_rate(jmat, kmat, jmask_base, kmask_base, Pout);
        scale_clv_states4_if_needed(D, site_scaler_ptr, (unsigned int)r, Pout, max_val);
    }
}

__device__ __forceinline__ void compute_tip_tip_site_4_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const size_t span     = (size_t)4 * (size_t)D.rate_cats;

    const unsigned char* left_tip  = D.d_tipchars + (size_t)op.left_tip_index  * D.sites;
    const unsigned char* right_tip = D.d_tipchars + (size_t)op.right_tip_index * D.sites;

    const unsigned int j = (unsigned int)left_tip[site];
    const unsigned int k = (unsigned int)right_tip[site];

    const unsigned int jmask_base = D.d_tipmap[j];
    const unsigned int kmask_base = D.d_tipmap[k];

    const fp_t* __restrict__ jmat_base =
        D.d_pmat + (size_t)op.left_id  * D.rate_cats * 4 * 4;
    const fp_t* __restrict__ kmat_base =
        D.d_pmat + (size_t)op.right_id * D.rate_cats * 4 * 4;

    fp_t* parent_clv = clv_write_ptr_for_node<fp_t>(D, op, op.parent_id);
    if (!parent_clv) return;
    fp_t* __restrict__ dst = parent_clv + (size_t)site * span;

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, (unsigned int)D.rate_cats);

    for (int r = 0; r < D.rate_cats; ++r) {
        write_clv_scaler_shift(D, site_scaler_ptr, r, 0u);
        const fp_t* __restrict__ jmat = jmat_base + (size_t)r * 4 * 4;
        const fp_t* __restrict__ kmat = kmat_base + (size_t)r * 4 * 4;
        fp_t* __restrict__ Pout = dst + (size_t)r * 4;
        const fp_t max_val =
            compute_tip_tip_states4_rate(jmat, kmat, jmask_base, kmask_base, Pout);
        scale_clv_states4_if_needed(D, site_scaler_ptr, (unsigned int)r, Pout, max_val);
    }
}

__device__ __forceinline__ void compute_tip_tip_site_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const unsigned int states = (unsigned int)D.states;
    const unsigned int rate_cats = (unsigned int)D.rate_cats;
    const size_t span     = (size_t)states * rate_cats;

    const unsigned char* left_tip  = D.d_tipchars + (size_t)op.left_tip_index  * D.sites;
    const unsigned char* right_tip = D.d_tipchars + (size_t)op.right_tip_index * D.sites;

    const unsigned int lmask = D.d_tipmap[left_tip[site]];
    const unsigned int rmask = D.d_tipmap[right_tip[site]];

    const fp_t* Lbase = D.d_pmat + (size_t)op.left_id  * rate_cats * states * states;
    const fp_t* Rbase = D.d_pmat + (size_t)op.right_id * rate_cats * states * states;

    fp_t* parent_clv = clv_write_ptr_for_node<fp_t>(D, op, op.parent_id);
    if (!parent_clv) return;
    fp_t* Pout = parent_clv + (size_t)site * span;

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, rate_cats);

    for (unsigned int r = 0; r < rate_cats; ++r) {
        const fp_t* Lmat = Lbase + (size_t)r * states * states;
        const fp_t* Rmat = Rbase + (size_t)r * states * states;
        fp_t* out_r = Pout + (size_t)r * states;

        fp_t maxv = fp_t(0);
        for (unsigned int j = 0; j < states; ++j) {
            fp_t left_term = fp_t(0);
            fp_t right_term = fp_t(0);
            for (unsigned int k = 0; k < states; ++k) {
                if (lmask & (1u << k)) left_term  += Lmat[j * states + k];
                if (rmask & (1u << k)) right_term += Rmat[j * states + k];
            }
            fp_t v = left_term * right_term;
            out_r[j] = v;
            if (v > maxv) maxv = v;
        }

        if (site_scaler_ptr) {
            unsigned int shift = clv_scale_shift(maxv);
            if (shift) {
                add_clv_scaler_shift(D, site_scaler_ptr, r, shift);
                for (unsigned int s = 0; s < states; ++s) {
                    fp_scale_pow2(out_r[s], shift);
                }
            }
        }
    }
}

template<int RATE_CATS>
__device__ __forceinline__ void compute_tip_inner_site_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const size_t span     = (size_t)4 * RATE_CATS;
    const size_t per_node = per_node_span(D);

    const bool tip_on_left = op.left_tip_index >= 0;
    const int  tip_index   = tip_on_left ? op.left_tip_index  : op.right_tip_index;
    const int  inner_id    = tip_on_left ? op.right_id : op.left_id;
    const int  tip_node_id = tip_on_left ? op.left_id  : op.right_id;

    const unsigned char* d_left_tip = D.d_tipchars + (size_t)tip_index * D.sites;
    const fp_t* d_right_clv = clv_read_ptr_for_node<const fp_t>(D, op, inner_id);
    fp_t* parent_clv = clv_write_ptr_for_node<fp_t>(D, op, op.parent_id);
    if (!d_right_clv || !parent_clv) return;

    const fp_t* d_Lmat = D.d_pmat + (size_t)tip_node_id * RATE_CATS * 4 * 4;
    const fp_t* d_Rmat = D.d_pmat + (size_t)inner_id * RATE_CATS * 4 * 4;

    const size_t site_off = (size_t)site * span;
    const unsigned int tmask = D.d_tipmap[d_left_tip[site]];

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, RATE_CATS);
    unsigned int* inner_scaler =
        scaler_ptr_for_pool(D, op.clv_pool, inner_id, site);

    for (int r = 0; r < RATE_CATS; ++r) {
        write_clv_scaler_shift(D, site_scaler_ptr, r, read_clv_scaler_shift(D, inner_scaler, r));
        const fp_t* Lmat = d_Lmat + (size_t)r * 4 * 4;
        const fp_t* Rmat = d_Rmat + (size_t)r * 4 * 4;
        const fp_t* Rclv = d_right_clv + site_off + (size_t)r * 4;
        fp_t* Pout = parent_clv + site_off + (size_t)r * 4;
        const fp_t max_val = compute_tip_inner_states4_rate(Lmat, Rmat, Rclv, tmask, Pout);
        scale_clv_states4_if_needed(D, site_scaler_ptr, (unsigned int)r, Pout, max_val);
    }
}

__device__ __forceinline__ void compute_tip_inner_site_4_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const size_t span     = (size_t)4 * (size_t)D.rate_cats;
    const size_t per_node = per_node_span(D);

    const bool tip_on_left = op.left_tip_index >= 0;
    const int tip_index    = tip_on_left ? op.left_tip_index : op.right_tip_index;
    const int inner_id     = tip_on_left ? op.right_id : op.left_id;
    const int tip_node_id  = tip_on_left ? op.left_id : op.right_id;

    const unsigned char* tip_chars = D.d_tipchars + (size_t)tip_index * D.sites;
    const fp_t* inner_clv = clv_read_ptr_for_node<const fp_t>(D, op, inner_id);
    fp_t* parent_clv = clv_write_ptr_for_node<fp_t>(D, op, op.parent_id);
    if (!inner_clv || !parent_clv) return;

    const fp_t* tip_mat_base = D.d_pmat + (size_t)tip_node_id * D.rate_cats * 4 * 4;
    const fp_t* inner_mat_base = D.d_pmat + (size_t)inner_id * D.rate_cats * 4 * 4;
    const size_t site_off = (size_t)site * span;
    const unsigned int tip_mask = D.d_tipmap[tip_chars[site]];

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, (unsigned int)D.rate_cats);
    unsigned int* inner_scaler =
        scaler_ptr_for_pool(D, op.clv_pool, inner_id, site);

    for (int r = 0; r < D.rate_cats; ++r) {
        write_clv_scaler_shift(D, site_scaler_ptr, r, read_clv_scaler_shift(D, inner_scaler, r));
        const fp_t* tip_mat = tip_mat_base + (size_t)r * 4 * 4;
        const fp_t* inner_mat = inner_mat_base + (size_t)r * 4 * 4;
        const fp_t* right_clv = inner_clv + site_off + (size_t)r * 4;
        fp_t* out = parent_clv + site_off + (size_t)r * 4;
        const fp_t max_val =
            compute_tip_inner_states4_rate(tip_mat, inner_mat, right_clv, tip_mask, out);
        scale_clv_states4_if_needed(D, site_scaler_ptr, (unsigned int)r, out, max_val);
    }
}

__device__ __forceinline__ void compute_tip_inner_site_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const unsigned int states = (unsigned int)D.states;
    const unsigned int rate_cats = (unsigned int)D.rate_cats;
    const size_t span     = (size_t)states * rate_cats;
    const size_t per_node = per_node_span(D);

    const bool tip_on_left = op.left_tip_index >= 0;
    const int  tip_index   = tip_on_left ? op.left_tip_index  : op.right_tip_index;
    const int  inner_id    = tip_on_left ? op.right_id : op.left_id;
    const int  tip_node_id = tip_on_left ? op.left_id  : op.right_id;

    const unsigned char* d_left_tip = D.d_tipchars + (size_t)tip_index * D.sites;
    const fp_t* d_right_clv = clv_read_ptr_for_node<const fp_t>(D, op, inner_id);
    fp_t* parent_clv = clv_write_ptr_for_node<fp_t>(D, op, op.parent_id);
    if (!d_right_clv || !parent_clv) return;

    const fp_t* d_Lmat = D.d_pmat + (size_t)tip_node_id * D.rate_cats * states * states;
    const fp_t* d_Rmat = D.d_pmat + (size_t)inner_id * D.rate_cats * states * states;

    const size_t site_off = (size_t)site * span;
    const unsigned int tmask = D.d_tipmap[d_left_tip[site]];

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, rate_cats);
    unsigned int* inner_scaler =
        scaler_ptr_for_pool(D, op.clv_pool, inner_id, site);

    for (unsigned int r = 0; r < rate_cats; ++r) {
        write_clv_scaler_shift(D, site_scaler_ptr, r, read_clv_scaler_shift(D, inner_scaler, r));
        fp_t col_scale_max_val = fp_t(0);
        const fp_t* Lmat = d_Lmat + (size_t)r * states * states;
        const fp_t* Rmat = d_Rmat + (size_t)r * states * states;
        const fp_t* Rclv = d_right_clv + site_off + (size_t)r * states;
        fp_t* Pout = parent_clv + site_off + (size_t)r * states;

        const fp_t* Lrow = Lmat;
        const fp_t* Rrow = Rmat;
        for (unsigned int i = 0; i < states; ++i) {
            fp_t lefterm = fp_t(0), righterm = fp_t(0);
            unsigned int lstate = tmask;
            for (unsigned int j = 0; j < states; ++j) {
                if (lstate & 1u) lefterm += Lrow[j];
                righterm += Rrow[j] * Rclv[j];
                lstate >>= 1;
            }
            Pout[i] = lefterm * righterm;
            if (Pout[i] > col_scale_max_val) col_scale_max_val = Pout[i];
            Lrow += states;
            Rrow += states;
        }
        if (site_scaler_ptr) {
            unsigned int shift = clv_scale_shift(col_scale_max_val);
            if (shift) {
                add_clv_scaler_shift(D, site_scaler_ptr, r, shift);
                for (unsigned int i = 0; i < states; ++i) {
                    fp_scale_pow2(Pout[i], shift);
                }
            }
        }
    }
}

template<int RATE_CATS>
__device__ __forceinline__ void compute_inner_inner_site_ratecat(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const size_t span     = (size_t)RATE_CATS * 4;
    const size_t site_off = (size_t)site * span;

    const fp_t* d_left_clv  = clv_read_ptr_for_node<const fp_t>(D, op, op.left_id);
    const fp_t* d_right_clv = clv_read_ptr_for_node<const fp_t>(D, op, op.right_id);

    fp_t* parent_clv = clv_write_ptr_for_node<fp_t>(D, op, op.parent_id);
    if (!d_left_clv || !d_right_clv || !parent_clv) return;
    const fp_t* d_left_mat  = D.d_pmat + (size_t)op.left_id  * RATE_CATS * 4 * 4;
    const fp_t* d_right_mat = D.d_pmat + (size_t)op.right_id * RATE_CATS * 4 * 4;

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, RATE_CATS);
    unsigned int* left_scaler =
        scaler_ptr_for_pool(D, op.clv_pool, op.left_id, site);
    unsigned int* right_scaler =
        scaler_ptr_for_pool(D, op.clv_pool, op.right_id, site);

    for (int r = 0; r < RATE_CATS; ++r) {
        write_clv_scaler_shift(
            D,
            site_scaler_ptr,
            r,
            read_clv_scaler_shift(D, left_scaler, r) +
            read_clv_scaler_shift(D, right_scaler, r));
        const fp_t* Lclv = d_left_clv  + site_off + (size_t)r * 4;
        const fp_t* Rclv = d_right_clv + site_off + (size_t)r * 4;
        const fp_t* Lmat = d_left_mat  + (size_t)r * 4 * 4;
        const fp_t* Rmat = d_right_mat + (size_t)r * 4 * 4;
        fp_t* Pout = parent_clv + site_off + (size_t)r * 4;
        const fp_t max_val = compute_inner_inner_states4_rate(Lmat, Rmat, Lclv, Rclv, Pout);
        scale_clv_states4_if_needed(D, site_scaler_ptr, (unsigned int)r, Pout, max_val);
    }
}

__device__ __forceinline__ void compute_inner_inner_site_4_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const size_t span = (size_t)4 * (size_t)D.rate_cats;
    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * span;

    const fp_t* left_clv = clv_read_ptr_for_node<const fp_t>(D, op, op.left_id);
    const fp_t* right_clv = clv_read_ptr_for_node<const fp_t>(D, op, op.right_id);
    fp_t* parent_clv = clv_write_ptr_for_node<fp_t>(D, op, op.parent_id);
    if (!left_clv || !right_clv || !parent_clv) return;

    const fp_t* left_mat_base = D.d_pmat + (size_t)op.left_id * (size_t)D.rate_cats * 4 * 4;
    const fp_t* right_mat_base = D.d_pmat + (size_t)op.right_id * (size_t)D.rate_cats * 4 * 4;
    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, (unsigned int)D.rate_cats);
    unsigned int* left_scaler =
        scaler_ptr_for_pool(D, op.clv_pool, op.left_id, site);
    unsigned int* right_scaler =
        scaler_ptr_for_pool(D, op.clv_pool, op.right_id, site);

    for (int r = 0; r < D.rate_cats; ++r) {
        write_clv_scaler_shift(
            D,
            site_scaler_ptr,
            r,
            read_clv_scaler_shift(D, left_scaler, r) +
            read_clv_scaler_shift(D, right_scaler, r));
        const fp_t* left_clv_r = left_clv + site_off + (size_t)r * 4;
        const fp_t* right_clv_r = right_clv + site_off + (size_t)r * 4;
        const fp_t* left_mat = left_mat_base + (size_t)r * 4 * 4;
        const fp_t* right_mat = right_mat_base + (size_t)r * 4 * 4;
        fp_t* out = parent_clv + site_off + (size_t)r * 4;
        const fp_t max_val =
            compute_inner_inner_states4_rate(left_mat, right_mat, left_clv_r, right_clv_r, out);
        scale_clv_states4_if_needed(D, site_scaler_ptr, (unsigned int)r, out, max_val);
    }
}

__device__ __forceinline__ void compute_inner_inner_site_generic(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    const unsigned int states = (unsigned int)D.states;
    const unsigned int rate_cats = (unsigned int)D.rate_cats;
    const size_t span = (size_t)states * (size_t)rate_cats;
    const size_t per_node = per_node_span(D);
    const size_t site_off = (size_t)site * span;

    const fp_t* d_left_clv  = clv_read_ptr_for_node<const fp_t>(D, op, op.left_id);
    const fp_t* d_right_clv = clv_read_ptr_for_node<const fp_t>(D, op, op.right_id);
    fp_t* parent_clv = clv_write_ptr_for_node<fp_t>(D, op, op.parent_id);
    if (!d_left_clv || !d_right_clv || !parent_clv) return;
    const fp_t* d_left_mat  = D.d_pmat + (size_t)op.left_id  * D.rate_cats * states * states;
    const fp_t* d_right_mat = D.d_pmat + (size_t)op.right_id * D.rate_cats * states * states;

    unsigned int* site_scaler_ptr =
        site_scaler_ptr_base(D, op, site, rate_cats);
    unsigned int* left_scaler =
        scaler_ptr_for_pool(D, op.clv_pool, op.left_id, site);
    unsigned int* right_scaler =
        scaler_ptr_for_pool(D, op.clv_pool, op.right_id, site);

    for (unsigned int r = 0; r < rate_cats; ++r) {
        write_clv_scaler_shift(
            D,
            site_scaler_ptr,
            r,
            read_clv_scaler_shift(D, left_scaler, r) +
            read_clv_scaler_shift(D, right_scaler, r));
        const fp_t* Lclv = d_left_clv  + site_off + (size_t)r * states;
        const fp_t* Rclv = d_right_clv + site_off + (size_t)r * states;

        const fp_t* Lmat = d_left_mat  + (size_t)r * states * states;
        const fp_t* Rmat = d_right_mat + (size_t)r * states * states;

        fp_t* Pout = parent_clv + site_off + (size_t)r * states;
        fp_t col_scale_max_val = fp_t(0);

        const fp_t* Lrow = Lmat;
        const fp_t* Rrow = Rmat;
        for (unsigned int j = 0; j < states; ++j) {
            fp_t lt = fp_t(0), rt = fp_t(0);
#pragma unroll
            for (unsigned int k = 0; k < states; ++k) {
                lt = fp_fma(Lrow[k], Lclv[k], lt);
                rt = fp_fma(Rrow[k], Rclv[k], rt);
            }
            Pout[j] = lt * rt;
            if (Pout[j] > col_scale_max_val) col_scale_max_val = Pout[j];
            Lrow += states;
            Rrow += states;
        }

        if (site_scaler_ptr) {
            unsigned int shift = clv_scale_shift(col_scale_max_val);
            if (shift) {
                add_clv_scaler_shift(D, site_scaler_ptr, r, shift);
                for (unsigned int j = 0; j < states; ++j) {
                    fp_scale_pow2(Pout[j], shift);
                }
            }
        }
    }
}

__device__ __forceinline__ void execute_upward_op_for_site(
    const DeviceTree& D, const NodeOpInfo& op, unsigned int site)
{
    switch (op.op_type) {
        case OP_TIP_TIP:
            if (D.states != 4) { compute_tip_tip_site_generic(D, op, site); break; }
            switch (D.rate_cats) {
                case 1: compute_tip_tip_site_ratecat<1>(D, op, site); break;
                case 4: compute_tip_tip_site_ratecat<4>(D, op, site); break;
                case 8: compute_tip_tip_site_ratecat<8>(D, op, site); break;
                default: compute_tip_tip_site_4_generic(D, op, site); break;
            }
            break;
        case OP_TIP_INNER:
            if (D.states != 4) { compute_tip_inner_site_generic(D, op, site); break; }
            switch (D.rate_cats) {
                case 1: compute_tip_inner_site_ratecat<1>(D, op, site); break;
                case 4: compute_tip_inner_site_ratecat<4>(D, op, site); break;
                case 8: compute_tip_inner_site_ratecat<8>(D, op, site); break;
                default: compute_tip_inner_site_4_generic(D, op, site); break;
            }
            break;
        case OP_INNER_INNER:
            if (D.states != 4) { compute_inner_inner_site_generic(D, op, site); break; }
            switch (D.rate_cats) {
                case 1: compute_inner_inner_site_ratecat<1>(D, op, site); break;
                case 4: compute_inner_inner_site_ratecat<4>(D, op, site); break;
                case 8: compute_inner_inner_site_ratecat<8>(D, op, site); break;
                default: compute_inner_inner_site_4_generic(D, op, site); break;
            }
            break;
        default: break;
    }
}

__global__ void UpdatePartialsUpwardKernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops
) {
    unsigned int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int step = blockDim.x * gridDim.x;

    for (unsigned int site = tid; site < D.sites; site += step) {
        for (int i = 0; i < num_ops; ++i) {
            execute_upward_op_for_site(D, ops[i], site);
        }
    }
}
__global__ void UpdatePartialsUpwardLevelKernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops)
{
    if (blockIdx.y >= static_cast<unsigned int>(num_ops)) return;
    const NodeOpInfo op = ops[blockIdx.y];
    const unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int step = blockDim.x * gridDim.x;
    for (unsigned int site = tid; site < D.sites; site += step) {
        execute_upward_op_for_site(D, op, site);
    }
}
__device__ __forceinline__ void execute_downward_op_for_site(
    const DeviceTree& D,
    const NodeOpInfo& op,
    unsigned int site)
{
    switch (op.op_type) {
        case OP_DOWN_INNER_INNER:
            if (D.states == 4) {
                switch (D.rate_cats) {
                    case 1:
                        compute_downward_inner_inner_ratecat<1>(D, op, site);
                        break;
                    case 4:
                        compute_downward_inner_inner_ratecat<4>(D, op, site);
                        break;
                    case 8:
                        compute_downward_inner_inner_ratecat<8>(D, op, site);
                        break;
                    default:
                        compute_downward_inner_inner_generic(D, op, site);
                        break;
                }
            } else {
                compute_downward_inner_inner_generic(D, op, site);
            }
            break;
        case OP_DOWN_INNER_TIP:
            if (D.states == 4) {
                switch (D.rate_cats) {
                    case 1:
                        compute_downward_inner_tip_ratecat<1>(D, op, site);
                        break;
                    case 4:
                        compute_downward_inner_tip_ratecat<4>(D, op, site);
                        break;
                    case 8:
                        compute_downward_inner_tip_ratecat<8>(D, op, site);
                        break;
                    default:
                        compute_downward_inner_tip_generic(D, op, site);
                        break;
                }
            } else {
                compute_downward_inner_tip_generic(D, op, site);
            }
            break;
        case OP_DOWN_TIP_INNER:
            if (D.states == 4) {
                switch (D.rate_cats) {
                    case 1:
                        compute_downward_tip_inner_ratecat<1>(D, op, site);
                        break;
                    case 4:
                        compute_downward_tip_inner_ratecat<4>(D, op, site);
                        break;
                    case 8:
                        compute_downward_tip_inner_ratecat<8>(D, op, site);
                        break;
                    default:
                        compute_downward_tip_inner_generic(D, op, site);
                        break;
                }
            } else {
                compute_downward_tip_inner_generic(D, op, site);
            }
            break;
        case OP_DOWN_TIP_TIP:
            if (D.states == 4) {
                switch (D.rate_cats) {
                    case 1:
                        compute_downward_tip_tip_ratecat<1>(D, op, site);
                        break;
                    case 4:
                        compute_downward_tip_tip_ratecat<4>(D, op, site);
                        break;
                    case 8:
                        compute_downward_tip_tip_ratecat<8>(D, op, site);
                        break;
                    default:
                        compute_downward_tip_tip_generic(D, op, site);
                        break;
                }
            } else {
                compute_downward_tip_tip_generic(D, op, site);
            }
            break;
        default:
            break;
    }
}

// The serial downward kernel preserves preorder dependency by applying the
// complete operation list for one site before advancing to another site.
__global__ void UpdatePartialsDownwardKernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops)
{
    unsigned int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int step = blockDim.x * gridDim.x;

    for (unsigned int site = tid; site < D.sites; site += step) {
        for (int i = 0; i < num_ops; ++i) {
            execute_downward_op_for_site(D, ops[i], site);
        }
    }
}

__global__ void UpdatePartialsDownwardLevelKernel(
    const DeviceTree D,
    const NodeOpInfo* ops,
    int num_ops)
{
    if (blockIdx.y >= static_cast<unsigned int>(num_ops)) {
        return;
    }

    const NodeOpInfo op = ops[blockIdx.y];
    unsigned int tid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int step = blockDim.x * gridDim.x;

    for (unsigned int site = tid; site < D.sites; site += step) {
        execute_downward_op_for_site(D, op, site);
    }
}

__device__ __forceinline__ fp_t warp_site_rate_max(
    fp_t value,
    unsigned int lane,
    unsigned int rate_mask)
{
    value = fp_fmax(
        value, __shfl_down_sync(rate_mask, value, 2, 4));
    value = fp_fmax(
        value, __shfl_down_sync(rate_mask, value, 1, 4));
    return __shfl_sync(rate_mask, value, lane & ~3u);
}

// The warp kernels specialize the common DNA4/four-rate layout. Each 16-lane
// group owns one site, with four consecutive state lanes per rate category;
// reductions and scaler updates must therefore remain confined to each group.
__global__ void BuildTreeEdgeOutsideWarpSiteKernel(
    const DeviceTree D,
    const NodeOpInfo* op_ptr)
{
    constexpr int kComponents = 16;
    if (!op_ptr || D.states != 4 || D.rate_cats != 4 ||
        !D.d_clv_down || !D.d_clv_up || !D.d_edge_outside_clv ||
        !D.d_pmat) {
        return;
    }
    const NodeOpInfo op = op_ptr[0];
    const bool target_is_left =
        op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT);
    const int target_id = target_is_left ? op.left_id : op.right_id;
    const int sibling_id = target_is_left ? op.right_id : op.left_id;
    if (target_id < 0 || sibling_id < 0 || op.parent_id < 0) return;

    extern __shared__ fp_t sibling_matrix[];
    const fp_t* matrix_source = D.d_pmat +
        static_cast<size_t>(sibling_id) * 4 * 16;
    for (int index = threadIdx.x; index < 64; index += blockDim.x) {
        sibling_matrix[index] = matrix_source[index];
    }
    __syncthreads();

    const unsigned int warp_lane = threadIdx.x & 31u;
    const unsigned int component_lane = warp_lane & 15u;
    const unsigned int state_lane = component_lane & 3u;
    const unsigned int rate_mask =
        0x0000000fu << (warp_lane & ~3u);
    const int rate = static_cast<int>(component_lane >> 2);
    const int state = static_cast<int>(state_lane);
    const size_t global_thread =
        static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const size_t site_group = global_thread / kComponents;
    const size_t site_step =
        static_cast<size_t>(gridDim.x) * blockDim.x / kComponents;

    for (size_t site = site_group; site < D.sites; site += site_step) {
        const size_t site_offset = site * 16;
        const fp_t* parent_down = down_clv_ptr<const fp_t>(
            D, op.parent_id, site_offset);
        const fp_t* sibling_up = up_clv_ptr<const fp_t>(
            D, sibling_id, site_offset);
        fp_t* edge_outside = edge_outside_clv_ptr<fp_t>(
            D, target_id, site_offset);
        unsigned int* parent_scaler =
            down_scaler_ptr(D, op.parent_id, site);
        unsigned int* sibling_scaler =
            up_scaler_ptr(D, sibling_id, site);
        unsigned int* base_scaler =
            edge_outside_scaler_ptr(D, target_id, site);
        if (!parent_down || !sibling_up || !edge_outside) continue;

        const fp_t* matrix = sibling_matrix + rate * 16 + state * 4;
        const fp_t* sibling = sibling_up + rate * 4;
        fp_t transformed = fp_t(0);
#pragma unroll
        for (int child_state = 0; child_state < 4; ++child_state) {
            transformed = fp_fma(
                matrix[child_state], sibling[child_state], transformed);
        }
        fp_t value = parent_down[rate * 4 + state] * transformed;
        const fp_t maximum = warp_site_rate_max(
            value, warp_lane, rate_mask);
        const unsigned int local_shift = clv_scale_shift(maximum);
        if (state_lane == 0) {
            write_clv_scaler_shift(
                D, base_scaler, static_cast<unsigned int>(rate),
                read_clv_scaler_shift(
                    D, parent_scaler, static_cast<unsigned int>(rate)) +
                read_clv_scaler_shift(
                    D, sibling_scaler, static_cast<unsigned int>(rate)) +
                local_shift);
        }
        if (local_shift) fp_scale_pow2(value, local_shift);
        edge_outside[rate * 4 + state] = value;
    }
}

__global__ void UpdateTreeUpwardWarpSiteKernel(
    const DeviceTree D,
    const NodeOpInfo* op_ptr)
{
    constexpr int kComponents = 16;
    if (!op_ptr || D.states != 4 || D.rate_cats != 4 ||
        !D.d_clv_up || !D.d_pmat) {
        return;
    }
    const NodeOpInfo op = op_ptr[0];
    if (op.parent_id < 0 || op.left_id < 0 || op.right_id < 0) return;

    extern __shared__ fp_t matrices[];
    fp_t* left_matrix = matrices;
    fp_t* right_matrix = matrices + 64;
    const fp_t* left_source = D.d_pmat +
        static_cast<size_t>(op.left_id) * 4 * 16;
    const fp_t* right_source = D.d_pmat +
        static_cast<size_t>(op.right_id) * 4 * 16;
    for (int index = threadIdx.x; index < 64; index += blockDim.x) {
        left_matrix[index] = left_source[index];
        right_matrix[index] = right_source[index];
    }
    __syncthreads();

    const unsigned int warp_lane = threadIdx.x & 31u;
    const unsigned int component_lane = warp_lane & 15u;
    const unsigned int state_lane = component_lane & 3u;
    const unsigned int rate_mask =
        0x0000000fu << (warp_lane & ~3u);
    const int rate = static_cast<int>(component_lane >> 2);
    const int state = static_cast<int>(state_lane);
    const size_t global_thread =
        static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const size_t site_group = global_thread / kComponents;
    const size_t site_step =
        static_cast<size_t>(gridDim.x) * blockDim.x / kComponents;

    for (size_t site = site_group; site < D.sites; site += site_step) {
        const size_t site_offset = site * 16;
        const fp_t* left_up = up_clv_ptr<const fp_t>(
            D, op.left_id, site_offset);
        const fp_t* right_up = up_clv_ptr<const fp_t>(
            D, op.right_id, site_offset);
        fp_t* parent_up = up_clv_ptr<fp_t>(
            D, op.parent_id, site_offset);
        unsigned int* left_scaler =
            up_scaler_ptr(D, op.left_id, site);
        unsigned int* right_scaler =
            up_scaler_ptr(D, op.right_id, site);
        unsigned int* parent_scaler =
            up_scaler_ptr(D, op.parent_id, site);
        if (!left_up || !right_up || !parent_up) continue;

        const fp_t* left_row = left_matrix + rate * 16 + state * 4;
        const fp_t* right_row = right_matrix + rate * 16 + state * 4;
        fp_t left_value = fp_t(0);
        fp_t right_value = fp_t(0);
#pragma unroll
        for (int child_state = 0; child_state < 4; ++child_state) {
            left_value = fp_fma(
                left_row[child_state],
                left_up[rate * 4 + child_state], left_value);
            right_value = fp_fma(
                right_row[child_state],
                right_up[rate * 4 + child_state], right_value);
        }
        fp_t value = left_value * right_value;
        const fp_t maximum = warp_site_rate_max(
            value, warp_lane, rate_mask);
        const unsigned int local_shift = clv_scale_shift(maximum);
        if (state_lane == 0) {
            write_clv_scaler_shift(
                D, parent_scaler, static_cast<unsigned int>(rate),
                read_clv_scaler_shift(
                    D, left_scaler, static_cast<unsigned int>(rate)) +
                read_clv_scaler_shift(
                    D, right_scaler, static_cast<unsigned int>(rate)) +
                local_shift);
        }
        if (local_shift) fp_scale_pow2(value, local_shift);
        parent_up[rate * 4 + state] = value;
    }
}

__global__ void RefreshTreeChildDownWarpSiteKernel(
    const DeviceTree D,
    int target_id)
{
    constexpr int kComponents = 16;
    if (target_id < 0 || target_id >= D.N || target_id == D.root_id ||
        D.states != 4 || D.rate_cats != 4 || !D.d_pmat ||
        !D.d_clv_down || !D.d_edge_outside_clv) {
        return;
    }
    extern __shared__ fp_t target_matrix[];
    const fp_t* matrix_source = D.d_pmat +
        static_cast<size_t>(target_id) * 4 * 16;
    for (int index = threadIdx.x; index < 64; index += blockDim.x) {
        target_matrix[index] = matrix_source[index];
    }
    __syncthreads();

    const unsigned int warp_lane = threadIdx.x & 31u;
    const unsigned int component_lane = warp_lane & 15u;
    const unsigned int state_lane = component_lane & 3u;
    const unsigned int rate_mask =
        0x0000000fu << (warp_lane & ~3u);
    const int rate = static_cast<int>(component_lane >> 2);
    const int state = static_cast<int>(state_lane);
    const size_t global_thread =
        static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const size_t site_group = global_thread / kComponents;
    const size_t site_step =
        static_cast<size_t>(gridDim.x) * blockDim.x / kComponents;

    for (size_t site = site_group; site < D.sites; site += site_step) {
        const size_t site_offset = site * 16;
        const fp_t* edge_outside = edge_outside_clv_ptr<const fp_t>(
            D, target_id, site_offset);
        fp_t* child_down = down_clv_ptr<fp_t>(
            D, target_id, site_offset);
        unsigned int* base_scaler =
            edge_outside_scaler_ptr(D, target_id, site);
        unsigned int* down_scaler =
            down_scaler_ptr(D, target_id, site);
        if (!edge_outside || !child_down) continue;

        fp_t value = fp_t(0);
#pragma unroll
        for (int parent_state = 0; parent_state < 4; ++parent_state) {
            const size_t matrix_index =
                D.downward_pmat_indexing == DownwardPmatIndexing::Rows
                ? static_cast<size_t>(state) * 4 + parent_state
                : static_cast<size_t>(parent_state) * 4 + state;
            value = fp_fma(
                target_matrix[rate * 16 + matrix_index],
                edge_outside[rate * 4 + parent_state], value);
        }
        const fp_t maximum = warp_site_rate_max(
            value, warp_lane, rate_mask);
        const unsigned int local_shift = clv_scale_shift(maximum);
        if (state_lane == 0) {
            write_clv_scaler_shift(
                D, down_scaler, static_cast<unsigned int>(rate),
                read_clv_scaler_shift(
                    D, base_scaler, static_cast<unsigned int>(rate)) +
                local_shift);
        }
        if (local_shift) fp_scale_pow2(value, local_shift);
        child_down[rate * 4 + state] = value;
    }
}

} // namespace mlipper::likelihood::partials
