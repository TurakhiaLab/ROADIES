#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <cuda_runtime.h>

#include "derivative.cuh"
#include "gpu/device_buffer.hpp"
#include "likelihood/placement_likelihood.cuh"
#include "placement.cuh"
#include "pmatrix/pmat.h"
#include "pmatrix/pmat_gpu.cuh"
#include "tree/tree.hpp"
#include "util/checked_size.hpp"
#include "util/mlipper_util.h"

namespace {
constexpr int kPlacementStateCount = 4;
constexpr int kMaximumPlacementRateCategories = 8;

template <typename Likelihood>
static std::vector<double> compute_like_weight_ratios(
    const std::vector<Likelihood>& log_likelihoods)
{
    std::vector<double> ratios(log_likelihoods.size(), 0.0);
    if (log_likelihoods.empty()) {
        return ratios;
    }

    // EPA-ng-compatible, numerically stable likelihood weight ratios:
    // exp(logL_i - max(logL)) / sum_j exp(logL_j - max(logL)).
    double max_ll = -std::numeric_limits<double>::infinity();
    for (const Likelihood value : log_likelihoods) {
        const double likelihood = static_cast<double>(value);
        if (!std::isfinite(likelihood)) {
            throw std::runtime_error("Placement ranking contains a non-finite likelihood.");
        }
        max_ll = std::max(max_ll, likelihood);
    }
    double sum_weights = 0.0;
    for (size_t i = 0; i < log_likelihoods.size(); ++i) {
        const double weight = std::exp(
            static_cast<double>(log_likelihoods[i]) - max_ll);
        ratios[i] = weight;
        sum_weights += weight;
    }
    if (sum_weights <= 0.0 || !std::isfinite(sum_weights)) {
        ratios.assign(log_likelihoods.size(), 0.0);
        ratios.front() = 1.0;
        return ratios;
    }
    for (double& value : ratios) {
        value /= sum_weights;
    }
    return ratios;
}

static void retain_accumulated_lwr(
    std::vector<RawPlacementResult::RankedPlacement>& placements,
    double threshold)
{
    if (placements.empty() || threshold <= 0.0) {
        return;
    }

    double accumulated_lwr = 0.0;
    size_t keep = 0;
    while (keep < placements.size() && accumulated_lwr < threshold) {
        accumulated_lwr += placements[keep].like_weight_ratio;
        ++keep;
    }
    placements.resize(std::max<size_t>(keep, 1));
}

struct LocalChildRefineFamilyOps {
    int selected_op = -1;
    int child_left_op = -1;
    int child_right_op = -1;
};

static std::vector<NodeOpInfo> load_host_ops_for_local_child_refine(
    const NodeOpInfo* d_ops,
    int num_ops,
    cudaStream_t stream)
{
    std::vector<NodeOpInfo> host_ops;
    if (!d_ops || num_ops <= 0) {
        return host_ops;
    }
    host_ops.resize(static_cast<size_t>(num_ops));
    CUDA_CHECK(cudaMemcpyAsync(
        host_ops.data(),
        d_ops,
        sizeof(NodeOpInfo) * static_cast<size_t>(num_ops),
        cudaMemcpyDeviceToHost,
        stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    return host_ops;
}

static LocalChildRefineFamilyOps find_local_child_refine_family_ops(
    const std::vector<NodeOpInfo>& host_ops,
    int selected_target_id)
{
    LocalChildRefineFamilyOps family;
    if (selected_target_id < 0) {
        return family;
    }
    for (size_t op_idx = 0; op_idx < host_ops.size(); ++op_idx) {
        const NodeOpInfo& op = host_ops[op_idx];
        const int target_id = node_op_target_id(op);
        if (target_id == selected_target_id && family.selected_op < 0) {
            family.selected_op = static_cast<int>(op_idx);
        }
        if (op.parent_id != selected_target_id) continue;
        if (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT)) {
            family.child_left_op = static_cast<int>(op_idx);
        } else if (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_RIGHT)) {
            family.child_right_op = static_cast<int>(op_idx);
        }
    }
    return family;
}

static unsigned int host_placement_scaler_shift_at(
    const std::vector<unsigned>& scaler_slice,
    const DeviceTree& D,
    size_t site_idx,
    size_t rate_idx)
{
    if (scaler_slice.empty()) return 0u;
    if (D.per_rate_scaling) {
        const size_t rate_count = static_cast<size_t>(D.rate_cats);
        return scaler_slice[site_idx * rate_count + rate_idx];
    }
    return scaler_slice[site_idx];
}

static unsigned host_pattern_weight_at(
    const std::vector<unsigned>& pattern_weights,
    size_t site_idx)
{
    return pattern_weights.empty() ? 1u : pattern_weights[site_idx];
}

struct HostPlacementEvalInputs {
    std::vector<fp_t> query_clv;
    std::vector<fp_t> rate_weights;
    std::vector<fp_t> frequencies;
    std::vector<unsigned> pattern_weights;
};

struct HostPlacementPostprocessCache {
    HostPlacementEvalInputs eval_inputs;
    bool eval_inputs_loaded = false;
    std::vector<NodeOpInfo> host_ops;
    bool host_ops_loaded = false;
};

static HostPlacementEvalInputs load_host_placement_eval_inputs(
    const DeviceTree& D,
    cudaStream_t stream)
{
    HostPlacementEvalInputs out;
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t per_site = mlipper::util::checked_mul_size(
        rate_count, state_count, "host placement CLV site stride");

    CUDA_CHECK(cudaStreamSynchronize(stream));

    out.query_clv.resize(mlipper::util::checked_mul_size(
        D.sites, per_site, "host placement query CLV"), fp_t(0));
    out.rate_weights.resize(rate_count, fp_t(0));
    out.frequencies.resize(state_count, fp_t(0));
    out.pattern_weights.assign(D.sites, 1u);

    if (!out.query_clv.empty()) {
        CUDA_CHECK(cudaMemcpy(
            out.query_clv.data(),
            D.d_query_clv,
            sizeof(fp_t) * out.query_clv.size(),
            cudaMemcpyDeviceToHost));
    }
    CUDA_CHECK(cudaMemcpy(
        out.rate_weights.data(),
        D.d_rate_weights,
        sizeof(fp_t) * out.rate_weights.size(),
        cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(
        out.frequencies.data(),
        D.d_frequencies,
        sizeof(fp_t) * out.frequencies.size(),
        cudaMemcpyDeviceToHost));
    if (D.d_pattern_weights_u && D.sites > 0) {
        CUDA_CHECK(cudaMemcpy(
            out.pattern_weights.data(),
            D.d_pattern_weights_u,
            sizeof(unsigned) * D.sites,
            cudaMemcpyDeviceToHost));
    }

    return out;
}

static const HostPlacementEvalInputs& ensure_host_placement_eval_inputs(
    const DeviceTree& D,
    cudaStream_t stream,
    HostPlacementPostprocessCache& cache)
{
    if (!cache.eval_inputs_loaded) {
        cache.eval_inputs = load_host_placement_eval_inputs(D, stream);
        cache.eval_inputs_loaded = true;
    }
    return cache.eval_inputs;
}

static const std::vector<NodeOpInfo>& ensure_host_ops_loaded(
    const NodeOpInfo* d_ops,
    int num_ops,
    cudaStream_t stream,
    HostPlacementPostprocessCache& cache)
{
    if (!cache.host_ops_loaded) {
        cache.host_ops = load_host_ops_for_local_child_refine(d_ops, num_ops, stream);
        cache.host_ops_loaded = true;
    }
    return cache.host_ops;
}

struct DoubleRerankCandidateBuffers {
    std::vector<fp_t> pendant_pmat;
    std::vector<fp_t> distal_pmat;
    std::vector<fp_t> proximal_pmat;
    std::vector<fp_t> distal_clv;
    std::vector<fp_t> proximal_clv;
    std::vector<unsigned> distal_scalers;
    std::vector<unsigned> proximal_scalers;
};

static DoubleRerankCandidateBuffers load_double_rerank_candidate_buffers(
    const DeviceTree& D,
    int op_index,
    int target_id)
{
    DoubleRerankCandidateBuffers out;
    const size_t per_query = D.pmat_per_node_elems();
    const size_t per_node_pmat = D.pmat_per_node_elems();
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t op_offset = static_cast<size_t>(op_index);
    const size_t target_offset = static_cast<size_t>(target_id);
    const size_t per_site = rate_count * state_count;
    const size_t per_node_clv = D.sites * per_site;
    const size_t scaler_span = D.per_rate_scaling
        ? (D.sites * rate_count)
        : D.sites;

    out.pendant_pmat.resize(per_query);
    out.distal_pmat.resize(per_node_pmat);
    out.proximal_pmat.resize(per_node_pmat);
    out.distal_clv.resize(per_node_clv);
    out.proximal_clv.resize(per_node_clv);
    if (scaler_span > 0) {
        out.distal_scalers.resize(scaler_span);
        out.proximal_scalers.resize(scaler_span);
    }

    CUDA_CHECK(cudaMemcpy(
        out.pendant_pmat.data(),
        D.d_query_pmat + op_offset * per_query,
        sizeof(fp_t) * per_query,
        cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(
        out.distal_pmat.data(),
        D.d_pmat_mid_dist + target_offset * per_node_pmat,
        sizeof(fp_t) * per_node_pmat,
        cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(
        out.proximal_pmat.data(),
        D.d_pmat_mid_prox + target_offset * per_node_pmat,
        sizeof(fp_t) * per_node_pmat,
        cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(
        out.distal_clv.data(),
        D.d_edge_outside_clv + target_offset * per_node_clv,
        sizeof(fp_t) * per_node_clv,
        cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(
        out.proximal_clv.data(),
        D.d_clv_up + target_offset * per_node_clv,
        sizeof(fp_t) * per_node_clv,
        cudaMemcpyDeviceToHost));
    if (scaler_span > 0) {
        if (D.d_edge_outside_scaler) {
            CUDA_CHECK(cudaMemcpy(
                out.distal_scalers.data(),
                D.d_edge_outside_scaler + target_offset * scaler_span,
                sizeof(unsigned) * scaler_span,
                cudaMemcpyDeviceToHost));
        }
        if (D.d_site_scaler_up) {
            CUDA_CHECK(cudaMemcpy(
                out.proximal_scalers.data(),
                D.d_site_scaler_up + target_offset * scaler_span,
                sizeof(unsigned) * scaler_span,
                cudaMemcpyDeviceToHost));
        }
    }
    return out;
}

static double recompute_candidate_loglikelihood_double(
    const DeviceTree& D,
    const HostPlacementEvalInputs& host_inputs,
    const DoubleRerankCandidateBuffers& candidate)
{
    const int state_count = D.states;
    const int rate_cat_count = D.rate_cats;
    const size_t states = static_cast<size_t>(state_count);
    const size_t rate_cats = static_cast<size_t>(rate_cat_count);
    const size_t matrix_elems = states * states;
    const size_t per_site = rate_cats * states;
    const double eps = 1e-300;
    constexpr double kLn2Host = 0.69314718055994530942;
    std::vector<double> rate_vals(rate_cats, 0.0);
    std::vector<unsigned int> rate_shifts(rate_cats, 0u);

    double total = 0.0;
    for (size_t site = 0; site < D.sites; ++site) {
        const size_t site_off = site * per_site;
        std::fill(rate_vals.begin(), rate_vals.end(), 0.0);
        std::fill(rate_shifts.begin(), rate_shifts.end(), 0u);
        unsigned int site_min_shift = 0u;
        bool have_positive = false;

        for (int rc = 0; rc < rate_cat_count; ++rc) {
            const size_t rate = static_cast<size_t>(rc);
            const size_t rc_off = rate * states;
            const size_t matrix_off = rate * matrix_elems;
            double rate_sum = 0.0;
            for (int s = 0; s < state_count; ++s) {
                const size_t state = static_cast<size_t>(s);
                double acc_pend = 0.0;
                double acc_dist = 0.0;
                double acc_prox = 0.0;
                const size_t row_off = matrix_off + state * states;
                for (int k = 0; k < state_count; ++k) {
                    const size_t column = static_cast<size_t>(k);
                    const size_t idx = row_off + column;
                    acc_pend += static_cast<double>(candidate.pendant_pmat[idx]) *
                                static_cast<double>(host_inputs.query_clv[site_off + rc_off + column]);
                    acc_dist += static_cast<double>(candidate.distal_pmat[idx]) *
                                static_cast<double>(candidate.distal_clv[site_off + rc_off + column]);
                    acc_prox += static_cast<double>(candidate.proximal_pmat[idx]) *
                                static_cast<double>(candidate.proximal_clv[site_off + rc_off + column]);
                }
                rate_sum += acc_pend * acc_dist * acc_prox *
                    static_cast<double>(host_inputs.frequencies[state]);
            }

            const unsigned int distal_shift =
                host_placement_scaler_shift_at(candidate.distal_scalers, D, site, rate);
            const unsigned int prox_shift =
                host_placement_scaler_shift_at(candidate.proximal_scalers, D, site, rate);
            const unsigned int total_shift = distal_shift + prox_shift;
            rate_vals[rate] = rate_sum;
            rate_shifts[rate] = total_shift;
            if (rate_sum > 0.0) {
                if (!have_positive || total_shift < site_min_shift) {
                    site_min_shift = total_shift;
                }
                have_positive = true;
            }
        }

        double site_lk = 0.0;
        for (int rc = 0; rc < rate_cat_count; ++rc) {
            const size_t rate = static_cast<size_t>(rc);
            double val = rate_vals[rate];
            if (val > 0.0) {
                const unsigned int diff = rate_shifts[rate] - site_min_shift;
                if (diff) val = std::ldexp(val, -static_cast<int>(diff));
                site_lk += static_cast<double>(host_inputs.rate_weights[rate]) * val;
            }
        }

        total += static_cast<double>(host_pattern_weight_at(host_inputs.pattern_weights, site)) *
            (std::log(site_lk > eps ? site_lk : eps) - static_cast<double>(site_min_shift) * kLn2Host);
    }
    return total;
}

#if !defined(MLIPPER_USE_DOUBLE)
static void maybe_apply_double_rerank(
    const DeviceTree& D,
    const NodeOpInfo* d_ops,
    const std::vector<int>& top_indices,
    RawPlacementResult& result,
    const HostPlacementEvalInputs& host_inputs,
    const PlacementTuningConfig& tuning)
{
    if (!tuning.enable_double_rerank) return;
    if (!d_ops) return;
    if (top_indices.size() < 2 || result.top_placements.size() < 2) return;

    const double best_ll = result.top_placements.front().loglikelihood;
    const double gap_top2 = result.top_placements[0].loglikelihood - result.top_placements[1].loglikelihood;
    const int ulp_factor = std::max(0, tuning.double_rerank_ulp_factor);
    const double gap_floor = std::max(0.0, tuning.double_rerank_gap_top2);
    const float best_ll_f = static_cast<float>(best_ll);
    const float best_ll_next =
        std::nextafter(best_ll_f, std::numeric_limits<float>::infinity());
    const double ulp_gap =
        std::fabs(static_cast<double>(best_ll_next) - static_cast<double>(best_ll_f)) *
        static_cast<double>(ulp_factor);
    const double trigger_gap = std::max(gap_floor, ulp_gap);
    if (!(gap_top2 <= trigger_gap)) return;

    size_t rerank_count = 1;
    while (rerank_count < result.top_placements.size()) {
        const double gap = result.top_placements[0].loglikelihood - result.top_placements[rerank_count].loglikelihood;
        if (gap > trigger_gap) break;
        ++rerank_count;
    }
    if (rerank_count < 2) return;

    struct RankedWithOriginal {
        RawPlacementResult::RankedPlacement placement;
        size_t original_rank = 0;
    };

    std::vector<RankedWithOriginal> reranked;
    reranked.reserve(result.top_placements.size());
    for (size_t i = 0; i < result.top_placements.size(); ++i) {
        reranked.push_back(RankedWithOriginal{result.top_placements[i], i});
    }

    for (size_t i = 0; i < rerank_count; ++i) {
        const int op_index = top_indices[i];
        if (op_index < 0) continue;

        NodeOpInfo host_op{};
        CUDA_CHECK(cudaMemcpy(
            &host_op,
            d_ops + op_index,
            sizeof(NodeOpInfo),
            cudaMemcpyDeviceToHost));

        const int target_id = node_op_target_id(host_op);
        if (target_id < 0 || target_id >= D.N) continue;

        const DoubleRerankCandidateBuffers candidate =
            load_double_rerank_candidate_buffers(D, op_index, target_id);
        reranked[i].placement.loglikelihood = recompute_candidate_loglikelihood_double(
            D,
            host_inputs,
            candidate);
    }

    std::stable_sort(
        reranked.begin(),
        reranked.end(),
        [](const RankedWithOriginal& lhs, const RankedWithOriginal& rhs) {
            if (lhs.placement.loglikelihood != rhs.placement.loglikelihood) {
                return lhs.placement.loglikelihood > rhs.placement.loglikelihood;
            }
            return lhs.original_rank < rhs.original_rank;
        });

    result.top_placements.clear();
    result.top_placements.reserve(reranked.size());
    std::vector<double> reranked_logliks;
    reranked_logliks.reserve(reranked.size());
    for (const RankedWithOriginal& entry : reranked) {
        result.top_placements.push_back(entry.placement);
        reranked_logliks.push_back(entry.placement.loglikelihood);
    }

    const std::vector<double> like_weight_ratios = compute_like_weight_ratios(reranked_logliks);
    for (size_t i = 0; i < result.top_placements.size() && i < like_weight_ratios.size(); ++i) {
        result.top_placements[i].like_weight_ratio = like_weight_ratios[i];
    }

    result.target_id = result.top_placements.front().target_id;
    result.loglikelihood = result.top_placements.front().loglikelihood;
    result.distal_length = result.top_placements.front().distal_length;
    result.pendant_length = result.top_placements.front().pendant_length;
}

#endif

static void rerank_selected_target_and_children(
    const DeviceTree& D,
    RawPlacementResult& result,
    const std::vector<NodeOpInfo>& host_ops,
    const HostPlacementEvalInputs& host_inputs,
    const PlacementTuningConfig& tuning)
{
    if (host_ops.empty()) return;
    if (result.target_id < 0 || result.target_id >= D.N) return;
    const LocalChildRefineFamilyOps family =
        find_local_child_refine_family_ops(host_ops, result.target_id);
    if (family.selected_op < 0) return;
    if (family.child_left_op < 0 && family.child_right_op < 0) return;

    struct LocalCandidate {
        RawPlacementResult::RankedPlacement placement;
        int op_index = -1;
    };

    std::vector<LocalCandidate> local_candidates;
    local_candidates.reserve(3);
    auto append_local_candidate = [&](int op_index) {
        if (op_index < 0) return;
        const NodeOpInfo& op = host_ops[static_cast<size_t>(op_index)];
        const int target_id = node_op_target_id(op);
        if (target_id < 0 || target_id >= D.N) return;

        LocalCandidate candidate;
        candidate.op_index = op_index;
        candidate.placement.target_id = target_id;

        fp_t pendant_length = fp_t(0);
        fp_t proximal_length = fp_t(0);
        CUDA_CHECK(cudaMemcpy(
            &pendant_length,
            D.d_prev_pendant_length + target_id,
            sizeof(fp_t),
            cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(
            &proximal_length,
            D.d_prev_proximal_length + target_id,
            sizeof(fp_t),
            cudaMemcpyDeviceToHost));
        candidate.placement.pendant_length = static_cast<double>(pendant_length);
        candidate.placement.distal_length = static_cast<double>(proximal_length);

        const DoubleRerankCandidateBuffers buffers =
            load_double_rerank_candidate_buffers(D, op_index, target_id);
        candidate.placement.loglikelihood = recompute_candidate_loglikelihood_double(
            D,
            host_inputs,
            buffers);
        local_candidates.push_back(candidate);
    };

    append_local_candidate(family.selected_op);
    append_local_candidate(family.child_left_op);
    append_local_candidate(family.child_right_op);
    if (local_candidates.size() < 2) return;

    std::stable_sort(
        local_candidates.begin(),
        local_candidates.end(),
        [](const LocalCandidate& lhs, const LocalCandidate& rhs) {
            if (lhs.placement.loglikelihood != rhs.placement.loglikelihood) {
                return lhs.placement.loglikelihood > rhs.placement.loglikelihood;
            }
            return lhs.op_index < rhs.op_index;
        });

    std::vector<RawPlacementResult::RankedPlacement> merged = result.top_placements;
    for (const LocalCandidate& local : local_candidates) {
        bool replaced = false;
        for (RawPlacementResult::RankedPlacement& existing : merged) {
            if (existing.target_id == local.placement.target_id) {
                existing = local.placement;
                replaced = true;
                break;
            }
        }
        if (!replaced) {
            merged.push_back(local.placement);
        }
    }

    auto existing_rank = [&](int target_id) -> size_t {
        for (size_t rank = 0; rank < result.top_placements.size(); ++rank) {
            if (result.top_placements[rank].target_id == target_id) return rank;
        }
        return result.top_placements.size();
    };

    std::stable_sort(
        merged.begin(),
        merged.end(),
        [&](const RawPlacementResult::RankedPlacement& lhs, const RawPlacementResult::RankedPlacement& rhs) {
            if (lhs.loglikelihood != rhs.loglikelihood) {
                return lhs.loglikelihood > rhs.loglikelihood;
            }
            return existing_rank(lhs.target_id) < existing_rank(rhs.target_id);
        });

    const size_t keep = tuning.export_accumulated_lwr_threshold > 0.0
        ? merged.size()
        : std::max<size_t>(
            static_cast<size_t>(std::max(1, tuning.export_placement_topk)),
            3);
    if (merged.size() > keep) {
        merged.resize(keep);
    }
    std::vector<fp_t> merged_logliks;
    merged_logliks.reserve(merged.size());
    for (const RawPlacementResult::RankedPlacement& placement : merged) {
        merged_logliks.push_back(static_cast<fp_t>(placement.loglikelihood));
    }
    const std::vector<double> like_weight_ratios = compute_like_weight_ratios(merged_logliks);
    for (size_t i = 0; i < merged.size() && i < like_weight_ratios.size(); ++i) {
        merged[i].like_weight_ratio = like_weight_ratios[i];
    }

    result.top_placements.swap(merged);
    result.target_id = result.top_placements.front().target_id;
    result.loglikelihood = result.top_placements.front().loglikelihood;
    result.distal_length = result.top_placements.front().distal_length;
    result.pendant_length = result.top_placements.front().pendant_length;
}

struct RankedPlacementDeviceMetadata {
    int target_id = -1;
    fp_t pendant_length = fp_t(0);
    fp_t proximal_length = fp_t(0);
};

struct PlacementKernelScratchBuffers {
    mlipper::gpu::DeviceBuffer<fp_t> previous_loglikelihoods;
    mlipper::gpu::DeviceBuffer<int> active_operations;
    mlipper::gpu::DeviceBuffer<int> ranked_indices;
    mlipper::gpu::DeviceBuffer<RankedPlacementDeviceMetadata> ranked_metadata;
    std::vector<fp_t> host_loglk_cache;
    std::vector<int> host_order_cache;
    std::vector<RankedPlacementDeviceMetadata> host_ranked_metadata;
};

__global__ void GatherRankedPlacementMetadataKernel(
    DeviceTree D,
    const NodeOpInfo* ops,
    const int* ranked_indices,
    int count,
    RankedPlacementDeviceMetadata* metadata)
{
    const int rank = blockIdx.x * blockDim.x + threadIdx.x;
    if (rank >= count || !ops || !ranked_indices || !metadata) return;
    const int op_index = ranked_indices[rank];
    if (op_index < 0) return;
    const NodeOpInfo op = ops[op_index];
    const int target_id = node_op_target_id(op);
    if (target_id < 0 || target_id >= D.N) return;
    metadata[rank].target_id = target_id;
    metadata[rank].pendant_length = D.d_prev_pendant_length[target_id];
    metadata[rank].proximal_length = D.d_prev_proximal_length[target_id];
}

static std::vector<RawPlacementResult::RankedPlacement> build_top_ranked_placements(
    const DeviceTree& D,
    const NodeOpInfo* d_ops,
    const std::vector<int>& top_indices,
    const std::vector<fp_t>& top_values,
    PlacementKernelScratchBuffers& scratch,
    cudaStream_t stream)
{
    // Only compact metadata crosses back after host ranking; the full
    // per-candidate CLV and PMAT buffers remain resident on the device.
    std::vector<RawPlacementResult::RankedPlacement> ranked;
    const size_t keep = std::min(top_indices.size(), top_values.size());
    ranked.reserve(keep);

    scratch.ranked_indices.ensureCapacity(keep);
    scratch.ranked_metadata.ensureCapacity(keep);
    scratch.host_ranked_metadata.assign(keep, RankedPlacementDeviceMetadata{});
    CUDA_CHECK(cudaMemcpyAsync(
        scratch.ranked_indices.get(), top_indices.data(),
        keep * sizeof(int), cudaMemcpyHostToDevice, stream));
    constexpr int kBlockSize = 256;
    const int grid_size = static_cast<int>((keep + kBlockSize - 1) / kBlockSize);
    GatherRankedPlacementMetadataKernel<<<grid_size, kBlockSize, 0, stream>>>(
        D, d_ops, scratch.ranked_indices.get(), static_cast<int>(keep),
        scratch.ranked_metadata.get());
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(
        scratch.host_ranked_metadata.data(), scratch.ranked_metadata.get(),
        keep * sizeof(RankedPlacementDeviceMetadata),
        cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));

    const std::vector<double> like_weight_ratios = compute_like_weight_ratios(top_values);
    for (size_t i = 0; i < keep; ++i) {
        const int op_index = top_indices[i];
        if (op_index < 0) {
            continue;
        }
        const RankedPlacementDeviceMetadata& metadata =
            scratch.host_ranked_metadata[i];
        const int target_id = metadata.target_id;
        if (target_id < 0 || target_id >= D.N) {
            continue;
        }

        RawPlacementResult::RankedPlacement candidate;
        candidate.target_id = target_id;
        candidate.loglikelihood = static_cast<double>(top_values[i]);
        candidate.distal_length = static_cast<double>(metadata.proximal_length);
        candidate.pendant_length = static_cast<double>(metadata.pendant_length);
        candidate.like_weight_ratio = like_weight_ratios[i];
        ranked.push_back(candidate);
    }
    return ranked;
}

template <typename CheckCudaFn>
static void fetch_ranked_loglikelihoods(
    const fp_t* d_source,
    int num_ops,
    int result_limit,
    PlacementKernelScratchBuffers& scratch,
    cudaStream_t stream,
    std::vector<int>& top_indices,
    std::vector<fp_t>& top_values,
    CheckCudaFn&& check_cuda)
{
    top_indices.clear();
    top_values.clear();
    if (!d_source || num_ops <= 0 || result_limit <= 0) return;

    scratch.host_loglk_cache.resize(static_cast<size_t>(num_ops));
    check_cuda("cudaMemcpyAsync host_loglk_cache", cudaMemcpyAsync(
        scratch.host_loglk_cache.data(),
        d_source,
        sizeof(fp_t) * static_cast<size_t>(num_ops),
        cudaMemcpyDeviceToHost,
        stream));
    check_cuda("cudaStreamSynchronize host_loglk_cache", cudaStreamSynchronize(stream));

    for (const fp_t value : scratch.host_loglk_cache) {
        if (!std::isfinite(static_cast<double>(value))) {
            throw std::runtime_error(
                "Placement likelihood evaluation produced a non-finite score.");
        }
    }

    const int result_count = std::min(num_ops, result_limit);
    scratch.host_order_cache.resize(static_cast<size_t>(num_ops));
    std::iota(scratch.host_order_cache.begin(), scratch.host_order_cache.end(), 0);
    std::partial_sort(
        scratch.host_order_cache.begin(),
        scratch.host_order_cache.begin() + result_count,
        scratch.host_order_cache.end(),
        [&](int lhs, int rhs) {
            return scratch.host_loglk_cache[static_cast<size_t>(lhs)] >
                   scratch.host_loglk_cache[static_cast<size_t>(rhs)];
        });

    top_indices.resize(static_cast<size_t>(result_count));
    top_values.resize(static_cast<size_t>(result_count));
    for (int i = 0; i < result_count; ++i) {
        const int op_idx = scratch.host_order_cache[static_cast<size_t>(i)];
        top_indices[static_cast<size_t>(i)] = op_idx;
        top_values[static_cast<size_t>(i)] =
            scratch.host_loglk_cache[static_cast<size_t>(op_idx)];
    }
}

static void include_best_target_children(
    const DeviceTree& D,
    const std::vector<NodeOpInfo>& host_ops,
    const std::vector<fp_t>& host_loglk_cache,
    int result_limit,
    std::vector<int>& final_top_indices,
    std::vector<fp_t>& final_top_values)
{
    if (host_ops.empty() ||
        final_top_indices.empty() ||
        final_top_values.empty() ||
        host_loglk_cache.size() != host_ops.size()) {
        return;
    }

    const int best_op_index = final_top_indices.front();
    if (best_op_index < 0 || best_op_index >= static_cast<int>(host_ops.size())) {
        return;
    }

    const int best_target_id = node_op_target_id(host_ops[static_cast<size_t>(best_op_index)]);
    if (best_target_id < 0 || best_target_id >= D.N) {
        return;
    }

    const LocalChildRefineFamilyOps family =
        find_local_child_refine_family_ops(host_ops, best_target_id);
    if (family.child_left_op < 0 && family.child_right_op < 0) {
        return;
    }

    std::vector<int> augmented_indices = final_top_indices;
    if (family.child_left_op >= 0) {
        augmented_indices.push_back(family.child_left_op);
    }
    if (family.child_right_op >= 0) {
        augmented_indices.push_back(family.child_right_op);
    }
    std::sort(
        augmented_indices.begin(),
        augmented_indices.end(),
        [&](int lhs, int rhs) {
            const fp_t lhs_ll = host_loglk_cache[static_cast<size_t>(lhs)];
            const fp_t rhs_ll = host_loglk_cache[static_cast<size_t>(rhs)];
            if (lhs_ll == rhs_ll) return lhs < rhs;
            return lhs_ll > rhs_ll;
        });
    augmented_indices.erase(
        std::unique(augmented_indices.begin(), augmented_indices.end()),
        augmented_indices.end());

    const int keep = std::min<int>(
        result_limit, static_cast<int>(augmented_indices.size()));
    augmented_indices.resize(static_cast<size_t>(keep));

    std::vector<fp_t> augmented_values(static_cast<size_t>(keep), fp_t(0));
    for (int i = 0; i < keep; ++i) {
        augmented_values[static_cast<size_t>(i)] =
            host_loglk_cache[static_cast<size_t>(augmented_indices[static_cast<size_t>(i)])];
    }

    final_top_indices.swap(augmented_indices);
    final_top_values.swap(augmented_values);
}

}

__global__ void BuildNodePendantLengthsKernel(
    const fp_t* node_lengths,
    fp_t* out_lengths,
    int total_nodes,
    int root_id,
    fp_t min_len,
    fp_t max_len,
    fp_t default_len)
{
    const int node_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (node_id >= total_nodes) return;
    if (!out_lengths) return;

    fp_t branch_length = default_len;
    if (node_lengths) {
        branch_length = node_lengths[node_id];
    }
    if (node_id == root_id) {
        branch_length = default_len;
    }
    if (branch_length < min_len) branch_length = min_len;
    if (branch_length > max_len) branch_length = max_len;
    out_lengths[node_id] = branch_length;
}

__global__ void BuildInitialProximalLengthsKernel(
    const fp_t* node_lengths,
    fp_t* out_lengths,
    int total_nodes,
    int root_id,
    fp_t min_len,
    fp_t max_len,
    fp_t default_len)
{
    const int node_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (node_id >= total_nodes) return;
    if (!out_lengths) return;

    fp_t branch_length = default_len;
    if (node_lengths) {
        branch_length = static_cast<fp_t>(0.5) * node_lengths[node_id];
    }
    if (node_id == root_id) {
        branch_length = default_len;
    }
    if (branch_length < min_len) branch_length = min_len;
    if (branch_length > max_len) branch_length = max_len;
    out_lengths[node_id] = branch_length;
}

// Keep per-op best log-likelihood; rollback branch lengths if current pass is worse.
__global__ void KeepBestBranchLengthsKernel(
    const NodeOpInfo* ops,
    const int* op_indices,
    fp_t* curr_loglk,
    fp_t* prev_loglk,
    fp_t* curr_pendant,
    fp_t* curr_proximal,
    fp_t* prev_pendant,
    fp_t* prev_proximal,
    int* active_ops,
    int num_ops,
    int total_nodes)
{
    const int op_local = blockIdx.x * blockDim.x + threadIdx.x;
    if (op_local >= num_ops) return;
    if (!ops || !curr_loglk || !prev_loglk ||
        !curr_pendant || !curr_proximal || !prev_pendant || !prev_proximal) return;
    const int op_idx = op_indices ? op_indices[op_local] : op_local;
    if (op_idx < 0) return;

    const NodeOpInfo op = ops[op_idx];
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_id = target_is_left ? op.left_id : op.right_id;
    if (target_id < 0 || target_id >= total_nodes) return;
    const fp_t curr = curr_loglk[op_idx];
    const fp_t prev = prev_loglk[op_idx];
    if (!isfinite(static_cast<double>(curr)) ||
        (isfinite(static_cast<double>(prev)) && curr < prev)) {
        curr_loglk[op_idx] = prev;
        curr_pendant[target_id] = prev_pendant[target_id];
        curr_proximal[target_id] = prev_proximal[target_id];
        if (active_ops) active_ops[op_local] = 0;
    } else {
        prev_loglk[op_idx] = curr;
        prev_pendant[target_id] = curr_pendant[target_id];
        prev_proximal[target_id] = curr_proximal[target_id];
    }
}

// Build per-placement pendant PMATs directly from the target branch lengths.
__global__ void BuildPendantPMATPerOpKernel(
    const NodeOpInfo* ops,
    const int* op_indices,
    const fp_t* node_lengths,
    const fp_t* Vinv,
    const fp_t* V,
    const fp_t* lambdas,
    fp_t p,
    fp_t* P,
    int states,
    int rate_cats,
    int num_ops,
    int total_nodes,
    fp_t min_len,
    fp_t max_len,
    fp_t default_len)
{
    const int flat_index = blockIdx.x * blockDim.x + threadIdx.x;
    const int total_entries = num_ops * rate_cats;
    if (flat_index >= total_entries) return;

    const int op_local = flat_index / rate_cats;
    const int rate_idx = flat_index - op_local * rate_cats;
    if (!ops || op_local >= num_ops || rate_idx >= rate_cats) return;
    const int op_idx = op_indices ? op_indices[op_local] : op_local;
    if (op_idx < 0) return;

    const NodeOpInfo op = ops[op_idx];
    const bool target_is_left = (op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const int target_id = target_is_left ? op.left_id : op.right_id;

    fp_t branch_length = default_len;
    if (node_lengths && target_id >= 0 && target_id < total_nodes) {
        branch_length = node_lengths[target_id];
    }
    if (branch_length < min_len) branch_length = min_len;
    if (branch_length > max_len) branch_length = max_len;

    const size_t state_count = static_cast<size_t>(states);
    const size_t rate_offset = static_cast<size_t>(rate_idx) * state_count;
    const size_t matrix_span = state_count * state_count;
    const fp_t* rate_lambdas = lambdas + rate_offset;
    fp_t* out_pmat = P + static_cast<size_t>(flat_index) * matrix_span;
    pmatrix_from_triple_gpu(Vinv, V, rate_lambdas, fp_t(1.0), branch_length, p, out_pmat, states);
}

// Build proximal PMATs for every node from the current proximal branch lengths.
__global__ void BuildNodeProximalPMATKernel(
    const fp_t* node_lengths,
    const fp_t* Vinv,
    const fp_t* V,
    const fp_t* lambdas,
    fp_t p,
    fp_t* P,
    int states,
    int rate_cats,
    int num_nodes,
    int root_id,
    fp_t min_len,
    fp_t max_len,
    fp_t default_len)
{
    const int flat_index = blockIdx.x * blockDim.x + threadIdx.x;
    const int total_entries = num_nodes * rate_cats;
    if (flat_index >= total_entries) return;

    const int node_idx = flat_index / rate_cats;
    const int rate_idx = flat_index - node_idx * rate_cats;
    if (node_idx >= num_nodes || rate_idx >= rate_cats) return;

    fp_t branch_length = default_len;
    if (node_lengths) {
        branch_length = node_lengths[node_idx];
    }
    if (node_idx == root_id) {
        branch_length = default_len;
    }
    if (branch_length < min_len) branch_length = min_len;
    if (branch_length > max_len) branch_length = max_len;

    const size_t state_count = static_cast<size_t>(states);
    const size_t rate_count = static_cast<size_t>(rate_cats);
    const size_t rate_offset = static_cast<size_t>(rate_idx) * state_count;
    const size_t matrix_span = state_count * state_count;
    const size_t output_base =
        (static_cast<size_t>(node_idx) * rate_count + static_cast<size_t>(rate_idx)) * matrix_span;
    const fp_t* rate_lambdas = lambdas + rate_offset;
    fp_t* out_pmat = P + output_base;
    pmatrix_from_triple_gpu(Vinv, V, rate_lambdas, fp_t(1.0), branch_length, p, out_pmat, states);
}

// Build distal PMATs for every node from total branch length minus proximal length.
__global__ void BuildNodeDistalPMATKernel(
    const fp_t* total_lengths,
    const fp_t* proximal_lengths,
    const fp_t* Vinv,
    const fp_t* V,
    const fp_t* lambdas,
    fp_t p,
    fp_t* P,
    int states,
    int rate_cats,
    int num_nodes,
    int root_id,
    fp_t min_len,
    fp_t max_len,
    fp_t default_len)
{
    const int flat_index = blockIdx.x * blockDim.x + threadIdx.x;
    const int total_entries = num_nodes * rate_cats;
    if (flat_index >= total_entries) return;

    const int node_idx = flat_index / rate_cats;
    const int rate_idx = flat_index - node_idx * rate_cats;
    if (node_idx >= num_nodes || rate_idx >= rate_cats) return;
    if (!total_lengths || !proximal_lengths) return;

    fp_t branch_length = total_lengths[node_idx] - proximal_lengths[node_idx];
    if (node_idx == root_id) {
        branch_length = default_len;
    }
    if (branch_length < min_len) branch_length = min_len;
    if (branch_length > max_len) branch_length = max_len;

    const size_t state_count = static_cast<size_t>(states);
    const size_t rate_count = static_cast<size_t>(rate_cats);
    const size_t rate_offset = static_cast<size_t>(rate_idx) * state_count;
    const size_t matrix_span = state_count * state_count;
    const size_t output_base =
        (static_cast<size_t>(node_idx) * rate_count + static_cast<size_t>(rate_idx)) * matrix_span;
    const fp_t* rate_lambdas = lambdas + rate_offset;
    fp_t* out_pmat = P + output_base;
    pmatrix_from_triple_gpu(Vinv, V, rate_lambdas, fp_t(1.0), branch_length, p, out_pmat, states);
}

// Placement is deliberately split into device-side scoring/optimization and
// host-side ranking. The returned node IDs still refer to the unchanged input
// tree; committing a result is the responsibility of the tree workflow layer.
RawPlacementResult EvaluatePlacementCandidates(
    const DeviceTree& D,
    const NodeOpInfo* d_ops,
    int num_ops,
    int smoothing,
    cudaStream_t stream,
    bool enable_local_child_refine,
    const PlacementTuningConfig& tuning,
    const PlacementScratchOverride* scratch_override
)
{
    RawPlacementResult result;
    assert(num_ops > 0 && "num_ops must be positive");
    if (num_ops <= 0) return result;
    assert(smoothing > 0 && "smoothing must be positive");
    if (scratch_override == nullptr) {
        throw std::invalid_argument(
            "EvaluatePlacementCandidates requires explicit scratch buffers.");
    }
    if (D.states != kPlacementStateCount ||
        D.rate_cats <= 0 ||
        D.rate_cats > kMaximumPlacementRateCategories) {
        throw std::invalid_argument(
            "Placement evaluation currently supports four-state models with 1-8 rate categories.");
    }

    // Stage 1: validate inputs and plan the runtime/launch configuration.
    const size_t num_ops_count = static_cast<size_t>(num_ops);
    const unsigned int num_ops_grid_y = static_cast<unsigned int>(num_ops_count);
    const size_t node_count = static_cast<size_t>(D.N);
    const size_t rate_count = static_cast<size_t>(D.rate_cats);
    const size_t state_count = static_cast<size_t>(D.states);
    const size_t sumtable_stride = D.sites * rate_count * state_count;
    if (scratch_override != nullptr) {
        const bool has_proximal_pmat =
            scratch_override->d_pmat_mid_prox != nullptr;
        const bool has_distal_pmat =
            scratch_override->d_pmat_mid_dist != nullptr;
        if (has_proximal_pmat != has_distal_pmat) {
            throw std::runtime_error(
                "EvaluatePlacementCandidates midpoint PMAT scratch is incomplete.");
        }
        const bool has_midpoint_clv =
            scratch_override->d_edge_midpoint_clv != nullptr;
        const bool has_midpoint_scaler =
            scratch_override->d_edge_midpoint_scaler != nullptr;
        if (has_midpoint_clv != has_midpoint_scaler) {
            throw std::runtime_error(
                "EvaluatePlacementCandidates midpoint CLV scratch is incomplete.");
        }
    }
    const size_t sumtable_capacity_ops = scratch_override->sumtable_capacity_ops;
    const size_t likelihood_capacity_ops =
        scratch_override->likelihood_capacity_ops;
    const size_t query_pmat_capacity_ops =
        scratch_override->query_pmat_capacity_ops;
    const size_t midpoint_pmat_capacity_nodes =
        (scratch_override != nullptr &&
         scratch_override->d_pmat_mid_prox != nullptr &&
         scratch_override->d_pmat_mid_dist != nullptr)
            ? scratch_override->midpoint_pmat_capacity_nodes
            : static_cast<size_t>(std::max(D.N, 0));
    const size_t midpoint_clv_capacity_nodes =
        (scratch_override != nullptr &&
         scratch_override->d_edge_midpoint_clv != nullptr &&
         scratch_override->d_edge_midpoint_scaler != nullptr)
            ? scratch_override->midpoint_clv_capacity_nodes
            : static_cast<size_t>(std::max(D.N, 0));
    if (num_ops_count > sumtable_capacity_ops ||
        num_ops_count > likelihood_capacity_ops ||
        num_ops_count > query_pmat_capacity_ops ||
        static_cast<size_t>(std::max(D.N, 0)) > midpoint_pmat_capacity_nodes ||
        static_cast<size_t>(std::max(D.N, 0)) > midpoint_clv_capacity_nodes) {
        throw std::runtime_error(
            "EvaluatePlacementCandidates scratch capacity is too small for " +
            std::to_string(num_ops_count) + " operations and " +
            std::to_string(node_count) + " nodes.");
    }
    auto grid_x = [](size_t work_items, unsigned int block_width) {
        return static_cast<unsigned int>((work_items + block_width - 1) / block_width);
    };
    auto check_launch = [&](const char* stage) {
        const cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string(stage) + ": " + cudaGetErrorString(err));
        }
    };
    auto check_cuda = [&](const char* stage, cudaError_t err) {
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string(stage) + ": " + cudaGetErrorString(err));
        }
    };
    PlacementKernelScratchBuffers scratch;
    fp_t* d_prev_loglk =
        scratch_override != nullptr && scratch_override->d_prev_loglk != nullptr
            ? scratch_override->d_prev_loglk
            : nullptr;
    int* d_active_ops =
        scratch_override != nullptr && scratch_override->d_active_ops != nullptr
            ? scratch_override->d_active_ops
            : nullptr;
    std::vector<fp_t>& host_loglk_cache = scratch.host_loglk_cache;
    HostPlacementPostprocessCache postprocess_cache;
    fp_t* d_likelihoods = scratch_override->d_likelihoods;
    fp_t* d_sumtable = scratch_override->d_sumtable;
    fp_t* d_query_pmat = scratch_override->d_query_pmat;
    fp_t* d_pmat_mid_prox =
        (scratch_override != nullptr && scratch_override->d_pmat_mid_prox != nullptr)
            ? scratch_override->d_pmat_mid_prox
            : D.d_pmat_mid_prox;
    fp_t* d_pmat_mid_dist =
        (scratch_override != nullptr && scratch_override->d_pmat_mid_dist != nullptr)
            ? scratch_override->d_pmat_mid_dist
            : D.d_pmat_mid_dist;
    if (d_likelihoods == nullptr) {
        throw std::runtime_error(
            "EvaluatePlacementCandidates requires a likelihood scratch buffer.");
    }
    if (d_sumtable == nullptr) {
        throw std::runtime_error(
            "EvaluatePlacementCandidates requires a sumtable scratch buffer.");
    }
    if (d_query_pmat == nullptr) {
        throw std::runtime_error(
            "EvaluatePlacementCandidates requires a query PMAT scratch buffer.");
    }
    if (d_pmat_mid_prox == nullptr || d_pmat_mid_dist == nullptr) {
        throw std::runtime_error(
            "EvaluatePlacementCandidates requires midpoint PMAT scratch buffers.");
    }
    // DeviceTree is a non-owning pointer view. Redirect scratch fields in a
    // shallow copy so this call cannot change the caller's buffer bindings.
    DeviceTree placement_view = D;
    placement_view.d_query_pmat = d_query_pmat;
    placement_view.d_pmat_mid_prox = d_pmat_mid_prox;
    placement_view.d_pmat_mid_dist = d_pmat_mid_dist;
    if (scratch_override != nullptr && scratch_override->d_edge_midpoint_clv != nullptr) {
        placement_view.d_edge_midpoint_clv = scratch_override->d_edge_midpoint_clv;
        placement_view.d_edge_midpoint_scaler =
            scratch_override->d_edge_midpoint_scaler;
    }
    if (scratch_override != nullptr &&
        (scratch_override->d_new_pendant_length != nullptr ||
         scratch_override->d_new_proximal_length != nullptr ||
         scratch_override->d_prev_pendant_length != nullptr ||
         scratch_override->d_prev_proximal_length != nullptr)) {
        if (scratch_override->d_new_pendant_length == nullptr ||
            scratch_override->d_new_proximal_length == nullptr ||
            scratch_override->d_prev_pendant_length == nullptr ||
            scratch_override->d_prev_proximal_length == nullptr ||
            scratch_override->branch_length_capacity_nodes < node_count) {
            throw std::runtime_error(
                "EvaluatePlacementCandidates branch-length scratch is incomplete or too small.");
        }
        placement_view.d_new_pendant_length =
            scratch_override->d_new_pendant_length;
        placement_view.d_new_proximal_length =
            scratch_override->d_new_proximal_length;
        placement_view.d_prev_pendant_length =
            scratch_override->d_prev_pendant_length;
        placement_view.d_prev_proximal_length =
            scratch_override->d_prev_proximal_length;
    }
    if (!placement_view.d_new_pendant_length ||
        !placement_view.d_new_proximal_length ||
        !placement_view.d_prev_pendant_length ||
        !placement_view.d_prev_proximal_length) {
        throw std::runtime_error(
            "EvaluatePlacementCandidates requires branch-length scratch buffers.");
    }

    const size_t diag_shared = rate_count * state_count * 4;
    const size_t midpoint_pmat_shared = rate_count * 16 * 2;
    size_t shmem_bytes = sizeof(fp_t) * diag_shared;
    shmem_bytes += sizeof(fp_t) * midpoint_pmat_shared;

    // Pendant and proximal derivative kernels have different register pressure.
    // Size them independently so one kernel does not inherit an invalid launch shape.
    int pendant_block_threads = 512;
    int max_blocks_per_sm = 0;
    cudaFuncAttributes attr{};
    check_cuda(
        "cudaFuncGetAttributes pendant",
        cudaFuncGetAttributes(&attr, LikelihoodDerivativePendantKernel));
    while (pendant_block_threads >= 32) {
        check_cuda(
            "cudaOccupancyMaxActiveBlocksPerMultiprocessor pendant",
            cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                &max_blocks_per_sm,
                LikelihoodDerivativePendantKernel,
                pendant_block_threads,
                shmem_bytes));
        if (max_blocks_per_sm > 0) break;
        pendant_block_threads /= 2;
    }
    if (max_blocks_per_sm == 0) {
        throw std::runtime_error("No valid block size for LikelihoodDerivativePendantKernel on this GPU.");
    }

    int proximal_block_threads = 512;
    max_blocks_per_sm = 0;
    check_cuda(
        "cudaFuncGetAttributes proximal",
        cudaFuncGetAttributes(&attr, LikelihoodDerivativeProximalKernel));
    while (proximal_block_threads >= 32) {
        check_cuda(
            "cudaOccupancyMaxActiveBlocksPerMultiprocessor proximal",
            cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                &max_blocks_per_sm,
                LikelihoodDerivativeProximalKernel,
                proximal_block_threads,
                shmem_bytes));
        if (max_blocks_per_sm > 0) break;
        proximal_block_threads /= 2;
    }
    if (max_blocks_per_sm == 0) {
        throw std::runtime_error("No valid block size for LikelihoodDerivativeProximalKernel on this GPU.");
    }

    const int pmat_block_threads = 128;
    dim3 pendant_block(pendant_block_threads);
    dim3 proximal_block(proximal_block_threads);
    dim3 pmat_block(pmat_block_threads);

    // Stage 2: allocate scratch buffers and initialize branch-length state.
    if (scratch_override != nullptr &&
        (scratch_override->d_prev_loglk != nullptr ||
         scratch_override->d_active_ops != nullptr) &&
        (scratch_override->d_prev_loglk == nullptr ||
         scratch_override->d_active_ops == nullptr ||
         scratch_override->ranking_state_capacity_ops < num_ops_count)) {
        throw std::runtime_error(
            "EvaluatePlacementCandidates ranking-state scratch is incomplete or too small.");
    }
    if (d_prev_loglk == nullptr) {
        scratch.previous_loglikelihoods.ensureCapacity(num_ops_count);
        d_prev_loglk = scratch.previous_loglikelihoods.get();
    }
    if (d_active_ops == nullptr) {
        scratch.active_operations.ensureCapacity(num_ops_count);
        d_active_ops = scratch.active_operations.get();
    }
    check_cuda(
        "cudaMemsetAsync d_active_ops",
        cudaMemsetAsync(d_active_ops, 1, sizeof(int) * num_ops_count, stream));
    {
        dim3 init_block(256);
        dim3 init_grid(grid_x(node_count, init_block.x));
        BuildNodePendantLengthsKernel<<<init_grid, init_block, 0, stream>>>(
            nullptr,
            placement_view.d_prev_pendant_length,
            D.N,
            D.root_id,
            OPT_BRANCH_LEN_MIN,
            OPT_BRANCH_LEN_MAX,
            DEFAULT_BRANCH_LENGTH);
        check_launch("BuildNodePendantLengthsKernel");
        BuildInitialProximalLengthsKernel<<<init_grid, init_block, 0, stream>>>(
            D.d_blen,
            placement_view.d_prev_proximal_length,
            D.N,
            D.root_id,
            OPT_BRANCH_LEN_MIN,
            OPT_BRANCH_LEN_MAX,
            DEFAULT_BRANCH_LENGTH);
        check_launch("BuildInitialProximalLengthsKernel");
    }

    // Stage 3: build the baseline placement state and score every op once.
    {
        const size_t op_rate_work = num_ops_count * rate_count;
        const size_t node_rate_work = node_count * rate_count;
        dim3 pmat_grid(grid_x(op_rate_work, pmat_block.x));
        BuildPendantPMATPerOpKernel<<<pmat_grid, pmat_block, 0, stream>>>(
            d_ops,
            nullptr,
            placement_view.d_prev_pendant_length,
            D.d_Vinv,
            D.d_V,
            D.d_lambdas,
            0.0,
            d_query_pmat,
            D.states,
            D.rate_cats,
            num_ops,
            D.N,
            OPT_BRANCH_LEN_MIN,
            OPT_BRANCH_LEN_MAX,
            DEFAULT_BRANCH_LENGTH);
        check_launch("BuildPendantPMATPerOpKernel baseline");

        dim3 node_grid(grid_x(node_rate_work, pmat_block.x));
        BuildNodeProximalPMATKernel<<<node_grid, pmat_block, 0, stream>>>(
            placement_view.d_prev_proximal_length,
            D.d_Vinv,
            D.d_V,
            D.d_lambdas,
            0.0,
            d_pmat_mid_prox,
            D.states,
            D.rate_cats,
            D.N,
            D.root_id,
            OPT_BRANCH_LEN_MIN,
            OPT_BRANCH_LEN_MAX,
            DEFAULT_BRANCH_LENGTH);
        check_launch("BuildNodeProximalPMATKernel baseline");

        BuildNodeDistalPMATKernel<<<node_grid, pmat_block, 0, stream>>>(
            D.d_blen,
            placement_view.d_prev_proximal_length,
            D.d_Vinv,
            D.d_V,
            D.d_lambdas,
            0.0,
            d_pmat_mid_dist,
            D.states,
            D.rate_cats,
            D.N,
            D.root_id,
            OPT_BRANCH_LEN_MIN,
            OPT_BRANCH_LEN_MAX,
            DEFAULT_BRANCH_LENGTH);
        check_launch("BuildNodeDistalPMATKernel baseline");

        mlipper::likelihood::placement::update_midpoint_partials(
            placement_view, d_ops, num_ops, stream);

        mlipper::likelihood::placement::compute_edge_loglikelihoods(
            placement_view,
            d_ops,
            num_ops,
            d_query_pmat,
            d_pmat_mid_dist,
            d_pmat_mid_prox,
            d_prev_loglk,
            stream);
    }

    // Stage 4: iteratively optimize pendant/proximal branch lengths per op.
    const int opt_passes = std::max(std::max(1, tuning.full_opt_passes), smoothing);
    const size_t op_rate_work = num_ops_count * rate_count;
    const size_t node_rate_work = node_count * rate_count;
    for (int pass = 0; pass < opt_passes; ++pass) {
        dim3 current_deriv_grid(num_ops_grid_y);
        dim3 current_pmat_grid(grid_x(op_rate_work, pmat_block.x));
        LikelihoodDerivativePendantKernel<<<current_deriv_grid, pendant_block, shmem_bytes, stream>>>(
            placement_view,
            d_ops,
            0,
            nullptr,
            nullptr,
            nullptr,
            0.0,
            d_sumtable,
            D.d_pattern_weights_u,
            30,
            placement_view.d_new_pendant_length,
            sumtable_stride,
            placement_view.d_prev_pendant_length,
            d_active_ops,
            tuning.pendant_branch_min);
        check_launch("LikelihoodDerivativePendantKernel");

        // Rebuild query-side PMATs from the updated pendant lengths.
        BuildPendantPMATPerOpKernel<<<current_pmat_grid, pmat_block, 0, stream>>>(
            d_ops,
            nullptr,
            placement_view.d_new_pendant_length,
            D.d_Vinv,
            D.d_V,
            D.d_lambdas,
            0.0,
            d_query_pmat,
            D.states,
            D.rate_cats,
            num_ops,
            D.N,
            OPT_BRANCH_LEN_MIN,
            OPT_BRANCH_LEN_MAX,
            DEFAULT_BRANCH_LENGTH);
        check_launch("BuildPendantPMATPerOpKernel refine");
        LikelihoodDerivativeProximalKernel<<<current_deriv_grid, proximal_block, shmem_bytes, stream>>>(
            placement_view,
            d_ops,
            0,
            nullptr,
            nullptr,
            nullptr,
            0.0,
            d_sumtable,
            D.d_pattern_weights_u,
            30,
            placement_view.d_new_proximal_length,
            sumtable_stride,
            placement_view.d_prev_proximal_length,
            d_active_ops,
            tuning.split_branch_min);
        check_launch("LikelihoodDerivativeProximalKernel");

        // Rebuild midpoint PMATs from the updated proximal lengths.
        {
            dim3 pmat_grid(grid_x(node_rate_work, pmat_block.x));
            BuildNodeProximalPMATKernel<<<pmat_grid, pmat_block, 0, stream>>>(
                placement_view.d_new_proximal_length,
                D.d_Vinv,
                D.d_V,
                D.d_lambdas,
                0.0,
                d_pmat_mid_prox,
                D.states,
                D.rate_cats,
                D.N,
                D.root_id,
                OPT_BRANCH_LEN_MIN,
                OPT_BRANCH_LEN_MAX,
                DEFAULT_BRANCH_LENGTH);
            check_launch("BuildNodeProximalPMATKernel refine");
        }

        // Rebuild distal PMATs from total branch length minus proximal length.
        {
            dim3 pmat_grid(grid_x(node_rate_work, pmat_block.x));
            BuildNodeDistalPMATKernel<<<pmat_grid, pmat_block, 0, stream>>>(
                D.d_blen,
                placement_view.d_new_proximal_length,
                D.d_Vinv,
                D.d_V,
                D.d_lambdas,
                0.0,
                d_pmat_mid_dist,
                D.states,
                D.rate_cats,
                D.N,
                D.root_id,
                OPT_BRANCH_LEN_MIN,
                OPT_BRANCH_LEN_MAX,
                DEFAULT_BRANCH_LENGTH);
            check_launch("BuildNodeDistalPMATKernel refine");
        }

        // Score each placement op after the pendant/proximal updates.
        mlipper::likelihood::placement::compute_edge_loglikelihoods(
            placement_view,
            d_ops,
            num_ops,
            d_query_pmat,
            d_pmat_mid_dist,
            d_pmat_mid_prox,
            d_likelihoods,
            stream);

        dim3 keep_block(256);
        dim3 keep_grid(grid_x(num_ops_count, keep_block.x));
        KeepBestBranchLengthsKernel<<<keep_grid, keep_block, 0, stream>>>(
            d_ops,
            nullptr,
            d_likelihoods,
            d_prev_loglk,
            placement_view.d_new_pendant_length,
            placement_view.d_new_proximal_length,
            placement_view.d_prev_pendant_length,
            placement_view.d_prev_proximal_length,
            d_active_ops,
            num_ops,
            D.N);
        check_launch("KeepBestBranchLengthsKernel");
    }

    // Stage 5: collect either the fixed-size internal ranking or every
    // candidate needed for accumulated-LWR output.
    std::vector<int> final_top_indices;
    std::vector<fp_t> final_top_values;
    const int result_limit =
        tuning.export_accumulated_lwr_threshold > 0.0
            ? num_ops
            : std::max(1, tuning.export_placement_topk);
    fetch_ranked_loglikelihoods(
        d_prev_loglk,
        num_ops,
        result_limit,
        scratch,
        stream,
        final_top_indices,
        final_top_values,
        check_cuda);
    if (enable_local_child_refine &&
        d_ops &&
        num_ops > 0 &&
        !final_top_indices.empty() &&
        !final_top_values.empty() &&
        host_loglk_cache.size() == static_cast<size_t>(num_ops)) {
        const std::vector<NodeOpInfo>& host_ops =
            ensure_host_ops_loaded(d_ops, num_ops, stream, postprocess_cache);
        include_best_target_children(
            D,
            host_ops,
            host_loglk_cache,
            result_limit,
            final_top_indices,
            final_top_values);
    }
    if (final_top_indices.empty() || final_top_values.empty()) {
        throw std::runtime_error("EvaluatePlacementCandidates: no placement candidates produced");
    }

    result.top_placements = build_top_ranked_placements(
        D,
        d_ops,
        final_top_indices,
        final_top_values,
        scratch,
        stream);
    if (result.top_placements.empty()) {
        throw std::runtime_error("EvaluatePlacementCandidates: failed to materialize ranked placements");
    }

    result.target_id = result.top_placements.front().target_id;
    result.loglikelihood = result.top_placements.front().loglikelihood;
    result.distal_length = result.top_placements.front().distal_length;
    result.pendant_length = result.top_placements.front().pendant_length;

    // Stage 6: apply host-side postprocessing to the assembled ranking.
#if !defined(MLIPPER_USE_DOUBLE)
    if (d_ops &&
        final_top_indices.size() > 1 &&
        result.top_placements.size() > 1) {
        maybe_apply_double_rerank(
            placement_view,
            d_ops,
            final_top_indices,
            result,
            ensure_host_placement_eval_inputs(placement_view, stream, postprocess_cache),
            tuning);
    }
#endif
    if (enable_local_child_refine &&
        d_ops &&
        num_ops > 0 &&
        result.target_id >= 0 &&
        result.target_id < D.N) {
        const std::vector<NodeOpInfo>& host_ops =
            ensure_host_ops_loaded(d_ops, num_ops, stream, postprocess_cache);
        if (!host_ops.empty()) {
            rerank_selected_target_and_children(
                placement_view,
                result,
                host_ops,
                ensure_host_placement_eval_inputs(placement_view, stream, postprocess_cache),
                tuning);
        }
    }
    retain_accumulated_lwr(
        result.top_placements,
        tuning.export_accumulated_lwr_threshold);
    return result;
}
