#include <cuda_runtime.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <utility>
#include <stdexcept>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <limits>
#include <chrono>
#include <functional>
#include "util/mlipper_util.h"
#include "tree.hpp" 
#include "pmat.h"
#include "likelihood/root_likelihood.cuh"
#include "placement/placement.cuh"
#include "likelihood/partial_likelihood.cuh"
#include <tbb/parallel_for.h>

namespace {

using SteadyClock = std::chrono::steady_clock;

static double elapsed_ms(
    const SteadyClock::time_point& start,
    const SteadyClock::time_point& end)
{
    return std::chrono::duration<double, std::milli>(end - start).count();
}

static double time_stream_stage_ms(
    cudaStream_t stream,
    const std::function<void()>& launch_fn)
{
    cudaEvent_t start = nullptr;
    cudaEvent_t stop = nullptr;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    CUDA_CHECK(cudaEventRecord(start, stream));
    launch_fn();
    CUDA_CHECK(cudaEventRecord(stop, stream));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float stage_ms = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&stage_ms, start, stop));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
    return static_cast<double>(stage_ms);
}

// ----- Downward-op construction helpers -----

static inline int down_op_type_for_target(bool target_is_tip, bool sibling_is_tip) {
    if (target_is_tip && sibling_is_tip) return static_cast<int>(OP_DOWN_TIP_TIP);
    if (target_is_tip) return static_cast<int>(OP_DOWN_TIP_INNER);
    if (sibling_is_tip) return static_cast<int>(OP_DOWN_INNER_TIP);
    return static_cast<int>(OP_DOWN_INNER_INNER);
}

static inline NodeOpInfo make_downward_op(
    int parent_id,
    int left_id,
    int right_id,
    bool left_is_tip,
    bool right_is_tip,
    const std::vector<int>& node_to_tip,
    uint8_t dir_tag)
{
    NodeOpInfo op{};
    op.parent_id = parent_id;
    op.left_id   = left_id;
    op.right_id  = right_id;
    op.left_tip_index  = left_is_tip  ? node_to_tip[left_id]  : -1;
    op.right_tip_index = right_is_tip ? node_to_tip[right_id] : -1;
    op.clv_pool = static_cast<uint8_t>(CLV_POOL_DOWN);
    op.dir_tag  = dir_tag;

    const bool target_is_left = (dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
    const bool target_is_tip  = target_is_left ? left_is_tip : right_is_tip;
    const bool sibling_is_tip = target_is_left ? right_is_tip : left_is_tip;
    op.op_type = down_op_type_for_target(target_is_tip, sibling_is_tip);
    return op;
}

} // namespace

static void rebuild_traversals(TreeBuildResult& T) {
    T.preorder.clear();
    T.postorder.clear();
    if (T.root_id < 0 || T.root_id >= (int)T.nodes.size()) return;

    std::vector<int> stack;
    stack.push_back(T.root_id);
    while (!stack.empty()) {
        int id = stack.back();
        stack.pop_back();
        T.preorder.push_back(id);
        const TreeNode& n = T.nodes[id];
        if (!n.is_tip) {
            if (n.right >= 0) stack.push_back(n.right);
            if (n.left >= 0) stack.push_back(n.left);
        }
    }

    std::vector<std::pair<int, bool>> st;
    st.emplace_back(T.root_id, false);
    while (!st.empty()) {
        const std::pair<int, bool> cur = st.back();
        st.pop_back();
        const int id = cur.first;
        const bool visited = cur.second;
        if (id < 0) continue;
        const TreeNode& n = T.nodes[id];
        if (visited) {
            T.postorder.push_back(id);
        } else {
            st.emplace_back(id, true);
            if (!n.is_tip) {
                if (n.right >= 0) st.emplace_back(n.right, false);
                if (n.left >= 0) st.emplace_back(n.left, false);
            }
        }
    }
}

struct InsertResult {
    int internal_id = -1;
    int tip_id = -1;
    std::string tip_name;
};

// ----- Host-side PMAT staging -----

// Build transition probability matrices for each node/rate category on host.
void fill_pmats_in_host_packing(
    const TreeBuildResult&       T,
    HostPacking&                 H,
    const EigResult&             er,
    const std::vector<double>&   rate_multipliers,   // len = rate_cats (per-category rate multipliers)
    int states,
    int rate_cats,
    const int* changed_nodes,
    int num_changed_nodes)
{
    const int N = (int)T.nodes.size();
    const size_t per_node = (size_t)rate_cats * states * states;

    const size_t required = (size_t)N * per_node;
    const bool want_incremental = (changed_nodes && num_changed_nodes > 0);

    auto reset_all = [&]() {
        H.pmats.assign(required, fp_t(0));
        H.pmats_mid.assign(required, fp_t(0));
        H.pmats_mid_prox.assign(required, fp_t(0));
        H.pmats_mid_dist.assign(required, fp_t(0));
    };

    auto buffers_look_compatible = [&]() -> bool {
        if (per_node == 0) return false;
        if (H.pmats.empty()) return false;
        if (H.pmats.size() % per_node != 0) return false;
        if (!H.pmats_mid.empty() && H.pmats_mid.size() % per_node != 0) return false;
        if (!H.pmats_mid_prox.empty() && H.pmats_mid_prox.size() % per_node != 0) return false;
        if (!H.pmats_mid_dist.empty() && H.pmats_mid_dist.size() % per_node != 0) return false;
        return true;
    };

    if (!want_incremental) {
        reset_all();
    } else if (!buffers_look_compatible()) {
        // No existing buffers to update incrementally; fall back to full rebuild.
        reset_all();
    } else {
        // Preserve existing pmats for unchanged nodes; extend buffers for newly added nodes.
        H.pmats.resize(required, fp_t(0));
        H.pmats_mid.resize(required, fp_t(0));
        if (!H.pmats_mid_prox.empty()) H.pmats_mid_prox.resize(required, fp_t(0));
        if (!H.pmats_mid_dist.empty()) H.pmats_mid_dist.resize(required, fp_t(0));
    }

    auto compute_node_pmats = [&](int nid) {
        if (nid < 0 || nid >= N) return;
        const TreeNode& nd = T.nodes[nid];
        if (nd.parent < 0) return;
        fp_t* base = H.pmats.data() + (size_t)nid * per_node;
        fp_t* base_mid = H.pmats_mid.data() + (size_t)nid * per_node;
        fp_t* base_mid_prox = H.pmats_mid_prox.empty() ? nullptr : (H.pmats_mid_prox.data() + (size_t)nid * per_node);
        fp_t* base_mid_dist = H.pmats_mid_dist.empty() ? nullptr : (H.pmats_mid_dist.data() + (size_t)nid * per_node);
        const double blen  = static_cast<double>(nd.branch_length_to_parent);
        std::vector<double> pbuf((size_t)states * (size_t)states);
        std::vector<double> pbuf_mid((size_t)states * (size_t)states);

        for (int rc = 0; rc < rate_cats; ++rc) {
            double r = rate_multipliers[rc];  // rate category multiplier
            double t = blen;                  // branch length
            double p = 0;

            fp_t* P = base + (size_t)rc * states * states;
            fp_t* Pmid = base_mid + (size_t)rc * states * states;
            fp_t* Pprox = base_mid_prox ? (base_mid_prox + (size_t)rc * states * states) : nullptr;
            fp_t* Pdist = base_mid_dist ? (base_mid_dist + (size_t)rc * states * states) : nullptr;
            pmatrix_from_triple(
                er.Vinv.data(), er.V.data(), er.lambdas.data(),
                            r, t, p, pbuf.data(), states);
            // half-branch PMAT for midpoint
            pmatrix_from_triple(
                er.Vinv.data(), er.V.data(), er.lambdas.data(),
                            r, t * 0.5, p, pbuf_mid.data(), states);
            for (size_t idx = 0; idx < pbuf.size(); ++idx) {
                P[idx] = static_cast<fp_t>(pbuf[idx]);
                Pmid[idx] = static_cast<fp_t>(pbuf_mid[idx]);
            }
            if (Pprox) std::copy(Pmid, Pmid + pbuf.size(), Pprox);
            if (Pdist) std::copy(Pmid, Pmid + pbuf.size(), Pdist);

        }
    };

    if (!want_incremental || !changed_nodes || num_changed_nodes <= 0) {
        tbb::parallel_for(0, N, [&](int nid) { compute_node_pmats(nid);});
    } else {
        for (int i = 0; i < num_changed_nodes; ++i) compute_node_pmats(changed_nodes[i]);
    }
}

// ----- Per-query CLV staging -----

// Select view for a specific query's PMAT chunk.
DeviceTree make_query_view(const DeviceTree& D, int query_idx) {
    DeviceTree view = D;
    if (query_idx < 0 || query_idx >= D.placement_queries) return view;
    const size_t clv_span =
        D.sites * static_cast<size_t>(D.rate_cats) * static_cast<size_t>(D.states);
    if (D.d_query_pmat) {
        // Shared query PMAT buffer across queries; rebuilt per-query.
        view.d_query_pmat = D.d_query_pmat;
    }
    if (D.d_query_clv) {
        view.d_query_clv = D.d_query_clv + static_cast<size_t>(query_idx) * clv_span;
    }
    return view;
}

__global__ void BuildQueryClvKernel(
    DeviceTree D,
    const uint8_t* query_chars,
    int query_idx)
{
    const size_t site =
        static_cast<size_t>(blockIdx.x) * static_cast<size_t>(blockDim.x) +
        static_cast<size_t>(threadIdx.x);
    if (site >= D.sites) return;
    const size_t base_char = static_cast<size_t>(query_idx) * D.sites + site;
    const uint8_t enc = query_chars ? query_chars[base_char] : (D.states == 4 ? 15 : 4);

    const size_t per_site = static_cast<size_t>(D.rate_cats) * static_cast<size_t>(D.states);
    const size_t clv_span = D.sites * per_site;
    fp_t* out = D.d_query_clv + static_cast<size_t>(query_idx) * clv_span + site * per_site;
    if (D.states == 4) {
        // Keep query decoding on the exact same contract as reference tips:
        // query_chars stores a DNA4 bitmask, and d_tipmap is identity in the
        // 4-state case. Using the shared table avoids any drift between the
        // reference-tip and query-tip paths.
        const unsigned int mask = D.d_tipmap
            ? D.d_tipmap[static_cast<unsigned int>(enc)]
            : static_cast<unsigned int>(enc);
        for (int rc = 0; rc < D.rate_cats; ++rc) {
            fp_t* row = out + static_cast<size_t>(rc) * static_cast<size_t>(D.states);
            for (int s = 0; s < D.states; ++s) {
                row[s] = (mask & (1u << s)) ? fp_t(1) : fp_t(0);
            }
        }
    } else {
        for (int rc = 0; rc < D.rate_cats; ++rc) {
            fp_t* row = out + static_cast<size_t>(rc) * static_cast<size_t>(D.states);
            for (int s = 0; s < D.states; ++s) {
                row[s] = (enc < D.states) ? (s == enc ? fp_t(1) : fp_t(0)) : fp_t(1);
            }
        }
    }
}

void build_query_clv(
    const DeviceTree& D,
    int query_idx,
    cudaStream_t stream)
{
    if (!D.d_query_clv || !D.d_query_chars) return;
    if (query_idx < 0 || query_idx >= D.placement_queries) return;
    dim3 block(256);
    dim3 grid(static_cast<unsigned int>((D.sites + block.x - 1) / block.x));
    BuildQueryClvKernel<<<grid, block, 0, stream>>>(D, D.d_query_chars, query_idx);
    CUDA_CHECK(cudaGetLastError());
}

__global__ void SeedRootDownClvKernel(
    fp_t* d_clv_down,
    size_t per_node,
    int root_id)
{
    const size_t idx = static_cast<size_t>(blockIdx.x) * static_cast<size_t>(blockDim.x) +
        static_cast<size_t>(threadIdx.x);
    if (idx >= per_node) return;
    d_clv_down[static_cast<size_t>(root_id) * per_node + idx] = fp_t(1);
}

__global__ void CopyUnscaledUpClvToQuerySlotKernel(
    DeviceTree src,
    int src_node_id,
    DeviceTree dst,
    int dst_query_idx)
{
    const size_t site =
        static_cast<size_t>(blockIdx.x) * static_cast<size_t>(blockDim.x) +
        static_cast<size_t>(threadIdx.x);
    if (site >= src.sites) return;
    if (!src.d_clv_up || !dst.d_query_clv) return;

    const size_t per_site = static_cast<size_t>(src.rate_cats) * static_cast<size_t>(src.states);
    const size_t src_clv_span = src.sites * per_site;
    const size_t dst_clv_span = dst.sites * per_site;
    const size_t src_site_base =
        static_cast<size_t>(src_node_id) * src_clv_span + site * per_site;
    const size_t dst_site_base =
        static_cast<size_t>(dst_query_idx) * dst_clv_span + site * per_site;

    const unsigned* scaler_base = nullptr;
    if (src.d_site_scaler_up) {
        const size_t scaler_span = src.per_rate_scaling
            ? src.sites * static_cast<size_t>(src.rate_cats)
            : src.sites;
        scaler_base = src.d_site_scaler_up + static_cast<size_t>(src_node_id) * scaler_span;
    }

    for (int rc = 0; rc < src.rate_cats; ++rc) {
        unsigned shift = 0u;
        if (scaler_base) {
            shift = src.per_rate_scaling
                ? scaler_base[site * static_cast<size_t>(src.rate_cats) + static_cast<size_t>(rc)]
                : scaler_base[site];
        }
        for (int state = 0; state < src.states; ++state) {
            const size_t offset =
                static_cast<size_t>(rc) * static_cast<size_t>(src.states) + static_cast<size_t>(state);
            fp_t value = src.d_clv_up[src_site_base + offset];
            if (shift) {
                value = fp_ldexp(value, -static_cast<int>(shift));
            }
            dst.d_query_clv[dst_site_base + offset] = value;
        }
    }
}

void copy_unscaled_up_clv_to_query_slot(
    const DeviceTree& src,
    int src_node_id,
    DeviceTree& dst,
    int dst_query_idx,
    cudaStream_t stream)
{
    if (src_node_id < 0 || src_node_id >= src.N) {
        throw std::runtime_error("copy_unscaled_up_clv_to_query_slot: source node id out of range.");
    }
    if (src.sites != dst.sites || src.states != dst.states || src.rate_cats != dst.rate_cats) {
        throw std::runtime_error("copy_unscaled_up_clv_to_query_slot: source/destination dimensions mismatch.");
    }
    if (!src.d_clv_up) {
        throw std::runtime_error("copy_unscaled_up_clv_to_query_slot: missing source upward CLV buffer.");
    }
    if (!dst.d_query_clv) {
        throw std::runtime_error("copy_unscaled_up_clv_to_query_slot: missing destination query CLV buffer.");
    }
    if (dst_query_idx < 0 || dst_query_idx >= dst.query_capacity) {
        throw std::runtime_error("copy_unscaled_up_clv_to_query_slot: destination query slot out of range.");
    }
    dim3 block(256);
    dim3 grid(static_cast<unsigned>((src.sites + block.x - 1) / block.x));
    CopyUnscaledUpClvToQuerySlotKernel<<<grid, block, 0, stream>>>(
        src,
        src_node_id,
        dst,
        dst_query_idx);
    CUDA_CHECK(cudaGetLastError());
}

void copy_upward_state(
    const DeviceTree& src,
    DeviceTree& dst,
    cudaStream_t stream)
{
    if (src.N != dst.N) {
        throw std::runtime_error("copy_upward_state: node count mismatch.");
    }
    if (src.per_node_elems() != dst.per_node_elems()) {
        throw std::runtime_error("copy_upward_state: CLV shape mismatch.");
    }
    if (src.scaler_elems() != dst.scaler_elems()) {
        throw std::runtime_error("copy_upward_state: scaler shape mismatch.");
    }

    const size_t clv_elems = static_cast<size_t>(src.N) * src.per_node_elems();
    if (clv_elems > 0 && src.d_clv_up && dst.d_clv_up) {
        CUDA_CHECK(cudaMemcpyAsync(
            dst.d_clv_up,
            src.d_clv_up,
            sizeof(fp_t) * clv_elems,
            cudaMemcpyDeviceToDevice,
            stream));
    }

    const size_t scaler_elems = static_cast<size_t>(src.N) * src.scaler_elems();
    if (scaler_elems > 0 && src.d_site_scaler_up && dst.d_site_scaler_up) {
        CUDA_CHECK(cudaMemcpyAsync(
            dst.d_site_scaler_up,
            src.d_site_scaler_up,
            sizeof(unsigned) * scaler_elems,
            cudaMemcpyDeviceToDevice,
            stream));
    }
}

// Pack topology and tip encodings from tree/MSA into HostPacking.
HostPacking pack_host_arrays_from_tree_and_msa(
        const TreeBuildResult& T,
        const std::vector<std::string>& msa_tip_names,  // len = tips
        const std::vector<std::string>& msa_rows,       // len = tips, each row length = sites
        size_t sites,
        int states)
{
    if (msa_rows.size() != msa_tip_names.size())
        throw std::runtime_error("MSA rows/names size mismatch.");
    if (msa_rows.empty()) throw std::runtime_error("Empty MSA.");
    if (msa_rows[0].size() != sites)
        throw std::runtime_error("Sites mismatch.");

    const int N = (int)T.nodes.size();

    // Topology
    HostPacking H;
    H.postorder = T.postorder;
    H.preorder  = T.preorder;
    H.parent.resize(N, -1);
    H.left.resize(N, -1);
    H.right.resize(N, -1);
    H.is_tip.resize(N, 0);
    H.blen.resize(N, fp_t(0));

    for (int i = 0; i < N; ++i) {
        const auto& nd = T.nodes[i];
        H.parent[i] = nd.parent;
        H.left[i]   = nd.left;
        H.right[i]  = nd.right;
        H.is_tip[i] = nd.is_tip ? 1 : 0;
        H.blen[i]   = nd.branch_length_to_parent;
    }
    
    // Build tip name -> MSA row lookup
    std::unordered_map<std::string,int> name_to_row;
    name_to_row.reserve(msa_tip_names.size()*2);
    for (int row_idx = 0; row_idx < (int)msa_tip_names.size(); ++row_idx) {
        const auto [_, inserted] = name_to_row.emplace(msa_tip_names[row_idx], row_idx);
        if (!inserted) {
            throw std::runtime_error("Duplicate MSA tip name: " + msa_tip_names[row_idx]);
        }
    }

    // Collect tip order in postorder and the corresponding node ids
    std::vector<int> tip_node_ids_host;
    tip_node_ids_host.reserve(msa_tip_names.size());
    for (int id : T.postorder) {
        if (T.nodes[id].is_tip) tip_node_ids_host.push_back(id);
    }
    const int tips = (int)tip_node_ids_host.size();

    if (tips != (int)msa_tip_names.size())
        throw std::runtime_error("Tip count in tree != MSA names.");

    H.tip_node_ids = tip_node_ids_host;

    // Encode MSA rows into tipchars following the fixed tipIndex order (0..tips-1)
    H.tipchars.resize((size_t)tips * sites);

    if (states == 4 || states == 5) {
        for (int t = 0; t < tips; ++t) {
            const int node_id = tip_node_ids_host[t];
            const auto& name  = T.nodes[node_id].name;
            auto it = name_to_row.find(name);
            if (it == name_to_row.end())
                throw std::runtime_error("Tip not found in MSA: " + name);
            const std::string& row = msa_rows[it->second];
            for (size_t s = 0; s < sites; ++s)
                H.tipchars[(size_t)t * sites + s] =
                    (states == 4) ? encode_state_DNA4_mask(row[s]) : encode_state_DNA5(row[s]);
        }
    } else {
        throw std::runtime_error("Unsupported states for tip encoding; expected 4 or 5.");
    }

    return H;
}

namespace {

static inline std::vector<fp_t> cast_to_fp(const std::vector<double>& src) {
    std::vector<fp_t> out(src.size());
    for (size_t i = 0; i < src.size(); ++i) {
        out[i] = static_cast<fp_t>(src[i]);
    }
    return out;
}

} // namespace

// ===== Copy host packing to GPU and build DeviceTree =====
// Upload host packing and model parameters to GPU, constructing DeviceTree.
DeviceTree upload_to_gpu(
    const TreeBuildResult& T,
    const HostPacking& H,
    const EigResult& er,
    const std::vector<double>& rate_weights,
    const std::vector<double>& rate_multipliers,
    const std::vector<double>& pi,
    size_t sites, int states, int rate_cats, bool per_rate_scaling,
    const PlacementQueryBatch* queries,
    bool commit_to_tree)
{
    DeviceTree D;
    const int node_count = static_cast<int>(T.nodes.size());
    const int tip_count = static_cast<int>(H.tip_node_ids.size());
    const int query_count = queries ? static_cast<int>(queries->size()) : 0;
    const size_t site_count = sites;
    const size_t state_count = static_cast<size_t>(states);
    const size_t rate_count = static_cast<size_t>(rate_cats);
    const size_t matrix_elems = state_count * state_count;

    D.N = node_count;
    D.tips = tip_count;
    D.inners = D.N - D.tips;
    D.placement_queries = query_count;
    D.root_id = T.root_id;

    // Pre-allocate device memory for tree growth during commit-to-tree:
    // Each query placement adds exactly 1 internal node + 1 tip node = 2 nodes
    int reserve_inserts = commit_to_tree ? D.placement_queries : 0;
    
    D.capacity_N = D.N + 2 * reserve_inserts;      // +2 nodes per query (1 internal + 1 tip)
    D.capacity_tips = D.tips + reserve_inserts;    // +1 tip per query
    D.query_capacity = D.placement_queries;
    D.sites = sites;
    D.states = states;
    D.rate_cats = rate_cats;
    D.log2_stride = ceil_log2_u32((unsigned int)(D.states + 1));
    D.per_rate_scaling = per_rate_scaling;

    if (static_cast<int>(rate_weights.size()) != rate_cats) {
        throw std::runtime_error("rate_weights size mismatch.");
    }
    if (static_cast<int>(pi.size()) != states) {
        throw std::runtime_error("pi size mismatch.");
    }

    // Host-side staging for model uploads.
    std::vector<fp_t> lambdas_scaled(rate_count * state_count, fp_t(0));
    for (int rate_idx = 0; rate_idx < D.rate_cats; ++rate_idx) {
        const double rate_multiplier = rate_multipliers[rate_idx];
        for (int state_idx = 0; state_idx < D.states; ++state_idx) {
            lambdas_scaled[static_cast<size_t>(rate_idx) * state_count + static_cast<size_t>(state_idx)] =
                static_cast<fp_t>(er.lambdas[state_idx] * rate_multiplier);
        }
    }

    const auto V_fp = cast_to_fp(er.V);
    const auto Vinv_fp = cast_to_fp(er.Vinv);
    const auto rate_weights_fp = cast_to_fp(rate_weights);
    const auto pi_fp = cast_to_fp(pi);

    const size_t model_matrix_bytes = sizeof(fp_t) * matrix_elems;
    const size_t lambdas_bytes = sizeof(fp_t) * lambdas_scaled.size();
    const size_t rate_weights_bytes = sizeof(fp_t) * rate_count;
    const size_t frequencies_bytes = sizeof(fp_t) * state_count;

    // --- model parameters ---
    CUDA_CHECK(cudaMalloc(&D.d_lambdas, lambdas_bytes));
    CUDA_CHECK(cudaMalloc(&D.d_V, model_matrix_bytes));
    CUDA_CHECK(cudaMalloc(&D.d_Vinv, model_matrix_bytes));
    CUDA_CHECK(cudaMalloc(&D.d_rate_weights, rate_weights_bytes));
    CUDA_CHECK(cudaMalloc(&D.d_frequencies, frequencies_bytes));

    CUDA_CHECK(cudaMemcpy(D.d_lambdas, lambdas_scaled.data(), lambdas_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_V, V_fp.data(), model_matrix_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_Vinv, Vinv_fp.data(), model_matrix_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_rate_weights, rate_weights_fp.data(), rate_weights_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(D.d_frequencies, pi_fp.data(), frequencies_bytes, cudaMemcpyHostToDevice));

    // --- branch lengths and tip characters ---
    const size_t capacity_node_bytes = sizeof(fp_t) * static_cast<size_t>(D.capacity_N);
    const size_t live_node_bytes = sizeof(fp_t) * static_cast<size_t>(D.N);
    const size_t tipchar_capacity_bytes = sizeof(uint8_t) * static_cast<size_t>(D.capacity_tips) * site_count;
    const size_t tipchar_live_bytes = sizeof(uint8_t) * static_cast<size_t>(D.tips) * site_count;
    const size_t tip_node_ids_capacity_bytes = sizeof(int) * static_cast<size_t>(D.capacity_tips);
    const size_t tip_node_ids_live_bytes = sizeof(int) * static_cast<size_t>(D.tips);

    CUDA_CHECK(cudaMalloc(&D.d_blen, capacity_node_bytes));
    CUDA_CHECK(cudaMalloc(&D.d_new_pendant_length, capacity_node_bytes));
    CUDA_CHECK(cudaMalloc(&D.d_new_proximal_length, capacity_node_bytes));
    CUDA_CHECK(cudaMalloc(&D.d_prev_pendant_length, capacity_node_bytes));
    CUDA_CHECK(cudaMalloc(&D.d_prev_proximal_length, capacity_node_bytes));
    CUDA_CHECK(cudaMemset(D.d_blen, 0, capacity_node_bytes));
    CUDA_CHECK(cudaMemcpy(D.d_blen, H.blen.data(), live_node_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(D.d_new_pendant_length, 0, capacity_node_bytes));
    CUDA_CHECK(cudaMemset(D.d_new_proximal_length, 0, capacity_node_bytes));
    CUDA_CHECK(cudaMemset(D.d_prev_pendant_length, 0, capacity_node_bytes));
    CUDA_CHECK(cudaMemset(D.d_prev_proximal_length, 0, capacity_node_bytes));

    CUDA_CHECK(cudaMalloc(&D.d_tipchars, tipchar_capacity_bytes));
    CUDA_CHECK(cudaMemcpy(D.d_tipchars, H.tipchars.data(), tipchar_live_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(&D.d_tip_node_ids, tip_node_ids_capacity_bytes));
    CUDA_CHECK(cudaMemcpy(D.d_tip_node_ids, H.tip_node_ids.data(), tip_node_ids_live_bytes, cudaMemcpyHostToDevice));

    // --- CLV and derivative workspaces ---
    const size_t per_node = site_count * rate_count * state_count;
    const size_t clv_capacity_elems = static_cast<size_t>(D.capacity_N) * per_node;
    const size_t clv_total = clv_capacity_elems * 2; // up + down
    const size_t clv_total_bytes = sizeof(fp_t) * clv_total;
    const size_t clv_capacity_bytes = sizeof(fp_t) * clv_capacity_elems;

    CUDA_CHECK(cudaMalloc(&D.d_clv_up, clv_total_bytes));
    D.d_clv_down = D.d_clv_up + clv_capacity_elems;
    D.clv_down_offset_elems = clv_capacity_elems;

    CUDA_CHECK(cudaMalloc(&D.d_clv_mid, clv_capacity_bytes));
    CUDA_CHECK(cudaMalloc(&D.d_clv_mid_base, clv_capacity_bytes));
    CUDA_CHECK(cudaMalloc(&D.d_root_loglik_total, sizeof(double)));
    CUDA_CHECK(cudaMemset(D.d_root_loglik_total, 0, sizeof(double)));

    const size_t sumtable_stride = site_count * rate_count * state_count;
    const size_t max_ops = static_cast<size_t>(D.capacity_N) * 2;
    D.sumtable_capacity_ops = max_ops;
    D.likelihood_capacity_ops = max_ops;
    if (sumtable_stride > 0 && max_ops > 0) {
        CUDA_CHECK(cudaMalloc(&D.d_sumtable, sizeof(fp_t) * sumtable_stride * max_ops));
    }
    if (max_ops > 0) {
        CUDA_CHECK(cudaMalloc(&D.d_likelihoods, sizeof(fp_t) * max_ops));
    }

    // Optionally zero out the CLV pool (or let kernels overwrite)
    CUDA_CHECK(cudaMemset(D.d_clv_up, 0, clv_total_bytes));
    CUDA_CHECK(cudaMemset(D.d_clv_mid, 0, clv_capacity_bytes));
    CUDA_CHECK(cudaMemset(D.d_clv_mid_base, 0, clv_capacity_bytes));

    // Initialize only the root slice of down pool to exact 1.0; others remain 0 until overwritten.
    if (per_node > 0) {
        std::vector<fp_t> ones(per_node, fp_t(1));
        CUDA_CHECK(cudaMemcpy(D.d_clv_down + (size_t)T.root_id * per_node,
                              ones.data(),
                              per_node * sizeof(fp_t),
                              cudaMemcpyHostToDevice));
    }

    // --- PMAT buffers ---
    const size_t pmat_elems_cur = static_cast<size_t>(D.N) * rate_count * matrix_elems;
    const size_t pmat_elems_cap = static_cast<size_t>(D.capacity_N) * rate_count * matrix_elems;
    const size_t pmat_bytes_cur = sizeof(fp_t) * pmat_elems_cur;
    const size_t pmat_bytes_cap = sizeof(fp_t) * pmat_elems_cap;
    CUDA_CHECK(cudaMalloc(&D.d_pmat, pmat_bytes_cap));
    CUDA_CHECK(cudaMemset(D.d_pmat, 0, pmat_bytes_cap));
    CUDA_CHECK(cudaMemcpy(D.d_pmat, H.pmats.data(), pmat_bytes_cur, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(&D.d_pmat_mid, pmat_bytes_cap));
    CUDA_CHECK(cudaMemset(D.d_pmat_mid, 0, pmat_bytes_cap));
    CUDA_CHECK(cudaMemcpy(D.d_pmat_mid, H.pmats_mid.data(), pmat_bytes_cur, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(&D.d_pmat_mid_prox, pmat_bytes_cap));
    if (!H.pmats_mid_prox.empty() && H.pmats_mid_prox.size() != pmat_elems_cur) {
        throw std::runtime_error("pmats_mid_prox size mismatch.");
    }
    const fp_t* pmat_mid_prox_src = H.pmats_mid_prox.empty() ? H.pmats_mid.data() : H.pmats_mid_prox.data();
    CUDA_CHECK(cudaMemset(D.d_pmat_mid_prox, 0, pmat_bytes_cap));
    CUDA_CHECK(cudaMemcpy(D.d_pmat_mid_prox, pmat_mid_prox_src, pmat_bytes_cur, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(&D.d_pmat_mid_dist, pmat_bytes_cap));
    if (!H.pmats_mid_dist.empty() && H.pmats_mid_dist.size() != pmat_elems_cur) {
        throw std::runtime_error("pmats_mid_dist size mismatch.");
    }
    const fp_t* pmat_mid_dist_src = H.pmats_mid_dist.empty() ? H.pmats_mid.data() : H.pmats_mid_dist.data();
    CUDA_CHECK(cudaMemset(D.d_pmat_mid_dist, 0, pmat_bytes_cap));
    CUDA_CHECK(cudaMemcpy(D.d_pmat_mid_dist, pmat_mid_dist_src, pmat_bytes_cur, cudaMemcpyHostToDevice));

    // Allocate query PMAT buffer sized for up to ~2*N placement ops (edges).
    {
        const size_t per_query = rate_count * matrix_elems;
        const size_t query_slots = static_cast<size_t>(D.capacity_N) * 2;
        const size_t query_pmat_bytes = sizeof(fp_t) * per_query * query_slots;
        CUDA_CHECK(cudaMalloc(&D.d_query_pmat, query_pmat_bytes));
        CUDA_CHECK(cudaMemset(D.d_query_pmat, 0, query_pmat_bytes));
    }

    // --- optional pattern weights ---
    if (!H.pattern_weights.empty()) {
        if (H.pattern_weights.size() != sites) {
            throw std::runtime_error("pattern_weights size mismatch.");
        }
        CUDA_CHECK(cudaMalloc(&D.d_pattern_weights_u, sizeof(unsigned) * sites));
        CUDA_CHECK(cudaMemcpy(
            D.d_pattern_weights_u,
            H.pattern_weights.data(),
            sizeof(unsigned) * sites,
            cudaMemcpyHostToDevice));
    }

    // --- optional query buffers ---
    if (queries && !queries->empty()) {
        const size_t qcount = queries->size();
        const size_t qcap = (size_t)D.query_capacity;
        const size_t chars_bytes = sizeof(uint8_t) * qcap * site_count;
        CUDA_CHECK(cudaMalloc(&D.d_query_chars, chars_bytes));
        CUDA_CHECK(cudaMemset(D.d_query_chars, 0, chars_bytes));
        CUDA_CHECK(cudaMemcpy(D.d_query_chars, queries->query_chars.data(), sizeof(uint8_t) * qcount * site_count, cudaMemcpyHostToDevice));

        const size_t query_clv_elems = qcap * site_count * rate_count * state_count;
        const size_t query_clv_bytes = sizeof(fp_t) * query_clv_elems;
        CUDA_CHECK(cudaMalloc(&D.d_query_clv, query_clv_bytes));
        CUDA_CHECK(cudaMemset(D.d_query_clv, 0, query_clv_bytes));
    }

    // --- scaler buffers ---
    const size_t scaler_span = per_rate_scaling ? (site_count * rate_count) : site_count;
    const size_t scaler_pool = static_cast<size_t>(D.capacity_N) * scaler_span;
    CUDA_CHECK(cudaMalloc(&D.d_site_scaler_storage, sizeof(unsigned) * scaler_pool * 4));
    CUDA_CHECK(cudaMemset(D.d_site_scaler_storage, 0, sizeof(unsigned) * scaler_pool * 4));
    D.d_site_scaler_up = D.d_site_scaler_storage;
    D.d_site_scaler_down = D.d_site_scaler_storage + scaler_pool;
    D.d_site_scaler_mid = D.d_site_scaler_storage + scaler_pool * 2;
    D.d_site_scaler_mid_base = D.d_site_scaler_storage + scaler_pool * 3;
    D.d_site_scaler = nullptr;

    // --- tip bitmask lookup ---
    const unsigned int tipmap_size = (D.states == 4) ? 16u : (unsigned int)D.states + 1u;
    std::vector<unsigned int> tipmap(tipmap_size);
    for (unsigned int j = 0; j < tipmap_size; ++j) {
        if (D.states == 4) {
            tipmap[j] = j;
        } else if (j == static_cast<unsigned int>(D.states)) {
            tipmap[j] = 15;
        } else {
            tipmap[j] = 1u << j;
        }
    }
    CUDA_CHECK(cudaMalloc(&D.d_tipmap, sizeof(unsigned) * tipmap_size));
    CUDA_CHECK(cudaMemcpy(D.d_tipmap, tipmap.data(), tipmap_size * sizeof(unsigned int), cudaMemcpyHostToDevice));

    return D;
}

namespace {

enum ReloadDebugSkipFlags : unsigned {
    RELOAD_SKIP_COPY_BLEN = 1u << 0,
    RELOAD_SKIP_ZERO_LENGTH_SCRATCH = 1u << 1,
    RELOAD_SKIP_COPY_TIPCHARS = 1u << 2,
    RELOAD_SKIP_ZERO_UP_DOWN_CLV = 1u << 3,
    RELOAD_SKIP_ZERO_MID_BUFFERS = 1u << 4,
    RELOAD_SKIP_ZERO_SCALERS = 1u << 5,
    RELOAD_SKIP_SEED_ROOT_DOWN = 1u << 6,
    RELOAD_SKIP_COPY_PMATS = 1u << 7,
    RELOAD_SKIP_COPY_PMAT = 1u << 8,
    RELOAD_SKIP_COPY_PMAT_MID = 1u << 9,
    RELOAD_SKIP_COPY_PMAT_MID_PROX = 1u << 10,
    RELOAD_SKIP_COPY_PMAT_MID_DIST = 1u << 11,
    RELOAD_SKIP_COPY_PATTERN_WEIGHTS = 1u << 12,
};

void reload_device_tree_live_data_impl(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const PlacementQueryBatch* queries,
    unsigned debug_skip_flags,
    cudaStream_t stream,
    DeviceTreeReloadTimingStats* timing)
{
    (void)timing;
    const int node_count = static_cast<int>(T.nodes.size());
    const int tip_count = static_cast<int>(H.tip_node_ids.size());
    const int query_count = queries ? static_cast<int>(queries->size()) : D.placement_queries;
    if (node_count <= 0) {
        throw std::runtime_error("reload_device_tree_live_data: empty tree.");
    }
    if (node_count > D.capacity_N) {
        throw std::runtime_error(
            "reload_device_tree_live_data: node count exceeds device capacity (" +
            std::to_string(node_count) + " > " + std::to_string(D.capacity_N) + ")");
    }
    if (tip_count > D.capacity_tips) {
        throw std::runtime_error(
            "reload_device_tree_live_data: tip count exceeds device capacity (" +
            std::to_string(tip_count) + " > " + std::to_string(D.capacity_tips) + ")");
    }
    if (query_count > D.query_capacity) {
        throw std::runtime_error(
            "reload_device_tree_live_data: query count exceeds device capacity (" +
            std::to_string(query_count) + " > " + std::to_string(D.query_capacity) + ")");
    }
    if (H.blen.size() != static_cast<size_t>(node_count)) {
        throw std::runtime_error("reload_device_tree_live_data: branch length host size mismatch.");
    }
    if (!H.pattern_weights.empty() && H.pattern_weights.size() != D.sites) {
        throw std::runtime_error("reload_device_tree_live_data: pattern_weights size mismatch.");
    }

    D.N = node_count;
    D.tips = tip_count;
    D.inners = D.N - D.tips;
    D.root_id = T.root_id;
    D.placement_queries = query_count;

    const size_t per_node = D.per_node_elems();
    const size_t matrix_per_node = D.pmat_per_node_elems();
    const size_t scaler_pool_bytes = sizeof(unsigned) * D.scaler_storage_elems();
    const size_t capacity_node_bytes = sizeof(fp_t) * static_cast<size_t>(D.capacity_N);
    const size_t live_node_bytes = sizeof(fp_t) * static_cast<size_t>(D.N);
    const size_t live_tipchar_bytes = sizeof(uint8_t) * static_cast<size_t>(D.tips) * D.sites;
    const size_t live_tip_node_ids_bytes = sizeof(int) * static_cast<size_t>(D.tips);
    const size_t clv_capacity_bytes = sizeof(fp_t) * static_cast<size_t>(D.capacity_N) * per_node;
    const size_t pmat_live_bytes = sizeof(fp_t) * static_cast<size_t>(D.N) * matrix_per_node;

    if ((debug_skip_flags & RELOAD_SKIP_COPY_BLEN) == 0u) {
        CUDA_CHECK(cudaMemcpyAsync(
            D.d_blen,
            H.blen.data(),
            live_node_bytes,
            cudaMemcpyHostToDevice,
            stream));
    }

    if ((debug_skip_flags & RELOAD_SKIP_ZERO_LENGTH_SCRATCH) == 0u) {
        CUDA_CHECK(cudaMemsetAsync(D.d_new_pendant_length, 0, capacity_node_bytes, stream));
        CUDA_CHECK(cudaMemsetAsync(D.d_new_proximal_length, 0, capacity_node_bytes, stream));
        CUDA_CHECK(cudaMemsetAsync(D.d_prev_pendant_length, 0, capacity_node_bytes, stream));
        CUDA_CHECK(cudaMemsetAsync(D.d_prev_proximal_length, 0, capacity_node_bytes, stream));
    }

    if ((debug_skip_flags & RELOAD_SKIP_COPY_TIPCHARS) == 0u &&
        D.tips > 0 && D.d_tipchars) {
        CUDA_CHECK(cudaMemcpyAsync(
            D.d_tipchars,
            H.tipchars.data(),
            live_tipchar_bytes,
            cudaMemcpyHostToDevice,
            stream));
        if (D.d_tip_node_ids) {
            CUDA_CHECK(cudaMemcpyAsync(
                D.d_tip_node_ids,
                H.tip_node_ids.data(),
                live_tip_node_ids_bytes,
                cudaMemcpyHostToDevice,
                stream));
        }
    }

    if ((debug_skip_flags & RELOAD_SKIP_ZERO_UP_DOWN_CLV) == 0u) {
        CUDA_CHECK(cudaMemsetAsync(D.d_clv_up, 0, clv_capacity_bytes, stream));
        CUDA_CHECK(cudaMemsetAsync(D.d_clv_down, 0, clv_capacity_bytes, stream));
    }
    if ((debug_skip_flags & RELOAD_SKIP_ZERO_MID_BUFFERS) == 0u) {
        CUDA_CHECK(cudaMemsetAsync(D.d_clv_mid, 0, clv_capacity_bytes, stream));
        CUDA_CHECK(cudaMemsetAsync(D.d_clv_mid_base, 0, clv_capacity_bytes, stream));
    }
    if ((debug_skip_flags & RELOAD_SKIP_ZERO_SCALERS) == 0u) {
        CUDA_CHECK(cudaMemsetAsync(D.d_site_scaler_storage, 0, scaler_pool_bytes, stream));
    }

    if ((debug_skip_flags & RELOAD_SKIP_SEED_ROOT_DOWN) == 0u && per_node > 0) {
        dim3 block(256);
        dim3 grid(static_cast<unsigned>((per_node + block.x - 1) / block.x));
        SeedRootDownClvKernel<<<grid, block, 0, stream>>>(
            D.d_clv_down,
            per_node,
            T.root_id);
        CUDA_CHECK(cudaGetLastError());
    }

    if (H.pmats.size() != static_cast<size_t>(D.N) * matrix_per_node ||
        H.pmats_mid.size() != static_cast<size_t>(D.N) * matrix_per_node) {
        throw std::runtime_error("reload_device_tree_live_data: PMAT host size mismatch.");
    }
    if ((debug_skip_flags & RELOAD_SKIP_COPY_PMATS) == 0u) {
        const fp_t* pmat_mid_prox_src =
            H.pmats_mid_prox.empty() ? H.pmats_mid.data() : H.pmats_mid_prox.data();
        const fp_t* pmat_mid_dist_src =
            H.pmats_mid_dist.empty() ? H.pmats_mid.data() : H.pmats_mid_dist.data();
        if ((debug_skip_flags & RELOAD_SKIP_COPY_PMAT) == 0u) {
            CUDA_CHECK(cudaMemcpyAsync(
                D.d_pmat, H.pmats.data(), pmat_live_bytes, cudaMemcpyHostToDevice, stream));
        }
        if ((debug_skip_flags & RELOAD_SKIP_COPY_PMAT_MID) == 0u) {
            CUDA_CHECK(cudaMemcpyAsync(
                D.d_pmat_mid,
                H.pmats_mid.data(),
                pmat_live_bytes,
                cudaMemcpyHostToDevice,
                stream));
        }
        if ((debug_skip_flags & RELOAD_SKIP_COPY_PMAT_MID_PROX) == 0u) {
            CUDA_CHECK(cudaMemcpyAsync(
                D.d_pmat_mid_prox,
                pmat_mid_prox_src,
                pmat_live_bytes,
                cudaMemcpyHostToDevice,
                stream));
        }
        if ((debug_skip_flags & RELOAD_SKIP_COPY_PMAT_MID_DIST) == 0u) {
            CUDA_CHECK(cudaMemcpyAsync(
                D.d_pmat_mid_dist,
                pmat_mid_dist_src,
                pmat_live_bytes,
                cudaMemcpyHostToDevice,
                stream));
        }
    }

    if (queries && query_count > 0) {
        if (!queries->query_chars.empty()) {
            const size_t query_chars_bytes = sizeof(uint8_t) * static_cast<size_t>(query_count) * D.sites;
            CUDA_CHECK(cudaMemcpyAsync(
                D.d_query_chars,
                queries->query_chars.data(),
                query_chars_bytes,
                cudaMemcpyHostToDevice,
                stream));
        }
        if (!queries->query_pmats.empty()) {
            const size_t query_pmat_bytes = sizeof(fp_t) * queries->query_pmats.size();
            CUDA_CHECK(cudaMemcpyAsync(
                D.d_query_pmat,
                queries->query_pmats.data(),
                query_pmat_bytes,
                cudaMemcpyHostToDevice,
                stream));
        }
    }

    if ((debug_skip_flags & RELOAD_SKIP_COPY_PATTERN_WEIGHTS) == 0u &&
        D.d_pattern_weights_u && !H.pattern_weights.empty()) {
        CUDA_CHECK(cudaMemcpyAsync(
            D.d_pattern_weights_u,
            H.pattern_weights.data(),
            sizeof(unsigned) * D.sites,
            cudaMemcpyHostToDevice,
            stream));
    }
}

} // namespace

void reload_device_tree_live_data(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const PlacementQueryBatch* queries,
    cudaStream_t stream,
    DeviceTreeReloadTimingStats* timing)
{
    reload_device_tree_live_data_impl(D, T, H, queries, 0u, stream, timing);
}

void reload_device_tree_live_data_local_spr(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const HostPacking& base_H,
    int current_main_pmat_node,
    int& previous_main_pmat_node,
    const PlacementQueryBatch* queries,
    cudaStream_t stream,
    DeviceTreeReloadTimingStats* timing)
{
    reload_device_tree_live_data_impl(
        D,
        T,
        H,
        queries,
        RELOAD_SKIP_COPY_PMAT |
            RELOAD_SKIP_COPY_PMAT_MID |
            RELOAD_SKIP_COPY_PMAT_MID_PROX |
            RELOAD_SKIP_COPY_PMAT_MID_DIST,
        stream,
        timing);

    const size_t per_node = D.pmat_per_node_elems();
    const size_t required = static_cast<size_t>(D.N) * per_node;
    if (!D.d_pmat || per_node == 0) {
        previous_main_pmat_node = current_main_pmat_node;
        return;
    }
    if (H.pmats.size() != required || base_H.pmats.size() != required) {
        throw std::runtime_error(
            "reload_device_tree_live_data_local_spr: PMAT host size mismatch.");
    }

    auto copy_main_pmat_slice = [&](const std::vector<fp_t>& src, int node_id) {
        if (node_id < 0 || node_id >= D.N) return;
        const size_t offset = static_cast<size_t>(node_id) * per_node;
        CUDA_CHECK(cudaMemcpyAsync(
            D.d_pmat + offset,
            src.data() + offset,
            sizeof(fp_t) * per_node,
            cudaMemcpyHostToDevice,
            stream));
    };

    if (previous_main_pmat_node >= 0 && previous_main_pmat_node != current_main_pmat_node) {
        copy_main_pmat_slice(base_H.pmats, previous_main_pmat_node);
    }
    copy_main_pmat_slice(H.pmats, current_main_pmat_node);
    previous_main_pmat_node = current_main_pmat_node;
}

// ===== Release GPU resources =====
// Release all device buffers held by DeviceTree.
void free_device_tree(DeviceTree& D) {
    auto F = [](void* p){ if(p) cudaFree(p); };
    F(D.d_blen);
    F(D.d_new_pendant_length); F(D.d_new_proximal_length);
    F(D.d_prev_pendant_length); F(D.d_prev_proximal_length);
    F(D.d_tipchars);
    F(D.d_tip_node_ids);
    // d_clv_down is an offset into d_clv_up; free only the base allocation.
    F(D.d_clv_up);
    F(D.d_clv_mid);
    F(D.d_clv_mid_base);
    F(D.d_site_scaler_storage);
    F(D.d_lambdas); F(D.d_V); F(D.d_Vinv); F(D.d_rate_weights); F(D.d_frequencies);
    F(D.d_pmat); F(D.d_pmat_mid); F(D.d_pmat_mid_prox); F(D.d_pmat_mid_dist);
    F(D.d_root_loglik_total);
    F(D.d_sumtable);
    F(D.d_likelihoods);
    F(D.d_pattern_weights_u);
    F(D.d_query_clv); F(D.d_query_chars);
    F(D.d_query_pmat);
    F(D.d_tipmap);
    D = DeviceTree{};
}


static void rebuild_node_to_tip_map(
    const TreeBuildResult& tree,
    const HostPacking& host,
    std::vector<int>& node_to_tip)
{
    const int node_count = static_cast<int>(tree.nodes.size());
    node_to_tip.assign(static_cast<size_t>(node_count), -1);
    for (int tip_idx = 0; tip_idx < static_cast<int>(host.tip_node_ids.size()); ++tip_idx) {
        const int node_id = host.tip_node_ids[tip_idx];
        if (node_id >= 0 && node_id < node_count) {
            node_to_tip[node_id] = tip_idx;
        }
    }
}

static bool make_upward_op_for_node(
    const TreeBuildResult& tree,
    const std::vector<int>& node_to_tip,
    int node_id,
    NodeOpInfo& op)
{
    const int node_count = static_cast<int>(tree.nodes.size());
    if (node_id < 0 || node_id >= node_count) return false;
    const TreeNode& node = tree.nodes[node_id];
    if (node.is_tip) return false;

    const int left_id = node.left;
    const int right_id = node.right;
    if (left_id < 0 || right_id < 0) return false;

    const bool left_is_tip = tree.nodes[left_id].is_tip;
    const bool right_is_tip = tree.nodes[right_id].is_tip;
    const int left_tip_idx = left_is_tip ? node_to_tip[left_id] : -1;
    const int right_tip_idx = right_is_tip ? node_to_tip[right_id] : -1;

    NodeOpType op_type = OP_INNER_INNER;
    if (left_is_tip && right_is_tip) {
        op_type = OP_TIP_TIP;
    } else if (left_is_tip || right_is_tip) {
        op_type = OP_TIP_INNER;
    }

    op = NodeOpInfo{};
    op.parent_id = node_id;
    op.left_id = left_id;
    op.right_id = right_id;
    op.left_tip_index = left_tip_idx;
    op.right_tip_index = right_tip_idx;
    op.op_type = static_cast<int>(op_type);
    op.clv_pool = static_cast<uint8_t>(CLV_POOL_UP);
    op.dir_tag = static_cast<uint8_t>(CLV_DIR_UP);
    return true;
}

static void build_upward_ops_host(
    const TreeBuildResult& tree,
    const std::vector<int>& node_to_tip,
    std::vector<NodeOpInfo>& upward_ops_host)
{
    const int node_count = static_cast<int>(tree.nodes.size());
    upward_ops_host.clear();
    upward_ops_host.reserve(static_cast<size_t>(node_count));

    for (int node_id : tree.postorder) {
        NodeOpInfo op{};
        if (make_upward_op_for_node(tree, node_to_tip, node_id, op)) {
            upward_ops_host.push_back(op);
        }
    }
}

static void build_upward_ops_host_for_path(
    const TreeBuildResult& tree,
    const std::vector<int>& node_to_tip,
    int start_node_id,
    std::vector<NodeOpInfo>& upward_ops_host)
{
    upward_ops_host.clear();
    if (start_node_id < 0 || start_node_id >= static_cast<int>(tree.nodes.size())) {
        return;
    }

    int depth = 0;
    for (int node_id = start_node_id; node_id >= 0; node_id = tree.nodes[node_id].parent) {
        ++depth;
    }
    upward_ops_host.reserve(static_cast<size_t>(depth));

    for (int node_id = start_node_id; node_id >= 0; node_id = tree.nodes[node_id].parent) {
        NodeOpInfo op{};
        if (make_upward_op_for_node(tree, node_to_tip, node_id, op)) {
            upward_ops_host.push_back(op);
        }
    }
}

static void build_downward_ops_host_levelized(
    const TreeBuildResult& tree,
    const std::vector<int>& node_to_tip,
    std::vector<NodeOpInfo>& downward_ops_host,
    std::vector<int>& level_offsets)
{
    // Downward dependencies only flow from a parent to the next child depth, so
    // children at the same depth can be launched together as one wavefront.
    const int node_count = static_cast<int>(tree.nodes.size());
    downward_ops_host.clear();
    level_offsets.clear();
    if (node_count <= 0 || tree.root_id < 0 || tree.root_id >= node_count) {
        return;
    }

    std::vector<int> node_depth(static_cast<size_t>(node_count), -1);
    std::vector<int> level_counts(1, 0);
    node_depth[static_cast<size_t>(tree.root_id)] = 0;

    for (int parent_id : tree.preorder) {
        if (parent_id < 0 || parent_id >= node_count) continue;
        const int parent_depth = node_depth[static_cast<size_t>(parent_id)];
        if (parent_depth < 0) continue;

        const TreeNode& node = tree.nodes[static_cast<size_t>(parent_id)];
        if (node.is_tip) continue;

        const int left_id = node.left;
        const int right_id = node.right;
        if (left_id < 0 || right_id < 0) continue;

        const int child_depth = parent_depth + 1;
        if (static_cast<size_t>(child_depth) >= level_counts.size()) {
            level_counts.resize(static_cast<size_t>(child_depth) + 1, 0);
        }
        level_counts[static_cast<size_t>(child_depth)] += 2;
        node_depth[static_cast<size_t>(left_id)] = child_depth;
        node_depth[static_cast<size_t>(right_id)] = child_depth;
    }

    level_offsets.assign(level_counts.size() + 1, 0);
    for (size_t depth = 0; depth < level_counts.size(); ++depth) {
        level_offsets[depth + 1] = level_offsets[depth] + level_counts[depth];
    }

    downward_ops_host.resize(static_cast<size_t>(level_offsets.back()));
    std::vector<int> write_offsets(level_offsets.begin(), level_offsets.end() - 1);

    for (int parent_id : tree.preorder) {
        if (parent_id < 0 || parent_id >= node_count) continue;
        const TreeNode& node = tree.nodes[static_cast<size_t>(parent_id)];
        if (node.is_tip) continue;

        const int left_id = node.left;
        const int right_id = node.right;
        if (left_id < 0 || right_id < 0) continue;

        const int child_depth = node_depth[static_cast<size_t>(left_id)];
        if (child_depth <= 0 || static_cast<size_t>(child_depth) >= write_offsets.size()) continue;

        const bool left_is_tip = tree.nodes[static_cast<size_t>(left_id)].is_tip;
        const bool right_is_tip = tree.nodes[static_cast<size_t>(right_id)].is_tip;

        downward_ops_host[static_cast<size_t>(write_offsets[static_cast<size_t>(child_depth)]++)] =
            make_downward_op(
                parent_id,
                left_id,
                right_id,
                left_is_tip,
                right_is_tip,
                node_to_tip,
                static_cast<uint8_t>(CLV_DIR_DOWN_LEFT));
        downward_ops_host[static_cast<size_t>(write_offsets[static_cast<size_t>(child_depth)]++)] =
            make_downward_op(
                parent_id,
                left_id,
                right_id,
                left_is_tip,
                right_is_tip,
                node_to_tip,
                static_cast<uint8_t>(CLV_DIR_DOWN_RIGHT));
    }
}

static void launch_upward_clv_update(
    const DeviceTree& D,
    NodeOpInfo* d_ops,
    int num_ops,
    int num_sms,
    cudaStream_t stream)
{
    if (num_ops <= 0 || !d_ops) {
        throw std::runtime_error("No upward ops to update");
    }
    launch_init_tip_clv(D, stream);

    struct LaunchConfig {
        int block = 256;
        int max_blocks_per_sm = 4;
        bool initialized = false;
    };
    static LaunchConfig cfg;
    if (!cfg.initialized) {
        cudaFuncAttributes attr{};
        CUDA_CHECK(cudaFuncGetAttributes(&attr, Rtree_Likelihood_Site_Parallel_Upward_Kernel));
        if (attr.maxThreadsPerBlock > 0 && cfg.block > attr.maxThreadsPerBlock) {
            cfg.block = attr.maxThreadsPerBlock;
        }
        CUDA_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &cfg.max_blocks_per_sm,
            Rtree_Likelihood_Site_Parallel_Upward_Kernel,
            cfg.block,
            0));
        cfg.initialized = true;
    }

    int max_blocks = num_sms * cfg.max_blocks_per_sm;
    int grid = static_cast<int>((D.sites + static_cast<size_t>(cfg.block) - 1) / static_cast<size_t>(cfg.block));
    if (max_blocks > 0 && grid > max_blocks) {
        grid = max_blocks;
    }

    Rtree_Likelihood_Site_Parallel_Upward_Kernel<<<grid, cfg.block, 0, stream>>>(
        D,
        d_ops,
        num_ops);
    CHECK_CUDA_LAST();
}

static void launch_downward_clv_update(
    const DeviceTree& D,
    NodeOpInfo* d_ops,
    int num_ops,
    int num_sms,
    cudaStream_t stream)
{
    if (num_ops <= 0 || !d_ops) return;

    struct LaunchConfig {
        int block = 256;
        int max_blocks_per_sm = 4;
        bool initialized = false;
    };
    static LaunchConfig cfg;
    if (!cfg.initialized) {
        cudaFuncAttributes attr{};
        CUDA_CHECK(cudaFuncGetAttributes(&attr, Rtree_Likelihood_Site_Parallel_Downward_Kernel));
        if (attr.maxThreadsPerBlock > 0 && cfg.block > attr.maxThreadsPerBlock) {
            cfg.block = attr.maxThreadsPerBlock;
        }
        CUDA_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &cfg.max_blocks_per_sm,
            Rtree_Likelihood_Site_Parallel_Downward_Kernel,
            cfg.block,
            0));
        cfg.initialized = true;
    }

    int max_blocks = num_sms * cfg.max_blocks_per_sm;
    int grid = static_cast<int>((D.sites + static_cast<size_t>(cfg.block) - 1) / static_cast<size_t>(cfg.block));
    if (max_blocks > 0 && grid > max_blocks) {
        grid = max_blocks;
    }

    Rtree_Likelihood_Site_Parallel_Downward_Kernel<<<grid, cfg.block, 0, stream>>>(D, d_ops, num_ops);
    CHECK_CUDA_LAST();
}

static void launch_downward_clv_update_levelized(
    const DeviceTree& D,
    NodeOpInfo* d_ops,
    const std::vector<int>& level_offsets,
    int num_sms,
    cudaStream_t stream)
{
    if (!d_ops || D.sites == 0 || level_offsets.size() < 2) return;

    struct LaunchConfig {
        int block = 256;
        int max_blocks_per_sm = 4;
        bool initialized = false;
    };
    static LaunchConfig cfg;
    if (!cfg.initialized) {
        cudaFuncAttributes attr{};
        CUDA_CHECK(cudaFuncGetAttributes(&attr, Rtree_Likelihood_Site_Parallel_Downward_Level_Kernel));
        if (attr.maxThreadsPerBlock > 0 && cfg.block > attr.maxThreadsPerBlock) {
            cfg.block = attr.maxThreadsPerBlock;
        }
        CUDA_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &cfg.max_blocks_per_sm,
            Rtree_Likelihood_Site_Parallel_Downward_Level_Kernel,
            cfg.block,
            0));
        cfg.initialized = true;
    }

    int max_blocks = num_sms * cfg.max_blocks_per_sm;
    unsigned int grid_x = static_cast<unsigned int>(
        (D.sites + static_cast<size_t>(cfg.block) - 1) / static_cast<size_t>(cfg.block));
    if (max_blocks > 0 && static_cast<int>(grid_x) > max_blocks) {
        grid_x = static_cast<unsigned int>(max_blocks);
    }
    if (grid_x == 0) {
        grid_x = 1;
    }

    constexpr int kMaxGridY = 65535;
    for (size_t level = 1; level + 1 < level_offsets.size(); ++level) {
        const int level_start = level_offsets[level];
        const int level_count = level_offsets[level + 1] - level_start;
        if (level_count <= 0) continue;

        for (int batch_start = 0; batch_start < level_count; batch_start += kMaxGridY) {
            const int batch_ops = std::min(level_count - batch_start, kMaxGridY);
            dim3 grid(grid_x, static_cast<unsigned int>(batch_ops));
            Rtree_Likelihood_Site_Parallel_Downward_Level_Kernel<<<grid, cfg.block, 0, stream>>>(
                D,
                d_ops + level_start + batch_start,
                batch_ops);
            CUDA_CHECK(cudaGetLastError());
        }
    }
}

static int get_device_sm_count()
{
    static int cached_device = -1;
    static int sm_count = 0;
    const int current_device = mlipper::gpu::current_device_or_throw();
    if (sm_count <= 0 || cached_device != current_device) {
        const cudaDeviceProp device_props =
            mlipper::gpu::current_device_properties_or_throw();
        cached_device = current_device;
        sm_count = device_props.multiProcessorCount;
    }
    return sm_count;
}

static void ensure_placement_op_capacity(
    PlacementOpBuffer& placement_ops,
    int required_ops,
    cudaStream_t stream)
{
    if (required_ops <= 0) return;
    if (placement_ops.capacity >= required_ops && placement_ops.d_ops) return;

    if (placement_ops.d_ops) {
        CUDA_CHECK(cudaStreamSynchronize(stream));
        CUDA_CHECK(cudaFree(placement_ops.d_ops));
        placement_ops.d_ops = nullptr;
    }

    CUDA_CHECK(cudaMalloc(
        &placement_ops.d_ops,
        sizeof(NodeOpInfo) * static_cast<size_t>(required_ops)));
    placement_ops.capacity = required_ops;
}

static void upload_ops_to_device(
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    cudaStream_t stream)
{
    const int num_ops = static_cast<int>(host_ops.size());
    if (num_ops <= 0) {
        placement_ops.num_ops = 0;
        return;
    }

    ensure_placement_op_capacity(placement_ops, num_ops, stream);
    CUDA_CHECK(cudaMemcpyAsync(
        placement_ops.d_ops,
        host_ops.data(),
        sizeof(NodeOpInfo) * static_cast<size_t>(num_ops),
        cudaMemcpyHostToDevice,
        stream));
    placement_ops.num_ops = num_ops;
}

template <typename LaunchFn>
static void run_clv_stage(
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    cudaStream_t stream,
    double* stage_ms,
    LaunchFn&& launch_fn)
{
    const auto stage = [&]() {
        upload_ops_to_device(placement_ops, host_ops, stream);
        launch_fn();
    };
    if (stage_ms) {
        *stage_ms += time_stream_stage_ms(stream, stage);
    } else {
        stage();
    }
}

static void run_upward_clv_stage(
    const DeviceTree& D,
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    int num_sms,
    cudaStream_t stream,
    double* stage_ms = nullptr)
{
    run_clv_stage(
        placement_ops,
        host_ops,
        stream,
        stage_ms,
        [&]() {
            launch_upward_clv_update(
                D,
                placement_ops.d_ops,
                placement_ops.num_ops,
                num_sms,
                stream);
        });
}

static void run_downward_clv_stage(
    const DeviceTree& D,
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    int num_sms,
    cudaStream_t stream,
    double* stage_ms = nullptr)
{
    run_clv_stage(
        placement_ops,
        host_ops,
        stream,
        stage_ms,
        [&]() {
            launch_downward_clv_update(
                D,
                placement_ops.d_ops,
                placement_ops.num_ops,
                num_sms,
                stream);
        });
}

static void run_downward_clv_stage_levelized(
    const DeviceTree& D,
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    const std::vector<int>& level_offsets,
    int num_sms,
    cudaStream_t stream,
    double* stage_ms = nullptr)
{
    run_clv_stage(
        placement_ops,
        host_ops,
        stream,
        stage_ms,
        [&]() {
            launch_downward_clv_update_levelized(
                D,
                placement_ops.d_ops,
                level_offsets,
                num_sms,
                stream);
        });
}

void UploadPlacementOps(
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    cudaStream_t stream)
{
    upload_ops_to_device(placement_ops, host_ops, stream);
}

void launch_init_tip_clv(const DeviceTree& D, cudaStream_t stream)
{
    if (D.tips <= 0 || D.sites == 0 || !D.d_tipchars || !D.d_tip_node_ids || !D.d_clv_up) {
        return;
    }

    const size_t total_tip_sites = static_cast<size_t>(D.tips) * D.sites;
    dim3 block(256);
    dim3 grid(static_cast<unsigned>((total_tip_sites + block.x - 1) / block.x));
    InitializeTipClvUpKernel<<<grid, block, 0, stream>>>(D);
    CUDA_CHECK(cudaGetLastError());
}

void free_placement_op_buffer(
    PlacementOpBuffer& placement_ops,
    cudaStream_t stream)
{
    if (placement_ops.d_ops) {
        CUDA_CHECK(cudaStreamSynchronize(stream));
        CUDA_CHECK(cudaFree(placement_ops.d_ops));
        placement_ops.d_ops = nullptr;
    }
    placement_ops.num_ops = 0;
    placement_ops.capacity = 0;
}

void DownloadClvDump(
    const DeviceTree& D,
    std::vector<fp_t>& clv_up,
    std::vector<unsigned>& scaler_up,
    cudaStream_t stream)
{
    clv_up.clear();
    scaler_up.clear();
    if (D.N <= 0) return;
    const size_t per_node = D.per_node_elems();
    if (per_node == 0) return;
    const size_t total = static_cast<size_t>(D.N) * per_node;
    if (!D.d_clv_up) {
        throw std::runtime_error("DownloadClvDump: missing d_clv_up.");
    }

    clv_up.resize(total);
    CUDA_CHECK(cudaMemcpyAsync(
        clv_up.data(),
        D.d_clv_up,
        sizeof(fp_t) * total,
        cudaMemcpyDeviceToHost,
        stream));

    const size_t scaler_span = D.per_rate_scaling
        ? static_cast<size_t>(D.sites) * static_cast<size_t>(D.rate_cats)
        : static_cast<size_t>(D.sites);
    if (D.d_site_scaler_up && scaler_span > 0) {
        const size_t scaler_total = static_cast<size_t>(D.N) * scaler_span;
        scaler_up.resize(scaler_total);
        CUDA_CHECK(cudaMemcpyAsync(
            scaler_up.data(),
            D.d_site_scaler_up,
            sizeof(unsigned) * scaler_total,
            cudaMemcpyDeviceToHost,
            stream));
    }

    CUDA_CHECK(cudaStreamSynchronize(stream));
}

// ----- Tree CLV update pipeline -----

// Full-tree CLV rebuild used after the initial topology/device upload.
void UpdateTreeClvs(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    cudaStream_t stream)
{
    const int sm_count = get_device_sm_count();
    const bool profile = placement_ops.profile_commit_timing;

    const auto upward_host_start = SteadyClock::now();
    rebuild_node_to_tip_map(T, H, placement_ops.node_to_tip);
    build_upward_ops_host(T, placement_ops.node_to_tip, placement_ops.upward_ops_host);
    const auto upward_host_end = SteadyClock::now();
    if (profile) {
        placement_ops.timing.initial_upward_host_ms += elapsed_ms(upward_host_start, upward_host_end);
        placement_ops.timing.initial_upward_ops += static_cast<long long>(placement_ops.upward_ops_host.size());
    }

    if (profile) {
        run_upward_clv_stage(
            D,
            placement_ops,
            placement_ops.upward_ops_host,
            sm_count,
            stream,
            &placement_ops.timing.initial_upward_stage_ms);
    } else {
        run_upward_clv_stage(
            D,
            placement_ops,
            placement_ops.upward_ops_host,
            sm_count,
            stream);
    }

    const auto downward_host_start = SteadyClock::now();
    build_downward_ops_host_levelized(
        T,
        placement_ops.node_to_tip,
        placement_ops.downward_ops_host,
        placement_ops.downward_level_offsets_host);
    const auto downward_host_end = SteadyClock::now();
    if (profile) {
        placement_ops.timing.initial_downward_host_ms += elapsed_ms(downward_host_start, downward_host_end);
        placement_ops.timing.initial_downward_ops += static_cast<long long>(placement_ops.downward_ops_host.size());
    }

    if (profile) {
        run_downward_clv_stage_levelized(
            D,
            placement_ops,
            placement_ops.downward_ops_host,
            placement_ops.downward_level_offsets_host,
            sm_count,
            stream,
            &placement_ops.timing.initial_downward_stage_ms);
        placement_ops.timing.initial_updates += 1;
    } else {
        run_downward_clv_stage_levelized(
            D,
            placement_ops,
            placement_ops.downward_ops_host,
            placement_ops.downward_level_offsets_host,
            sm_count,
            stream);
    }
}

// Local SPR uses this after pruning to refresh only the affected CLV region.
// The upward pass is limited to the path from upward_start_node to the root;
// passing a negative start skips upward work and applies only downward updates.
void UpdateTreeClvsAfterPrune(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    int upward_start_node,
    const std::vector<NodeOpInfo>& required_downward_ops,
    cudaStream_t stream)
{
    const int sm_count = get_device_sm_count();
    rebuild_node_to_tip_map(T, H, placement_ops.node_to_tip);

    build_upward_ops_host_for_path(
        T,
        placement_ops.node_to_tip,
        upward_start_node,
        placement_ops.upward_ops_host);
    if (!placement_ops.upward_ops_host.empty()) {
        run_upward_clv_stage(
            D,
            placement_ops,
            placement_ops.upward_ops_host,
            sm_count,
            stream);
    }

    if (required_downward_ops.empty()) {
        placement_ops.num_ops = 0;
        return;
    }

    // Local SPR already supplies the exact dependent downward ops to refresh,
    // so this path keeps the simple sequential launch instead of rebuilding
    // level batches from the full tree topology.
    run_downward_clv_stage(
        D,
        placement_ops,
        required_downward_ops,
        sm_count,
        stream);
}

// After a committed placement inserts a new internal node and query tip, this
// incrementally refreshes the modified upward path and then rebuilds downward CLVs.
static void UpdateTreeClvsAfterInsertion(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    int upward_start_node,
    cudaStream_t stream)
{
    const int sm_count = get_device_sm_count();
    const bool profile = placement_ops.profile_commit_timing;

    const auto upward_host_start = SteadyClock::now();
    rebuild_node_to_tip_map(T, H, placement_ops.node_to_tip);
    build_upward_ops_host_for_path(
        T,
        placement_ops.node_to_tip,
        upward_start_node,
        placement_ops.upward_ops_host);
    const auto upward_host_end = SteadyClock::now();
    if (profile) {
        placement_ops.timing.insertion_upward_host_ms += elapsed_ms(upward_host_start, upward_host_end);
        placement_ops.timing.insertion_upward_ops += static_cast<long long>(placement_ops.upward_ops_host.size());
    }

    if (!placement_ops.upward_ops_host.empty()) {
        if (profile) {
            run_upward_clv_stage(
                D,
                placement_ops,
                placement_ops.upward_ops_host,
                sm_count,
                stream,
                &placement_ops.timing.insertion_upward_stage_ms);
        } else {
            run_upward_clv_stage(
                D,
                placement_ops,
                placement_ops.upward_ops_host,
                sm_count,
                stream);
        }
    }

    const auto downward_host_start = SteadyClock::now();
    build_downward_ops_host_levelized(
        T,
        placement_ops.node_to_tip,
        placement_ops.downward_ops_host,
        placement_ops.downward_level_offsets_host);
    const auto downward_host_end = SteadyClock::now();
    if (profile) {
        placement_ops.timing.insertion_downward_host_ms += elapsed_ms(downward_host_start, downward_host_end);
        placement_ops.timing.insertion_downward_ops += static_cast<long long>(placement_ops.downward_ops_host.size());
    }

    if (profile) {
        run_downward_clv_stage_levelized(
            D,
            placement_ops,
            placement_ops.downward_ops_host,
            placement_ops.downward_level_offsets_host,
            sm_count,
            stream,
            &placement_ops.timing.insertion_downward_stage_ms);
        placement_ops.timing.insertion_updates += 1;
    } else {
        run_downward_clv_stage_levelized(
            D,
            placement_ops,
            placement_ops.downward_ops_host,
            placement_ops.downward_level_offsets_host,
            sm_count,
            stream);
    }
}

// ----- Query placement evaluation -----
static InsertResult insert_query_with_intermediate(
    TreeBuildResult& T,
    const std::string& raw_name,
    int target_id,
    double pendant,
    double proximal)
{
    InsertResult out{};
    if (target_id < 0 || target_id >= (int)T.nodes.size()) {
        throw std::runtime_error("insert_query_with_intermediate: invalid target_id.");
    }
    int parent_id = T.nodes[target_id].parent;
    const double total = T.nodes[target_id].branch_length_to_parent;
    double proximal_len = 0.0;
    double distal_len = 0.0;
    normalize_split_branch_lengths(total, proximal, OPT_BRANCH_LEN_MIN, proximal_len, distal_len);
    double pendant_len = sanitize_branch_length(pendant);

    int new_internal_id = (int)T.nodes.size();
    int new_tip_id = new_internal_id + 1;

    std::string name = raw_name.empty()
        ? ("query_" + std::to_string(new_tip_id))
        : raw_name;
    if (T.tip_node_by_name.count(name)) {
        name += "_" + std::to_string(new_tip_id);
    }

    TreeNode internal{};
    internal.id = new_internal_id;
    internal.is_tip = false;
    internal.parent = parent_id;
    internal.left = target_id;
    internal.right = new_tip_id;
    internal.branch_length_to_parent = proximal_len;

    TreeNode tip{};
    tip.id = new_tip_id;
    tip.is_tip = true;
    tip.parent = new_internal_id;
    tip.left = -1;
    tip.right = -1;
    tip.branch_length_to_parent = pendant_len;
    tip.name = name;

    if (parent_id >= 0) {
        if (T.nodes[parent_id].left == target_id) {
            T.nodes[parent_id].left = new_internal_id;
        } else if (T.nodes[parent_id].right == target_id) {
            T.nodes[parent_id].right = new_internal_id;
        } else {
            throw std::runtime_error(
                "insert_query_with_intermediate: parent does not reference target.");
        }
    } else {
        // Target was root; new internal becomes root.
        T.root_id = new_internal_id;
    }

    T.nodes[target_id].parent = new_internal_id;
    T.nodes[target_id].branch_length_to_parent = distal_len;

    T.nodes.push_back(internal);
    T.nodes.push_back(tip);
    T.tip_node_by_name[name] = new_tip_id;

    rebuild_traversals(T);
    out.internal_id = new_internal_id;
    out.tip_id = new_tip_id;
    out.tip_name = name;
    return out;
}

static void rebuild_host_topology_from_tree(const TreeBuildResult& T, HostPacking& H)
{
    const int N = (int)T.nodes.size();
    H.postorder = T.postorder;
    H.preorder  = T.preorder;
    H.parent.assign(N, -1);
    H.left.assign(N, -1);
    H.right.assign(N, -1);
    H.is_tip.assign(N, 0);
    H.blen.assign(N, 0.0);
    for (int i = 0; i < N; ++i) {
        const TreeNode& nd = T.nodes[i];
        H.parent[i] = nd.parent;
        H.left[i]   = nd.left;
        H.right[i]  = nd.right;
        H.is_tip[i] = nd.is_tip;
        H.blen[i]   = nd.branch_length_to_parent;
    }
}
static void append_query_tip_to_host_packing(
    HostPacking& H,
    const PlacementQueryBatch& Q,
    int query_idx,
    int new_tip_node_id,
    size_t sites)
{
    if (query_idx < 0 || (size_t)query_idx >= Q.count) {
        throw std::runtime_error("append_query_tip_to_host_packing: query_idx out of range.");
    }
    const size_t needed = ((size_t)query_idx + 1) * sites;
    if (Q.query_chars.size() < needed) {
        throw std::runtime_error("append_query_tip_to_host_packing: query_chars buffer too small.");
    }

    H.tip_node_ids.push_back(new_tip_node_id);
    const size_t old = H.tipchars.size();
    H.tipchars.resize(old + sites);
    std::memcpy(
        H.tipchars.data() + old,
        Q.query_chars.data() + (size_t)query_idx * sites,
        sites * sizeof(uint8_t));
}

static void update_insertion_device(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    size_t sites,
    int states,
    int rate_cats,
    int target_id,
    int internal_id,
    int tip_id,
    cudaStream_t stream)
{
    const int old_tips = D.tips;
    const int new_tips = (int)H.tip_node_ids.size();
    const int newN = (int)T.nodes.size();
    if (newN > D.capacity_N) {
        throw std::runtime_error("update_insertion_device: node capacity exceeded.");
    }
    if (new_tips > D.capacity_tips) {
        throw std::runtime_error("update_insertion_device: tip capacity exceeded.");
    }

    D.N = newN;
    D.tips = new_tips;
    D.inners = D.N - D.tips;
    D.root_id = T.root_id;

    CUDA_CHECK(cudaMemcpyAsync(
        D.d_tipchars + (size_t)old_tips * sites,
        H.tipchars.data() + (size_t)old_tips * sites,
        sizeof(uint8_t) * sites,
        cudaMemcpyHostToDevice,
        stream));
    CUDA_CHECK(cudaMemcpyAsync(
        D.d_tip_node_ids + (size_t)old_tips,
        H.tip_node_ids.data() + (size_t)old_tips,
        sizeof(int),
        cudaMemcpyHostToDevice,
        stream));

    const size_t stride = (size_t)rate_cats * (size_t)states * (size_t)states;
    const fp_t* mid_prox_src = H.pmats_mid_prox.empty() ? H.pmats_mid.data() : H.pmats_mid_prox.data();
    const fp_t* mid_dist_src = H.pmats_mid_dist.empty() ? H.pmats_mid.data() : H.pmats_mid_dist.data();

    auto sync_node = [&](int nid) {
        if (nid < 0 || nid >= newN) {
            throw std::runtime_error("update_insertion_device: node id out of range.");
        }
        const size_t off = (size_t)nid * stride;
        CUDA_CHECK(cudaMemcpyAsync(
            D.d_blen + (size_t)nid,
            &H.blen[(size_t)nid],
            sizeof(fp_t),
            cudaMemcpyHostToDevice,
            stream));

        CUDA_CHECK(cudaMemcpyAsync(
            D.d_pmat + off,
            H.pmats.data() + off,
            sizeof(fp_t) * stride,
            cudaMemcpyHostToDevice,
            stream));
        CUDA_CHECK(cudaMemcpyAsync(
            D.d_pmat_mid + off,
            H.pmats_mid.data() + off,
            sizeof(fp_t) * stride,
            cudaMemcpyHostToDevice,
            stream));
        CUDA_CHECK(cudaMemcpyAsync(
            D.d_pmat_mid_prox + off,
            mid_prox_src + off,
            sizeof(fp_t) * stride,
            cudaMemcpyHostToDevice,
            stream));
        CUDA_CHECK(cudaMemcpyAsync(
            D.d_pmat_mid_dist + off,
            mid_dist_src + off,
            sizeof(fp_t) * stride,
            cudaMemcpyHostToDevice,
            stream));
    };

    sync_node(target_id);
    sync_node(internal_id);
    sync_node(tip_id);
}
static std::string query_name_for_commit(
    const PlacementCommitContext& commit_ctx,
    int query_idx)
{
    if (!commit_ctx.query_names) return std::string{};
    if (query_idx < 0 || query_idx >= static_cast<int>(commit_ctx.query_names->size())) {
        return std::string{};
    }
    return (*commit_ctx.query_names)[query_idx];
}

static void validate_commit_context(const PlacementCommitContext& commit_ctx)
{
    if (!commit_ctx.tree) {
        throw std::runtime_error("EvaluatePlacementQueries: commit tree is null.");
    }
    if (!commit_ctx.host) {
        throw std::runtime_error("EvaluatePlacementQueries: commit host packing is null.");
    }
    if (!commit_ctx.queries) {
        throw std::runtime_error("EvaluatePlacementQueries: commit query batch is null.");
    }
    if (!commit_ctx.placement_ops) {
        throw std::runtime_error("EvaluatePlacementQueries: commit placement ops are null.");
    }
}

static void commit_placement_result(
    DeviceTree& D,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    PlacementCommitContext& commit_ctx,
    const PlacementResult& placement,
    int query_idx,
    cudaStream_t stream)
{
    if (placement.target_id < 0) {
        throw std::runtime_error("commit_placement_result: invalid placement target.");
    }

    TreeBuildResult& tree = *commit_ctx.tree;
    HostPacking& host = *commit_ctx.host;
    PlacementQueryBatch& queries = *commit_ctx.queries;
    PlacementOpBuffer& placement_ops = *commit_ctx.placement_ops;
    const bool profile = placement_ops.profile_commit_timing;
    const auto pre_clv_start = SteadyClock::now();

    const double total_branch_length = tree.nodes[placement.target_id].branch_length_to_parent;
    // PlacementEvaluationKernel currently stores the jplace distal coordinate
    // in placement.proximal_length. Tree insertion expects the parent-side
    // split length, so convert distal -> proximal before mutating the tree.
    const double commit_proximal_length = total_branch_length - placement.proximal_length;

    const InsertResult insert_result = insert_query_with_intermediate(
        tree,
        query_name_for_commit(commit_ctx, query_idx),
        placement.target_id,
        placement.pendant_length,
        commit_proximal_length);
    if (insert_result.internal_id < 0 || insert_result.tip_id < 0) {
        throw std::runtime_error("commit_placement_result: insertion failed.");
    }
    if (commit_ctx.inserted_query_names &&
        query_idx >= 0 &&
        query_idx < static_cast<int>(commit_ctx.inserted_query_names->size())) {
        (*commit_ctx.inserted_query_names)[query_idx] = insert_result.tip_name;
    }

    rebuild_host_topology_from_tree(tree, host);
    append_query_tip_to_host_packing(host, queries, query_idx, insert_result.tip_id, D.sites);

    const int changed_nodes[3] = {
        placement.target_id,
        insert_result.internal_id,
        insert_result.tip_id,
    };
    fill_pmats_in_host_packing(
        tree,
        host,
        er,
        rate_multipliers,
        D.states,
        D.rate_cats,
        changed_nodes,
        3);

    update_insertion_device(
        D,
        tree,
        host,
        D.sites,
        D.states,
        D.rate_cats,
        placement.target_id,
        insert_result.internal_id,
        insert_result.tip_id,
        stream);

    if (profile) {
        const auto pre_clv_end = SteadyClock::now();
        placement_ops.timing.insertion_pre_clv_ms += elapsed_ms(pre_clv_start, pre_clv_end);
    }

    UpdateTreeClvsAfterInsertion(
        D,
        tree,
        host,
        placement_ops,
        insert_result.internal_id,
        stream);
}

static int resolve_query_count(const DeviceTree& D)
{
    int query_count = D.placement_queries;
    if (const char* env_max_queries = std::getenv("MLIPPER_MAX_QUERIES")) {
        const int parsed = std::atoi(env_max_queries);
        if (parsed > 0) {
            query_count = std::min(D.placement_queries, parsed);
        }
    }
    return query_count;
}

static void reset_placement_results(
    std::vector<PlacementResult>* placement_results_out,
    int query_count)
{
    if (!placement_results_out) return;
    placement_results_out->clear();
    placement_results_out->reserve(static_cast<size_t>(query_count));
}

static PlacementResult evaluate_single_placement_query(
    DeviceTree& D,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    const PlacementOpBuffer& placement_ops,
    int query_idx,
    int smoothing,
    cudaStream_t stream)
{
    const size_t node_bytes = sizeof(fp_t) * static_cast<size_t>(D.N);
    if (placement_ops.profile_commit_timing) {
        PlacementOpBuffer& timing_ops = const_cast<PlacementOpBuffer&>(placement_ops);
        timing_ops.timing.query_evals += 1;
        timing_ops.timing.query_reset_stage_ms += time_stream_stage_ms(stream, [&]() {
            CUDA_CHECK(cudaMemset(D.d_new_pendant_length, 0, node_bytes));
            CUDA_CHECK(cudaMemset(D.d_new_proximal_length, 0, node_bytes));
            CUDA_CHECK(cudaMemset(D.d_prev_pendant_length, 0, node_bytes));
            CUDA_CHECK(cudaMemset(D.d_prev_proximal_length, 0, node_bytes));
        });
        timing_ops.timing.query_build_clv_stage_ms += time_stream_stage_ms(stream, [&]() {
            build_query_clv(D, query_idx, stream);
            CHECK_CUDA_LAST();
        });
        const auto kernel_start = SteadyClock::now();
        DeviceTree query_view = make_query_view(D, query_idx);
        PlacementResult result = PlacementEvaluationKernel(
            query_view,
            placement_ops.d_ops,
            placement_ops.num_ops,
            smoothing,
            stream,
            true);
        timing_ops.timing.query_kernel_total_ms +=
            elapsed_ms(kernel_start, SteadyClock::now());
        return result;
    }

    CUDA_CHECK(cudaMemset(D.d_new_pendant_length, 0, node_bytes));
    CUDA_CHECK(cudaMemset(D.d_new_proximal_length, 0, node_bytes));
    CUDA_CHECK(cudaMemset(D.d_prev_pendant_length, 0, node_bytes));
    CUDA_CHECK(cudaMemset(D.d_prev_proximal_length, 0, node_bytes));

    build_query_clv(D, query_idx, stream);
    CHECK_CUDA_LAST();

    DeviceTree query_view = make_query_view(D, query_idx);
    return PlacementEvaluationKernel(
        query_view,
        placement_ops.d_ops,
        placement_ops.num_ops,
        smoothing,
        stream,
        true);
}

static void evaluate_queries(
    DeviceTree& D,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    const PlacementOpBuffer& placement_ops,
    std::vector<PlacementResult>* placement_results_out,
    int smoothing,
    int query_count,
    bool commit_to_tree,
    PlacementCommitContext* commit_ctx,
    cudaStream_t stream)
{
    if (commit_to_tree && !commit_ctx) {
        throw std::runtime_error("evaluate_queries: commit_ctx is null but commit_to_tree is true.");
    }
    if (commit_to_tree) {
        validate_commit_context(*commit_ctx);
    }

    for (int query_idx = 0; query_idx < query_count; ++query_idx) {
        PlacementResult placement = evaluate_single_placement_query(
            D,
            er,
            rate_multipliers,
            placement_ops,
            query_idx,
            smoothing,
            stream);
        if (placement_results_out) {
            placement_results_out->push_back(placement);
        }
        if (commit_to_tree) {
            commit_placement_result(
                D,
                er,
                rate_multipliers,
                *commit_ctx,
                placement,
                query_idx,
                stream);
        }
    }
}

void EvaluatePlacementQueries(
    DeviceTree& D,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    PlacementCommitContext& commit_ctx,
    std::vector<PlacementResult>* placement_results_out,
    int smoothing,
    bool commit_to_tree,
    cudaStream_t stream)
{
    const int query_count = resolve_query_count(D);
    reset_placement_results(placement_results_out, query_count);

    if (!commit_ctx.placement_ops) {
        throw std::runtime_error("EvaluatePlacementQueries: placement ops are null.");
    }

    evaluate_queries(
        D,
        er,
        rate_multipliers,
        *commit_ctx.placement_ops,
        placement_results_out,
        smoothing,
        query_count,
        commit_to_tree,
        commit_to_tree ? &commit_ctx : nullptr,
        stream);
}
