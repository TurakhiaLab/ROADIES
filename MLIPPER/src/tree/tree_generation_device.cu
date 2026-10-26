#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cuda_runtime.h>
#include <tbb/parallel_for.h>

#include "gpu/device_buffer.hpp"
#include "gpu/gpu_admission.hpp"
#include "likelihood/partials.cuh"
#include "pmat.h"
#include "placement/placement.cuh"
#include "tree.hpp"
#include "tree_topology_utils.hpp"
#include "util/checked_size.hpp"
#include "util/mlipper_util.h"

namespace {

static bool supports_dna_g4_fast_path(const DeviceTree& D)
{
    return D.states == 4 && D.rate_cats == 4;
}

static void copy_host_to_device_async(void* dst, const void* src, size_t bytes, cudaStream_t stream)
{
    CUDA_CHECK(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, stream));
}

static void copy_host_to_device(void* dst, const void* src, size_t bytes)
{
    CUDA_CHECK(cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice));
}

static void copy_device_to_device_async(void* dst, const void* src, size_t bytes, cudaStream_t stream)
{
    CUDA_CHECK(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToDevice, stream));
}

template <typename T>
static void cuda_malloc_bytes(T*& ptr, size_t bytes)
{
    const cudaError_t status =
        cudaMalloc(reinterpret_cast<void**>(&ptr), bytes);
    if (status != cudaSuccess) {
        size_t free_bytes = 0;
        size_t total_bytes = 0;
        cudaMemGetInfo(&free_bytes, &total_bytes);
        throw CudaRuntimeError(
            status,
            "[CUDA] cudaMalloc failed: requested_bytes=" +
            std::to_string(bytes) +
            " free_bytes=" + std::to_string(free_bytes) +
            " total_bytes=" + std::to_string(total_bytes) +
            " error=" + cudaGetErrorString(status));
    }
}

template <typename T>
static void cuda_malloc_zeroed(T*& ptr, size_t bytes)
{
    cuda_malloc_bytes(ptr, bytes);
    CUDA_CHECK(cudaMemset(ptr, 0, bytes));
}

template <typename T>
static void cuda_free_if_allocated(T*& ptr) noexcept
{
    if (ptr != nullptr) {
        cudaFree(ptr);
        ptr = nullptr;
    }
}

static void zero_length_scratch_async(const DeviceTree& D, size_t bytes, cudaStream_t stream)
{
    CUDA_CHECK(cudaMemsetAsync(D.d_new_pendant_length, 0, bytes, stream));
    CUDA_CHECK(cudaMemsetAsync(D.d_new_proximal_length, 0, bytes, stream));
    CUDA_CHECK(cudaMemsetAsync(D.d_prev_pendant_length, 0, bytes, stream));
    CUDA_CHECK(cudaMemsetAsync(D.d_prev_proximal_length, 0, bytes, stream));
}

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

struct LaunchConfig {
    int block = 256, max_blocks_per_sm = 4;
    int device = -1;
};

template <typename KernelFn>
static LaunchConfig& initialized_launch_config(LaunchConfig& cfg, KernelFn kernel)
{
    const int current_device = mlipper::gpu::current_device_or_throw();
    if (cfg.device != current_device) {
        cfg.block = 256;
        cfg.max_blocks_per_sm = 4;
        cudaFuncAttributes attr{};
        CUDA_CHECK(cudaFuncGetAttributes(&attr, kernel));
        if (attr.maxThreadsPerBlock > 0 && cfg.block > attr.maxThreadsPerBlock) {
            cfg.block = attr.maxThreadsPerBlock;
        }
        CUDA_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &cfg.max_blocks_per_sm, kernel, cfg.block, 0));
        cfg.device = current_device;
    }
    return cfg;
}

static int capped_site_grid(size_t sites, const LaunchConfig& cfg, int num_sms)
{
    int grid = static_cast<int>(
        (sites + static_cast<size_t>(cfg.block) - 1) /
        static_cast<size_t>(cfg.block));
    const int max_blocks = num_sms * cfg.max_blocks_per_sm;
    return (max_blocks > 0 && grid > max_blocks) ? max_blocks : grid;
}

static size_t clv_pool_elements(const DeviceTree& D)
{
    return mlipper::util::checked_mul_size(
        static_cast<size_t>(D.capacity_N), D.per_node_elems(),
        "DeviceTree CLV pool");
}

// Upward and downward CLVs share one allocation. Scalers similarly occupy four
// consecutive pools; the individual DeviceTree pointers are borrowed slices.
static void configure_resident_buffer_views(DeviceTree& D)
{
    D.d_clv_down = D.d_clv_up != nullptr
        ? D.d_clv_up + clv_pool_elements(D)
        : nullptr;
    if (D.d_scaler_storage == nullptr) {
        D.d_site_scaler_up = nullptr;
        D.d_site_scaler_down = nullptr;
        D.d_edge_midpoint_scaler = nullptr;
        D.d_edge_outside_scaler = nullptr;
        return;
    }

    const size_t scaler_pool = D.scaler_pool_elems();
    D.d_site_scaler_up = D.d_scaler_storage;
    D.d_site_scaler_down = D.d_scaler_storage + scaler_pool;
    D.d_edge_midpoint_scaler = D.d_scaler_storage + scaler_pool * 2;
    D.d_edge_outside_scaler = D.d_scaler_storage + scaler_pool * 3;
}

static void validate_resident_buffer_views(
    const DeviceTree& D,
    const char* context)
{
    const fp_t* expected_clv_down = D.d_clv_up != nullptr
        ? D.d_clv_up + clv_pool_elements(D)
        : nullptr;
    if (D.d_clv_down != expected_clv_down) {
        throw std::runtime_error(
            std::string(context) +
            ": inconsistent upward/downward CLV layout.");
    }

    const size_t scaler_pool = D.scaler_pool_elems();
    const unsigned* expected_scaler_down = D.d_scaler_storage != nullptr
        ? D.d_scaler_storage + scaler_pool
        : nullptr;
    const unsigned* expected_midpoint = D.d_scaler_storage != nullptr
        ? D.d_scaler_storage + scaler_pool * 2
        : nullptr;
    const unsigned* expected_outside = D.d_scaler_storage != nullptr
        ? D.d_scaler_storage + scaler_pool * 3
        : nullptr;
    if (D.d_site_scaler_up != D.d_scaler_storage ||
        D.d_site_scaler_down != expected_scaler_down ||
        D.d_edge_midpoint_scaler != expected_midpoint ||
        D.d_edge_outside_scaler != expected_outside) {
        throw std::runtime_error(
            std::string(context) + ": inconsistent scaler slice layout.");
    }
}

static void allocate_query_buffers(
    const DeviceTree& D,
    int query_capacity,
    uint8_t*& query_chars,
    fp_t*& query_clv)
{
    query_chars = nullptr;
    query_clv = nullptr;
    if (query_capacity == 0) return;

    const size_t capacity = static_cast<size_t>(query_capacity);
    try {
        cuda_malloc_zeroed(
            query_chars,
            mlipper::util::checked_product(
                "query character buffer", sizeof(uint8_t), capacity, D.sites));
        cuda_malloc_zeroed(
            query_clv,
            mlipper::util::checked_product(
                "query CLV buffer", sizeof(fp_t), capacity, D.per_node_elems()));
    } catch (...) {
        cuda_free_if_allocated(query_chars);
        cuda_free_if_allocated(query_clv);
        throw;
    }
}

} // namespace

static void release_device_tree_buffers(DeviceTree& device_tree) noexcept;

struct InsertResult {
    int internal_id = -1;
    int tip_id = -1;
    std::string tip_name;
};

void fill_pmats_in_host_packing(
    const TreeBuildResult&       T,
    HostPacking&                 H,
    const EigResult&             er,
    const std::vector<double>&   rate_multipliers,   // len = rate_cats (per-category rate multipliers)
    int states,
    int rate_cats,
    const int* changed_nodes,
    int num_changed_nodes,
    bool include_midpoint_pmats)
{
    if (states <= 0 || rate_cats <= 0 ||
        rate_multipliers.size() < static_cast<size_t>(rate_cats)) {
        throw std::invalid_argument(
            "fill_pmats_in_host_packing requires positive dimensions and "
            "one multiplier per rate category");
    }
    const int N = (int)T.nodes.size();
    const size_t per_node = mlipper::util::checked_product(
        "host transition matrix per-node size",
        static_cast<size_t>(rate_cats), static_cast<size_t>(states),
        static_cast<size_t>(states));

    const size_t required = mlipper::util::checked_mul_size(
        static_cast<size_t>(N), per_node, "host transition matrix storage");
    const bool want_incremental = (changed_nodes && num_changed_nodes > 0);

    auto reset_all = [&]() {
        H.pmats.assign(required, fp_t(0));
        if (include_midpoint_pmats) {
            H.pmats_mid.assign(required, fp_t(0));
            H.pmats_mid_prox.assign(required, fp_t(0));
            H.pmats_mid_dist.assign(required, fp_t(0));
        }
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
        reset_all();
    } else {
        H.pmats.resize(required, fp_t(0));
        if (include_midpoint_pmats) {
            H.pmats_mid.resize(required, fp_t(0));
            if (!H.pmats_mid_prox.empty()) H.pmats_mid_prox.resize(required, fp_t(0));
            if (!H.pmats_mid_dist.empty()) H.pmats_mid_dist.resize(required, fp_t(0));
        }
    }

    auto compute_node_pmats = [&](int nid) {
        if (nid < 0 || nid >= N) return;
        const TreeNode& nd = T.nodes[nid];
        if (nd.parent < 0) return;
        fp_t* base = H.pmats.data() + (size_t)nid * per_node;
        fp_t* base_mid = include_midpoint_pmats
            ? H.pmats_mid.data() + (size_t)nid * per_node : nullptr;
        fp_t* base_mid_prox = include_midpoint_pmats && !H.pmats_mid_prox.empty()
            ? H.pmats_mid_prox.data() + (size_t)nid * per_node : nullptr;
        fp_t* base_mid_dist = include_midpoint_pmats && !H.pmats_mid_dist.empty()
            ? H.pmats_mid_dist.data() + (size_t)nid * per_node : nullptr;
        const double blen  = static_cast<double>(nd.branch_length_to_parent);
        std::vector<double> pbuf((size_t)states * (size_t)states);
        std::vector<double> pbuf_mid((size_t)states * (size_t)states);

        for (int rc = 0; rc < rate_cats; ++rc) {
            double r = rate_multipliers[rc];  // rate category multiplier
            double t = blen;                  // branch length
            double p = 0;

            fp_t* P = base + (size_t)rc * states * states;
            fp_t* Pmid = base_mid ? base_mid + (size_t)rc * states * states : nullptr;
            fp_t* Pprox = base_mid_prox ? (base_mid_prox + (size_t)rc * states * states) : nullptr;
            fp_t* Pdist = base_mid_dist ? (base_mid_dist + (size_t)rc * states * states) : nullptr;
            pmatrix_from_triple(
                er.Vinv.data(), er.V.data(), er.lambdas.data(),
                            r, t, p, pbuf.data(), states);
            if (include_midpoint_pmats) {
                pmatrix_from_triple(
                    er.Vinv.data(), er.V.data(), er.lambdas.data(),
                                r, t * 0.5, p, pbuf_mid.data(), states);
            }
            for (size_t idx = 0; idx < pbuf.size(); ++idx) {
                P[idx] = static_cast<fp_t>(pbuf[idx]);
                if (Pmid) Pmid[idx] = static_cast<fp_t>(pbuf_mid[idx]);
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

DeviceTree make_query_view(const DeviceTree& D, int query_idx)
{
    // This is a non-owning slice: all non-query pointers still alias D, while
    // query pointers are rebased so kernels can address the selected row as 0.
    DeviceTree view = D;
    if (query_idx < 0 || query_idx >= D.placement_queries) return view;
    const size_t clv_span = D.per_node_elems();
    view.query_capacity = 1;
    view.placement_queries = 1;
    if (D.d_query_chars) {
        view.d_query_chars = D.d_query_chars + static_cast<size_t>(query_idx) * D.sites;
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

static void build_query_clv(
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

__global__ void SeedRootDownClvKernel(DeviceTree D, int root_id)
{
    fp_t* root_down = down_clv_ptr<fp_t>(D, root_id);
    if (!root_down) return;
    const size_t per_node = per_node_span(D);
    const size_t idx = static_cast<size_t>(blockIdx.x) * static_cast<size_t>(blockDim.x) +
        static_cast<size_t>(threadIdx.x);
    if (idx >= per_node) return;
    root_down[idx] = fp_t(1);
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
        copy_device_to_device_async(dst.d_clv_up, src.d_clv_up, sizeof(fp_t) * clv_elems, stream);
    }

    const size_t scaler_elems = static_cast<size_t>(src.N) * src.scaler_elems();
    if (scaler_elems > 0 && src.d_site_scaler_up && dst.d_site_scaler_up) {
        copy_device_to_device_async(
            dst.d_site_scaler_up,
            src.d_site_scaler_up,
            sizeof(unsigned) * scaler_elems,
            stream);
    }
}

__global__ void CopySelectedUpwardStateKernel(
    DeviceTree src,
    DeviceTree dst,
    const int* node_ids,
    int node_count)
{
    const int selected_idx = static_cast<int>(blockIdx.y);
    if (selected_idx >= node_count) return;
    const int node_id = node_ids[selected_idx];
    if (node_id < 0 || node_id >= src.N) return;
    const size_t node_elems = src.sites * static_cast<size_t>(src.rate_cats) *
        static_cast<size_t>(src.states);
    for (size_t elem = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
         elem < node_elems;
         elem += static_cast<size_t>(blockDim.x) * gridDim.x) {
        const size_t offset = static_cast<size_t>(node_id) * node_elems + elem;
        dst.d_clv_up[offset] = src.d_clv_up[offset];
    }
    const size_t scaler_elems = src.per_rate_scaling
        ? src.sites * static_cast<size_t>(src.rate_cats)
        : src.sites;
    for (size_t elem = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
         elem < scaler_elems;
         elem += static_cast<size_t>(blockDim.x) * gridDim.x) {
        const size_t offset = static_cast<size_t>(node_id) * scaler_elems + elem;
        dst.d_site_scaler_up[offset] = src.d_site_scaler_up[offset];
    }
}

void copy_selected_upward_state(
    const DeviceTree& src,
    DeviceTree& dst,
    const int* d_node_ids,
    int node_count,
    cudaStream_t stream)
{
    if (node_count <= 0) return;
    if (src.N != dst.N || src.per_node_elems() != dst.per_node_elems() ||
        src.scaler_elems() != dst.scaler_elems()) {
        throw std::runtime_error("copy_selected_upward_state: tree shape mismatch.");
    }
    const unsigned blocks_x = static_cast<unsigned>(std::min<size_t>(
        1024, (src.per_node_elems() + 255) / 256));
    const dim3 grid(std::max(1u, blocks_x), static_cast<unsigned>(node_count));
    CopySelectedUpwardStateKernel<<<grid, 256, 0, stream>>>(
        src, dst, d_node_ids, node_count);
    CUDA_CHECK(cudaGetLastError());
}

// Populate compact NNI scoring targets from resident full-tree state. Each
// block-row handles one mapping and each thread handles one site. The target
// upward CLV is copied directly; the edge-outside CLV either takes the cached
// fast path or is reconstructed from parent-down and sibling-up messages while
// preserving the source scaler convention.
__global__ void BuildDirectNNITargetContextsKernel(
    DeviceTree src,
    DeviceTree dst,
    const DirectNNIContextOp* ops,
    int op_count)
{
    const int op_idx = static_cast<int>(blockIdx.y);
    const size_t site = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (op_idx >= op_count || site >= src.sites) return;
    const DirectNNIContextOp op = ops[op_idx];
    if (op.target_src < 0 || op.target_src >= src.N ||
        op.dst_target < 0 || op.dst_target >= dst.N) return;

    const size_t per_site = static_cast<size_t>(src.rate_cats) * src.states;
    const size_t src_node_span = src.sites * per_site;
    const size_t dst_node_span = dst.sites * per_site;
    const size_t scaler_span = src.per_rate_scaling
        ? src.sites * static_cast<size_t>(src.rate_cats) : src.sites;
    const size_t src_target_base = static_cast<size_t>(op.target_src) * src_node_span + site * per_site;
    const size_t dst_target_base = static_cast<size_t>(op.dst_target) * dst_node_span + site * per_site;
    const size_t src_target_scaler = static_cast<size_t>(op.target_src) * scaler_span;
    const size_t dst_target_scaler = static_cast<size_t>(op.dst_target) * scaler_span;

    for (size_t elem = 0; elem < per_site; ++elem)
        dst.d_clv_up[dst_target_base + elem] = src.d_clv_up[src_target_base + elem];
    if (src.per_rate_scaling) {
        for (int rate = 0; rate < src.rate_cats; ++rate)
            dst.d_site_scaler_up[dst_target_scaler + site * src.rate_cats + rate] =
                src.d_site_scaler_up[src_target_scaler + site * src.rate_cats + rate];
    } else {
        dst.d_site_scaler_up[dst_target_scaler + site] =
            src.d_site_scaler_up[src_target_scaler + site];
    }

    const size_t dst_edge_outside = static_cast<size_t>(op.dst_target) * dst_node_span + site * per_site;
    const size_t dst_mid_scaler = static_cast<size_t>(op.dst_target) * scaler_span;
    if (op.direct_edge_outside_src >= 0) {
        const size_t source_base = static_cast<size_t>(op.direct_edge_outside_src) * src_node_span + site * per_site;
        const size_t source_scaler = static_cast<size_t>(op.direct_edge_outside_src) * scaler_span;
        for (size_t elem = 0; elem < per_site; ++elem)
            dst.d_edge_outside_clv[dst_edge_outside + elem] = src.d_edge_outside_clv[source_base + elem];
        if (src.per_rate_scaling) {
            for (int rate = 0; rate < src.rate_cats; ++rate)
                dst.d_edge_outside_scaler[dst_mid_scaler + site * src.rate_cats + rate] =
                    src.d_edge_outside_scaler[source_scaler + site * src.rate_cats + rate];
        } else {
            dst.d_edge_outside_scaler[dst_mid_scaler + site] =
                src.d_edge_outside_scaler[source_scaler + site];
        }
        return;
    }

    if (op.parent_down_src < 0 || op.sibling_up_src < 0 ||
        op.second_pmat_src < 0) return;
    const size_t parent_base = static_cast<size_t>(op.parent_down_src) * src_node_span + site * per_site;
    const size_t sibling_base = static_cast<size_t>(op.sibling_up_src) * src_node_span + site * per_site;
    const size_t parent_scaler = static_cast<size_t>(op.parent_down_src) * scaler_span;
    const size_t sibling_scaler = static_cast<size_t>(op.sibling_up_src) * scaler_span;
    const size_t matrix_span = static_cast<size_t>(src.states) * src.states;
    const size_t pmat_span = static_cast<size_t>(src.rate_cats) * matrix_span;
    const fp_t* first_pmats = src.d_pmat + static_cast<size_t>(op.sibling_up_src) * pmat_span;
    const fp_t* second_pmats = src.d_pmat + static_cast<size_t>(op.second_pmat_src) * pmat_span;

    for (int rate = 0; rate < src.rate_cats; ++rate) {
        fp_t maximum = fp_t(0);
        for (int state = 0; state < src.states; ++state) {
            fp_t transformed = fp_t(0);
            for (int middle = 0; middle < src.states; ++middle) {
                fp_t first = fp_t(0);
                for (int child = 0; child < src.states; ++child) {
                    first = fp_fma(
                        first_pmats[static_cast<size_t>(rate) * matrix_span +
                            static_cast<size_t>(middle) * src.states + child],
                        src.d_clv_up[sibling_base + static_cast<size_t>(rate) * src.states + child],
                        first);
                }
                transformed = fp_fma(
                    second_pmats[static_cast<size_t>(rate) * matrix_span +
                        static_cast<size_t>(state) * src.states + middle],
                    first, transformed);
            }
            const size_t parent_index =
                parent_base + static_cast<size_t>(rate) * src.states + state;
            const fp_t value = src.d_clv_down[parent_index] * transformed;
            dst.d_edge_outside_clv[dst_edge_outside + static_cast<size_t>(rate) * src.states + state] = value;
            maximum = fp_fmax(maximum, value);
        }
        const fp_t scale_threshold = fp_ldexp(fp_t(1), kClvScaleThresholdExponent);
        const unsigned local_shift = maximum < scale_threshold
            ? static_cast<unsigned>(-kClvScaleThresholdExponent) : 0u;
        const unsigned inherited = src.per_rate_scaling
            ? src.d_site_scaler_down[parent_scaler + site * src.rate_cats + rate] +
              src.d_site_scaler_up[sibling_scaler + site * src.rate_cats + rate]
            : src.d_site_scaler_down[parent_scaler + site] +
              src.d_site_scaler_up[sibling_scaler + site];
        if (src.per_rate_scaling)
            dst.d_edge_outside_scaler[dst_mid_scaler + site * src.rate_cats + rate] = inherited + local_shift;
        else if (rate == 0)
            dst.d_edge_outside_scaler[dst_mid_scaler + site] = inherited + local_shift;
        if (local_shift) {
            for (int state = 0; state < src.states; ++state)
                fp_scale_pow2(
                    dst.d_edge_outside_clv[dst_edge_outside + static_cast<size_t>(rate) * src.states + state],
                    local_shift);
        }
    }
}

void build_direct_nni_target_contexts(
    const DeviceTree& src,
    DeviceTree& dst,
    const std::vector<DirectNNIContextOp>& ops,
    cudaStream_t stream)
{
    if (ops.empty()) return;
    if (ops.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
        ops.size() > static_cast<size_t>(std::numeric_limits<unsigned>::max())) {
        throw std::length_error(
            "direct NNI operation count exceeds CUDA indexing limits");
    }
    DirectNNIContextOp* d_ops = nullptr;
    const size_t ops_bytes = mlipper::util::checked_allocation_bytes<DirectNNIContextOp>(
        ops.size(), "direct NNI operations");
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_ops), ops_bytes));
    CUDA_CHECK(cudaMemcpyAsync(
        d_ops, ops.data(), ops_bytes, cudaMemcpyHostToDevice, stream));
    const dim3 grid(static_cast<unsigned>((src.sites + 255) / 256), static_cast<unsigned>(ops.size()));
    BuildDirectNNITargetContextsKernel<<<grid, 256, 0, stream>>>(src, dst, d_ops, static_cast<int>(ops.size()));
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaFreeAsync(d_ops, stream));
}

static void populate_host_topology(
    const TreeBuildResult& tree,
    HostPacking& host)
{
    const int node_count = static_cast<int>(tree.nodes.size());
    host.postorder = tree.postorder;
    host.preorder = tree.preorder;
    host.parent.assign(node_count, -1);
    host.left.assign(node_count, -1);
    host.right.assign(node_count, -1);
    host.is_tip.assign(node_count, 0);
    host.blen.assign(node_count, fp_t(0));
    for (int node_id = 0; node_id < node_count; ++node_id) {
        const TreeNode& node = tree.nodes[static_cast<size_t>(node_id)];
        host.parent[static_cast<size_t>(node_id)] = node.parent;
        host.left[static_cast<size_t>(node_id)] = node.left;
        host.right[static_cast<size_t>(node_id)] = node.right;
        host.is_tip[static_cast<size_t>(node_id)] = node.is_tip ? 1 : 0;
        host.blen[static_cast<size_t>(node_id)] =
            node.branch_length_to_parent;
    }
}

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
    for (const std::string& row : msa_rows) {
        if (row.size() != sites) {
            throw std::runtime_error("MSA row length does not match sites.");
        }
    }

    HostPacking H;
    populate_host_topology(T, H);

    std::unordered_map<std::string,int> name_to_row;
    name_to_row.reserve(msa_tip_names.size()*2);
    for (int row_idx = 0; row_idx < (int)msa_tip_names.size(); ++row_idx) {
        const auto [_, inserted] = name_to_row.emplace(msa_tip_names[row_idx], row_idx);
        if (!inserted) {
            throw std::runtime_error("Duplicate MSA tip name: " + msa_tip_names[row_idx]);
        }
    }

    std::vector<int> tip_node_ids_host;
    tip_node_ids_host.reserve(msa_tip_names.size());
    for (int id : T.postorder) {
        if (T.nodes[id].is_tip) tip_node_ids_host.push_back(id);
    }
    const int tips = (int)tip_node_ids_host.size();

    if (tips != (int)msa_tip_names.size())
        throw std::runtime_error("Tip count in tree != MSA names.");

    H.tip_node_ids = tip_node_ids_host;

    H.tipchars.resize(mlipper::util::checked_mul_size(
        static_cast<size_t>(tips), sites, "host tip character matrix"));

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

OwnedDeviceTree::~OwnedDeviceTree() noexcept
{
    reset();
}

OwnedDeviceTree::OwnedDeviceTree(OwnedDeviceTree&& other) noexcept
{
    static_cast<DeviceTree&>(*this) = static_cast<DeviceTree&>(other);
    static_cast<DeviceTree&>(other) = DeviceTree{};
}

OwnedDeviceTree& OwnedDeviceTree::operator=(OwnedDeviceTree&& other) noexcept
{
    if (this != &other) {
        reset();
        static_cast<DeviceTree&>(*this) = static_cast<DeviceTree&>(other);
        static_cast<DeviceTree&>(other) = DeviceTree{};
    }
    return *this;
}

void OwnedDeviceTree::reset() noexcept
{
    release_device_tree_buffers(*this);
}

// Materialize the complete resident likelihood state in four phases: validate
// host dimensions, compute capacities/strides, allocate every owned buffer,
// then upload immutable model/topology/alignment data. Construction uses a
// temporary owner so a partial CUDA failure cannot leak or publish half-built
// state. CLV and scaler fields inside DeviceTree are borrowed slices of the
// larger allocations released by OwnedDeviceTree::reset().
void allocate_device_tree_on_current_gpu(
    OwnedDeviceTree& target,
    const TreeBuildResult& T,
    const HostPacking& H,
    const EigResult& er,
    const std::vector<double>& rate_weights,
    const std::vector<double>& rate_multipliers,
    const std::vector<double>& pi,
    size_t sites, int states, int rate_cats, bool per_rate_scaling,
    const PlacementQueryBatch* queries,
    bool commit_to_tree,
    int insert_capacity,
    bool allocate_directional_clvs)
{
    if (target.device_id >= 0 || target.d_tipchars != nullptr) {
        throw std::runtime_error(
            "allocate_device_tree_on_current_gpu requires an empty DeviceTree");
    }
    // Build into a temporary owner. A failed allocation/copy releases only the
    // candidate; target is changed atomically by the final move assignment.
    OwnedDeviceTree candidate;
    DeviceTree& D = candidate;
    if (T.nodes.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
        H.tip_node_ids.size() >
            static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("tree size exceeds integer indexing limits");
    }
    const int node_count = static_cast<int>(T.nodes.size());
    const int tip_count = static_cast<int>(H.tip_node_ids.size());
    if (queries &&
        queries->size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("placement query count exceeds integer indexing limits");
    }
    const int query_count = queries
        ? static_cast<int>(queries->size())
        : 0;
    const size_t site_count = sites;
    const size_t state_count = static_cast<size_t>(states);
    const size_t rate_count = static_cast<size_t>(rate_cats);
    const size_t matrix_elems = mlipper::util::checked_mul_size(
        state_count, state_count, "model transition matrix");

    if (node_count <= 0 || tip_count <= 0 || sites == 0 ||
        states <= 0 || rate_cats <= 0 || T.root_id < 0 ||
        T.root_id >= node_count) {
        throw std::runtime_error("invalid DeviceTree shape");
    }
    if (H.blen.size() != static_cast<size_t>(node_count) ||
        H.tip_node_ids.size() != static_cast<size_t>(tip_count) ||
        H.tipchars.size() != mlipper::util::checked_mul_size(
            static_cast<size_t>(tip_count), sites, "host tip character matrix")) {
        throw std::runtime_error("host tree buffers do not match DeviceTree shape");
    }
    if (queries && queries->query_chars.size() !=
            mlipper::util::checked_mul_size(
                queries->size(), sites, "placement query character matrix")) {
        throw std::runtime_error(
            "placement query character matrix does not match query count and sites");
    }
    if (static_cast<int>(rate_weights.size()) != rate_cats ||
        static_cast<int>(rate_multipliers.size()) != rate_cats ||
        static_cast<int>(pi.size()) != states ||
        static_cast<int>(er.lambdas.size()) != states ||
        er.V.size() != matrix_elems || er.Vinv.size() != matrix_elems) {
        throw std::runtime_error("model buffers do not match DeviceTree shape");
    }

    D.N = node_count;
    D.device_id = mlipper::gpu::current_device_or_throw();
    D.tips = tip_count;
    D.inners = D.N - D.tips;
    D.placement_queries = query_count;
    D.root_id = T.root_id;

    if (insert_capacity < 0) {
        throw std::runtime_error("insert_capacity must be non-negative");
    }
    int reserve_inserts = 0;
    if (commit_to_tree) {
        reserve_inserts = std::max(insert_capacity, D.placement_queries);
    }

    if (reserve_inserts > (std::numeric_limits<int>::max() - D.N) / 2 ||
        reserve_inserts > std::numeric_limits<int>::max() - D.tips) {
        throw std::length_error("DeviceTree insertion capacity exceeds integer indexing limits");
    }
    D.capacity_N = D.N + 2 * reserve_inserts;
    D.capacity_tips = D.tips + reserve_inserts;
    D.query_capacity = D.placement_queries;
    D.sites = sites;
    D.states = states;
    D.rate_cats = rate_cats;
    D.log2_stride = ceil_log2_u32((unsigned int)(D.states + 1));
    D.per_rate_scaling = per_rate_scaling;

    // Rate multipliers are folded into the eigenvalues once so kernels can use
    // branch length directly without another per-site multiplication.
    std::vector<fp_t> lambdas_scaled(mlipper::util::checked_mul_size(
        rate_count, state_count, "scaled model eigenvalues"), fp_t(0));
    for (int rate_idx = 0; rate_idx < D.rate_cats; ++rate_idx) {
        const double rate_multiplier = rate_multipliers[rate_idx];
        for (int state_idx = 0; state_idx < D.states; ++state_idx) {
            lambdas_scaled[static_cast<size_t>(rate_idx) * state_count + static_cast<size_t>(state_idx)] =
                static_cast<fp_t>(er.lambdas[state_idx] * rate_multiplier);
        }
    }

    const std::vector<fp_t> V_fp(er.V.begin(), er.V.end());
    const std::vector<fp_t> Vinv_fp(er.Vinv.begin(), er.Vinv.end());
    const std::vector<fp_t> rate_weights_fp(rate_weights.begin(), rate_weights.end());
    const std::vector<fp_t> pi_fp(pi.begin(), pi.end());

    const size_t model_matrix_bytes = mlipper::util::checked_allocation_bytes<fp_t>(
        matrix_elems, "model transition matrix");
    const size_t lambdas_bytes = mlipper::util::checked_allocation_bytes<fp_t>(
        lambdas_scaled.size(), "scaled model eigenvalues");
    const size_t rate_weights_bytes = mlipper::util::checked_allocation_bytes<fp_t>(
        rate_count, "model rate weights");
    const size_t frequencies_bytes = mlipper::util::checked_allocation_bytes<fp_t>(
        state_count, "model frequencies");

    cuda_malloc_bytes(D.d_lambdas, lambdas_bytes);
    cuda_malloc_bytes(D.d_V, model_matrix_bytes);
    cuda_malloc_bytes(D.d_Vinv, model_matrix_bytes);
    cuda_malloc_bytes(D.d_rate_weights, rate_weights_bytes);
    cuda_malloc_bytes(D.d_frequencies, frequencies_bytes);

    copy_host_to_device(D.d_lambdas, lambdas_scaled.data(), lambdas_bytes);
    copy_host_to_device(D.d_V, V_fp.data(), model_matrix_bytes);
    copy_host_to_device(D.d_Vinv, Vinv_fp.data(), model_matrix_bytes);
    copy_host_to_device(D.d_rate_weights, rate_weights_fp.data(), rate_weights_bytes);
    copy_host_to_device(D.d_frequencies, pi_fp.data(), frequencies_bytes);

    const size_t capacity_node_bytes = mlipper::util::checked_allocation_bytes<fp_t>(
        static_cast<size_t>(D.capacity_N), "tree node buffers");
    const size_t live_node_bytes = mlipper::util::checked_allocation_bytes<fp_t>(
        static_cast<size_t>(D.N), "live tree nodes");
    const size_t tipchar_capacity_bytes = mlipper::util::checked_product(
        "tip character buffer", sizeof(uint8_t),
        static_cast<size_t>(D.capacity_tips), site_count);
    const size_t tipchar_live_bytes = mlipper::util::checked_product(
        "live tip characters", sizeof(uint8_t), static_cast<size_t>(D.tips), site_count);
    const size_t tip_node_ids_capacity_bytes = mlipper::util::checked_allocation_bytes<int>(
        static_cast<size_t>(D.capacity_tips), "tip node id buffer");
    const size_t tip_node_ids_live_bytes = mlipper::util::checked_allocation_bytes<int>(
        static_cast<size_t>(D.tips), "live tip node ids");

    cuda_malloc_zeroed(D.d_blen, capacity_node_bytes);
    cuda_malloc_zeroed(D.d_new_pendant_length, capacity_node_bytes);
    cuda_malloc_zeroed(D.d_new_proximal_length, capacity_node_bytes);
    cuda_malloc_zeroed(D.d_prev_pendant_length, capacity_node_bytes);
    cuda_malloc_zeroed(D.d_prev_proximal_length, capacity_node_bytes);
    copy_host_to_device(D.d_blen, H.blen.data(), live_node_bytes);

    cuda_malloc_bytes(D.d_tipchars, tipchar_capacity_bytes);
    copy_host_to_device(D.d_tipchars, H.tipchars.data(), tipchar_live_bytes);
    cuda_malloc_bytes(D.d_tip_node_ids, tip_node_ids_capacity_bytes);
    copy_host_to_device(D.d_tip_node_ids, H.tip_node_ids.data(), tip_node_ids_live_bytes);

    const size_t clv_capacity_elems = clv_pool_elements(D);
    const size_t per_node = D.per_node_elems();
    const size_t clv_total = mlipper::util::checked_mul_size(
        clv_capacity_elems, 2, "upward/downward CLV pools");
    const size_t clv_total_bytes = mlipper::util::checked_allocation_bytes<fp_t>(
        clv_total, "upward/downward CLV pools");
    const size_t clv_capacity_bytes = mlipper::util::checked_allocation_bytes<fp_t>(
        clv_capacity_elems, "directional CLV pool");

    if (allocate_directional_clvs) {
        cuda_malloc_zeroed(D.d_clv_up, clv_total_bytes);
        cuda_malloc_zeroed(D.d_edge_midpoint_clv, clv_capacity_bytes);
        cuda_malloc_zeroed(D.d_edge_outside_clv, clv_capacity_bytes);
        const size_t scaler_pool = D.scaler_pool_elems();
        cuda_malloc_zeroed(
            D.d_scaler_storage,
            mlipper::util::checked_product(
                "directional scaler storage", sizeof(unsigned), scaler_pool, 4));
        configure_resident_buffer_views(D);
    }

    if (allocate_directional_clvs && per_node > 0) {
        std::vector<fp_t> ones(per_node, fp_t(1));
        copy_host_to_device(
            D.d_clv_down + static_cast<size_t>(T.root_id) * per_node,
            ones.data(),
            mlipper::util::checked_allocation_bytes<fp_t>(
                per_node, "root downward CLV"));
    }

    const size_t pmat_elems_cur = mlipper::util::checked_product(
        "live transition matrices", static_cast<size_t>(D.N), rate_count, matrix_elems);
    const size_t pmat_elems_cap = mlipper::util::checked_product(
        "transition matrix capacity", static_cast<size_t>(D.capacity_N),
        rate_count, matrix_elems);
    const size_t pmat_bytes_cur = mlipper::util::checked_allocation_bytes<fp_t>(
        pmat_elems_cur, "live transition matrices");
    const size_t pmat_bytes_cap = mlipper::util::checked_allocation_bytes<fp_t>(
        pmat_elems_cap, "transition matrix capacity");
    if (H.pmats.size() != pmat_elems_cur ||
        H.pmats_mid.size() != pmat_elems_cur) {
        throw std::runtime_error(
            "host transition matrices do not match the live DeviceTree shape");
    }
    cuda_malloc_zeroed(D.d_pmat, pmat_bytes_cap);
    copy_host_to_device(D.d_pmat, H.pmats.data(), pmat_bytes_cur);
    cuda_malloc_zeroed(D.d_pmat_mid, pmat_bytes_cap);
    copy_host_to_device(D.d_pmat_mid, H.pmats_mid.data(), pmat_bytes_cur);
    cuda_malloc_zeroed(D.d_pmat_mid_prox, pmat_bytes_cap);
    if (!H.pmats_mid_prox.empty() && H.pmats_mid_prox.size() != pmat_elems_cur) {
        throw std::runtime_error("pmats_mid_prox size mismatch.");
    }
    const fp_t* pmat_mid_prox_src = H.pmats_mid_prox.empty() ? H.pmats_mid.data() : H.pmats_mid_prox.data();
    copy_host_to_device(D.d_pmat_mid_prox, pmat_mid_prox_src, pmat_bytes_cur);
    cuda_malloc_zeroed(D.d_pmat_mid_dist, pmat_bytes_cap);
    if (!H.pmats_mid_dist.empty() && H.pmats_mid_dist.size() != pmat_elems_cur) {
        throw std::runtime_error("pmats_mid_dist size mismatch.");
    }
    const fp_t* pmat_mid_dist_src = H.pmats_mid_dist.empty() ? H.pmats_mid.data() : H.pmats_mid_dist.data();
    copy_host_to_device(D.d_pmat_mid_dist, pmat_mid_dist_src, pmat_bytes_cur);

    if (!H.pattern_weights.empty()) {
        if (H.pattern_weights.size() != sites) {
            throw std::runtime_error("pattern_weights size mismatch.");
        }
        cuda_malloc_bytes(D.d_pattern_weights_u,
            mlipper::util::checked_allocation_bytes<unsigned>(
                sites, "pattern weights"));
        copy_host_to_device(D.d_pattern_weights_u, H.pattern_weights.data(), sizeof(unsigned) * sites);
    }

    if (queries && !queries->empty()) {
        const size_t qcount = queries->size();
        allocate_query_buffers(
            D,
            D.query_capacity,
            D.d_query_chars,
            D.d_query_clv);
        copy_host_to_device(
            D.d_query_chars,
            queries->query_chars.data(),
            mlipper::util::checked_product(
                "live query characters", sizeof(uint8_t), qcount,
                site_count));
    }

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
    cuda_malloc_bytes(D.d_tipmap,
        mlipper::util::checked_allocation_bytes<unsigned>(
            tipmap_size, "tip decode map"));
    copy_host_to_device(D.d_tipmap, tipmap.data(), tipmap_size * sizeof(unsigned int));
    validate_resident_buffer_views(
        D,
        "allocate_device_tree_on_current_gpu");
    target = std::move(candidate);
}

namespace {

enum ReloadSkipFlags : unsigned {
    RELOAD_SKIP_ZERO_LENGTH_SCRATCH = 1u << 0,
    RELOAD_SKIP_COPY_TIPCHARS = 1u << 1,
    RELOAD_SKIP_CLEAR_LIKELIHOOD_BUFFERS = 1u << 2,
    RELOAD_SKIP_COPY_PMATS = 1u << 3,
    RELOAD_SKIP_COPY_PATTERN_WEIGHTS = 1u << 4,
};

void reload_device_tree_live_data_impl(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const PlacementQueryBatch* queries,
    unsigned skip_flags,
    cudaStream_t stream)
{
    ensure_device_tree_current_device(D, "reload_device_tree_live_data_impl");
    validate_resident_buffer_views(
        D,
        "reload_device_tree_live_data_impl");
    if (T.nodes.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
        H.tip_node_ids.size() >
            static_cast<size_t>(std::numeric_limits<int>::max()) ||
        (queries && queries->size() >
            static_cast<size_t>(std::numeric_limits<int>::max()))) {
        throw std::length_error(
            "reload_device_tree_live_data: live counts exceed integer limits");
    }
    const int node_count = static_cast<int>(T.nodes.size());
    const int tip_count = static_cast<int>(H.tip_node_ids.size());
    const int query_count = queries
        ? static_cast<int>(queries->size())
        : D.placement_queries;
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
    if (queries && queries->query_chars.size() !=
            mlipper::util::checked_mul_size(
                queries->size(), D.sites,
                "reload query character matrix")) {
        throw std::runtime_error(
            "reload_device_tree_live_data: query character size mismatch");
    }

    D.N = node_count;
    D.tips = tip_count;
    D.inners = D.N - D.tips;
    D.root_id = T.root_id;
    D.placement_queries = query_count;

    const size_t per_node = D.per_node_elems();
    const size_t matrix_per_node = D.pmat_per_node_elems();
    const size_t scaler_pool_bytes =
        mlipper::util::checked_allocation_bytes<unsigned>(
            D.scaler_storage_elems(), "reload scaler storage");
    const size_t capacity_node_bytes =
        mlipper::util::checked_allocation_bytes<fp_t>(
            static_cast<size_t>(D.capacity_N), "reload node capacity");
    const size_t live_node_bytes =
        mlipper::util::checked_allocation_bytes<fp_t>(
            static_cast<size_t>(D.N), "reload live nodes");
    const size_t live_tipchar_bytes = mlipper::util::checked_product(
        "reload live tip characters", sizeof(uint8_t),
        static_cast<size_t>(D.tips), D.sites);
    const size_t live_tip_node_ids_bytes =
        mlipper::util::checked_allocation_bytes<int>(
            static_cast<size_t>(D.tips), "reload live tip IDs");
    const size_t clv_capacity_bytes =
        mlipper::util::checked_allocation_bytes<fp_t>(
            mlipper::util::checked_mul_size(
                static_cast<size_t>(D.capacity_N), per_node,
                "reload CLV capacity"),
            "reload CLV capacity");
    const size_t pmat_live_bytes =
        mlipper::util::checked_allocation_bytes<fp_t>(
            mlipper::util::checked_mul_size(
                static_cast<size_t>(D.N), matrix_per_node,
                "reload live transition matrices"),
            "reload live transition matrices");

    copy_host_to_device_async(D.d_blen, H.blen.data(), live_node_bytes, stream);

    if ((skip_flags & RELOAD_SKIP_ZERO_LENGTH_SCRATCH) == 0u) {
        zero_length_scratch_async(D, capacity_node_bytes, stream);
    }

    if ((skip_flags & RELOAD_SKIP_COPY_TIPCHARS) == 0u &&
        D.tips > 0 && D.d_tipchars) {
        copy_host_to_device_async(D.d_tipchars, H.tipchars.data(), live_tipchar_bytes, stream);
        if (D.d_tip_node_ids) {
            copy_host_to_device_async(
                D.d_tip_node_ids,
                H.tip_node_ids.data(),
                live_tip_node_ids_bytes,
                stream);
        }
    }

    if ((skip_flags & RELOAD_SKIP_CLEAR_LIKELIHOOD_BUFFERS) == 0u) {
        if (D.d_clv_up) {
            CUDA_CHECK(cudaMemsetAsync(D.d_clv_up, 0, clv_capacity_bytes * 2, stream));
        }
        if (D.d_edge_midpoint_clv) {
            CUDA_CHECK(cudaMemsetAsync(D.d_edge_midpoint_clv, 0, clv_capacity_bytes, stream));
        }
        if (D.d_edge_outside_clv) {
            CUDA_CHECK(cudaMemsetAsync(D.d_edge_outside_clv, 0, clv_capacity_bytes, stream));
        }
        if (D.d_scaler_storage) {
            CUDA_CHECK(cudaMemsetAsync(D.d_scaler_storage, 0, scaler_pool_bytes, stream));
        }
    }

    if (per_node > 0) {
        dim3 block(256);
        dim3 grid(static_cast<unsigned>((per_node + block.x - 1) / block.x));
        SeedRootDownClvKernel<<<grid, block, 0, stream>>>(D, T.root_id);
        CUDA_CHECK(cudaGetLastError());
    }

    const size_t expected_pmat_elements = mlipper::util::checked_mul_size(
        static_cast<size_t>(D.N), matrix_per_node,
        "reload transition matrix elements");
    if (H.pmats.size() != expected_pmat_elements ||
        H.pmats_mid.size() != expected_pmat_elements) {
        throw std::runtime_error("reload_device_tree_live_data: PMAT host size mismatch.");
    }
    if ((!H.pmats_mid_prox.empty() &&
         H.pmats_mid_prox.size() != expected_pmat_elements) ||
        (!H.pmats_mid_dist.empty() &&
         H.pmats_mid_dist.size() != expected_pmat_elements)) {
        throw std::runtime_error(
            "reload_device_tree_live_data: split PMAT host size mismatch");
    }
    if ((skip_flags & RELOAD_SKIP_COPY_PMATS) == 0u) {
        const fp_t* pmat_mid_prox_src =
            H.pmats_mid_prox.empty() ? H.pmats_mid.data() : H.pmats_mid_prox.data();
        const fp_t* pmat_mid_dist_src =
            H.pmats_mid_dist.empty() ? H.pmats_mid.data() : H.pmats_mid_dist.data();
        copy_host_to_device_async(
            D.d_pmat,
            H.pmats.data(),
            pmat_live_bytes,
            stream);
        copy_host_to_device_async(
            D.d_pmat_mid,
            H.pmats_mid.data(),
            pmat_live_bytes,
            stream);
        if (D.d_pmat_mid_prox != nullptr) {
            copy_host_to_device_async(
                D.d_pmat_mid_prox,
                pmat_mid_prox_src,
                pmat_live_bytes,
                stream);
        }
        if (D.d_pmat_mid_dist != nullptr) {
            copy_host_to_device_async(
                D.d_pmat_mid_dist,
                pmat_mid_dist_src,
                pmat_live_bytes,
                stream);
        }
    }

    if (queries && query_count > 0) {
        if (!queries->query_chars.empty()) {
            const size_t query_chars_bytes = mlipper::util::checked_product(
                "reload query characters", sizeof(uint8_t),
                static_cast<size_t>(query_count), D.sites);
            copy_host_to_device_async(
                D.d_query_chars,
                queries->query_chars.data(),
                query_chars_bytes,
                stream);
        }
    }

    if ((skip_flags & RELOAD_SKIP_COPY_PATTERN_WEIGHTS) == 0u &&
        D.d_pattern_weights_u && !H.pattern_weights.empty()) {
        copy_host_to_device_async(
            D.d_pattern_weights_u,
            H.pattern_weights.data(),
            sizeof(unsigned) * D.sites,
            stream);
    }
}

} // namespace

void reload_device_tree_live_data(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const PlacementQueryBatch* queries,
    cudaStream_t stream)
{
    reload_device_tree_live_data_impl(D, T, H, queries, 0u, stream);
}

void reload_device_tree_live_data_local_spr(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const HostPacking& base_H,
    int current_main_pmat_node,
    int& previous_main_pmat_node,
    const PlacementQueryBatch* queries,
    cudaStream_t stream)
{
    const unsigned reload_skip_flags =
        RELOAD_SKIP_COPY_PMATS |
        RELOAD_SKIP_ZERO_LENGTH_SCRATCH |
        RELOAD_SKIP_COPY_TIPCHARS |
        RELOAD_SKIP_CLEAR_LIKELIHOOD_BUFFERS |
        RELOAD_SKIP_COPY_PATTERN_WEIGHTS;
    // Ranking restores the immutable upward state after this reload and
    // recomputes every downward/midpoint value it consumes. Tip characters
    // and pattern weights are unchanged across prune roots.
    reload_device_tree_live_data_impl(
        D,
        T,
        H,
        queries,
        reload_skip_flags,
        stream);

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
        copy_host_to_device_async(
            D.d_pmat + offset,
            src.data() + offset,
            sizeof(fp_t) * per_node,
            stream);
    };

    if (previous_main_pmat_node >= 0 && previous_main_pmat_node != current_main_pmat_node) {
        copy_main_pmat_slice(base_H.pmats, previous_main_pmat_node);
    }
    copy_main_pmat_slice(H.pmats, current_main_pmat_node);
    previous_main_pmat_node = current_main_pmat_node;
}

void reload_device_tree_live_data_preserving_clvs(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const PlacementQueryBatch* queries,
    cudaStream_t stream)
{
    reload_device_tree_live_data_impl(
        D,
        T,
        H,
        queries,
        RELOAD_SKIP_CLEAR_LIKELIHOOD_BUFFERS,
        stream);
}

bool subtree_workspace_requires_rebuild(
    const OwnedSubtreeWorkspace& workspace,
    const TreeBuildResult& tree,
    const HostPacking& host,
    size_t sites,
    int states,
    int rate_cats,
    bool per_rate_scaling)
{
    const DeviceTree& D = workspace.dev;
    if (tree.nodes.size() >
            static_cast<size_t>(std::numeric_limits<int>::max()) ||
        host.tip_node_ids.size() >
            static_cast<size_t>(std::numeric_limits<int>::max())) {
        return true;
    }
    const int required_nodes = static_cast<int>(tree.nodes.size());
    const int required_tips = static_cast<int>(host.tip_node_ids.size());
    return D.device_id < 0 ||
           D.capacity_N < required_nodes ||
           D.capacity_tips < required_tips ||
           D.sites != sites ||
           D.states != states ||
           D.rate_cats != rate_cats ||
           D.per_rate_scaling != per_rate_scaling;
}

void ensure_device_tree_query_capacity(
    DeviceTree& D,
    int required_queries,
    const char* context)
{
    if (required_queries < 0) {
        throw std::runtime_error(
            std::string(context) +
            ": required_queries must be >= 0.");
    }

    ensure_device_tree_current_device(D, context);
    validate_resident_buffer_views(D, context);

    if (required_queries == 0) {
        D.placement_queries = 0;
        return;
    }

    if (D.query_capacity >= required_queries &&
        D.d_query_chars != nullptr &&
        D.d_query_clv != nullptr) {
        D.placement_queries = required_queries;
        return;
    }

    uint8_t* new_query_chars = nullptr;
    fp_t* new_query_clv = nullptr;
    allocate_query_buffers(
        D,
        required_queries,
        new_query_chars,
        new_query_clv);

    cuda_free_if_allocated(D.d_query_chars);
    cuda_free_if_allocated(D.d_query_clv);
    D.d_query_chars = new_query_chars;
    D.d_query_clv = new_query_clv;

    D.query_capacity = required_queries;
    D.placement_queries = required_queries;
}

void release_subtree_workspace(OwnedSubtreeWorkspace& workspace) noexcept
{
    workspace.dev.reset();
}

void load_subtree_workspace(
    OwnedSubtreeWorkspace& workspace,
    const SubtreeWorkspaceLoadConfig& config,
    cudaStream_t stream,
    const char* context)
{
    if (config.queries != nullptr &&
        config.queries->size() >
            static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error(
            std::string(context) +
            ": query count exceeds integer indexing limits");
    }
    const int query_count =
        config.queries != nullptr
            ? static_cast<int>(config.queries->size())
            : 0;
    const int required_query_capacity = config.query_capacity;
    if (required_query_capacity < query_count) {
        throw std::runtime_error(
            std::string(context) +
            ": query_capacity is smaller than the provided query batch.");
    }
    if (config.insert_capacity < 0) {
        throw std::runtime_error(
            std::string(context) + ": insert_capacity must be non-negative.");
    }

    if (subtree_workspace_requires_rebuild(
            workspace,
            config.tree,
            config.host,
            config.sites,
            config.states,
            config.rate_cats,
            config.per_rate_scaling)) {
        release_subtree_workspace(workspace);
        allocate_device_tree_on_current_gpu(
            workspace.dev,
            config.tree,
            config.host,
            config.eig,
            config.rate_weights,
            config.rate_multipliers,
            config.pi,
            config.sites,
            config.states,
            config.rate_cats,
            config.per_rate_scaling,
            config.queries,
            config.commit_to_tree,
            config.insert_capacity);
    }

    ensure_device_tree_current_device(workspace.dev, context);
    ensure_device_tree_query_capacity(
        workspace.dev,
        required_query_capacity,
        context);
    reload_device_tree_live_data(
        workspace.dev,
        config.tree,
        config.host,
        config.queries,
        stream);
}

static void release_device_tree_buffers(DeviceTree& D) noexcept
{
    if (D.device_id >= 0) {
        // Cleanup must remain safe during stack unwinding. CUDA cleanup errors
        // cannot be surfaced from a destructor, so switch devices best-effort
        // and let each cudaFree below independently release its allocation.
        cudaSetDevice(D.device_id);
    }
    cuda_free_if_allocated(D.d_blen);
    cuda_free_if_allocated(D.d_new_pendant_length);
    cuda_free_if_allocated(D.d_new_proximal_length);
    cuda_free_if_allocated(D.d_prev_pendant_length);
    cuda_free_if_allocated(D.d_prev_proximal_length);
    cuda_free_if_allocated(D.d_tipchars);
    cuda_free_if_allocated(D.d_tip_node_ids);
    cuda_free_if_allocated(D.d_clv_up);
    cuda_free_if_allocated(D.d_edge_midpoint_clv);
    cuda_free_if_allocated(D.d_edge_outside_clv);
    cuda_free_if_allocated(D.d_scaler_storage);
    cuda_free_if_allocated(D.d_lambdas);
    cuda_free_if_allocated(D.d_V);
    cuda_free_if_allocated(D.d_Vinv);
    cuda_free_if_allocated(D.d_rate_weights);
    cuda_free_if_allocated(D.d_frequencies);
    cuda_free_if_allocated(D.d_pmat);
    cuda_free_if_allocated(D.d_pmat_mid);
    cuda_free_if_allocated(D.d_pmat_mid_prox);
    cuda_free_if_allocated(D.d_pmat_mid_dist);
    cuda_free_if_allocated(D.d_pattern_weights_u);
    cuda_free_if_allocated(D.d_query_clv);
    cuda_free_if_allocated(D.d_query_chars);
    cuda_free_if_allocated(D.d_query_pmat);
    cuda_free_if_allocated(D.d_tipmap);
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

static void build_upward_ops_host_levelized(
    const TreeBuildResult& tree,
    const std::vector<int>& node_to_tip,
    std::vector<NodeOpInfo>& upward_ops_host,
    std::vector<int>& level_offsets)
{
    const int node_count = static_cast<int>(tree.nodes.size());
    upward_ops_host.clear();
    level_offsets.clear();
    if (node_count <= 0) return;

    std::vector<int> height(static_cast<size_t>(node_count), 0);
    int max_height = 0;
    for (int node_id : tree.postorder) {
        if (node_id < 0 || node_id >= node_count) continue;
        const TreeNode& node = tree.nodes[static_cast<size_t>(node_id)];
        if (!node.is_tip && node.left >= 0 && node.right >= 0) {
            height[static_cast<size_t>(node_id)] = 1 + std::max(
                height[static_cast<size_t>(node.left)],
                height[static_cast<size_t>(node.right)]);
            max_height = std::max(max_height, height[static_cast<size_t>(node_id)]);
        }
    }
    // Level zero contains only tips and therefore no operations. Each later
    // [offset[level], offset[level + 1]) range has no intra-level dependency.
    std::vector<int> counts(static_cast<size_t>(max_height) + 1, 0);
    for (int node_id : tree.postorder) {
        NodeOpInfo op{};
        if (make_upward_op_for_node(tree, node_to_tip, node_id, op))
            ++counts[static_cast<size_t>(height[static_cast<size_t>(node_id)])];
    }
    level_offsets.assign(counts.size() + 1, 0);
    for (size_t level = 0; level < counts.size(); ++level)
        level_offsets[level + 1] = level_offsets[level] + counts[level];
    upward_ops_host.resize(static_cast<size_t>(level_offsets.back()));
    std::vector<int> write_offsets(level_offsets.begin(), level_offsets.end() - 1);
    for (int node_id : tree.postorder) {
        NodeOpInfo op{};
        if (!make_upward_op_for_node(tree, node_to_tip, node_id, op)) continue;
        const int level = height[static_cast<size_t>(node_id)];
        upward_ops_host[static_cast<size_t>(write_offsets[static_cast<size_t>(level)]++)] = op;
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
    const int node_count = static_cast<int>(tree.nodes.size());
    downward_ops_host.clear();
    level_offsets.clear();
    if (node_count <= 0 || tree.root_id < 0 || tree.root_id >= node_count) {
        return;
    }

    // Depth zero is the seeded root. Both child operations at a later depth
    // depend only on the fully completed preceding depth.
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

    static thread_local LaunchConfig cfg_storage;
    const LaunchConfig& cfg = initialized_launch_config(
        cfg_storage,
        mlipper::likelihood::partials::UpdatePartialsUpwardKernel);
    const int grid = capped_site_grid(D.sites, cfg, num_sms);

    mlipper::likelihood::partials::UpdatePartialsUpwardKernel<<<grid, cfg.block, 0, stream>>>(
        D,
        d_ops,
        num_ops);
    CUDA_CHECK(cudaGetLastError());
}

static void launch_upward_clv_update_levelized(
    const DeviceTree& D,
    NodeOpInfo* d_ops,
    const std::vector<int>& level_offsets,
    int num_sms,
    cudaStream_t stream,
    bool initialize_tips = true)
{
    if (!d_ops || D.sites == 0 || level_offsets.size() < 2) return;
    if (initialize_tips) {
        launch_init_tip_clv(D, stream);
    }
    static thread_local LaunchConfig cfg_storage;
    const LaunchConfig& cfg = initialized_launch_config(
        cfg_storage, mlipper::likelihood::partials::UpdatePartialsUpwardLevelKernel);
    const unsigned grid_x = static_cast<unsigned>(
        std::max(1, capped_site_grid(D.sites, cfg, num_sms)));
    constexpr int kMaxGridY = 65535;
    for (size_t level = 1; level + 1 < level_offsets.size(); ++level) {
        const int start = level_offsets[level];
        const int count = level_offsets[level + 1] - start;
        for (int batch_start = 0; batch_start < count; batch_start += kMaxGridY) {
            const int batch_ops = std::min(count - batch_start, kMaxGridY);
            dim3 grid(grid_x, static_cast<unsigned>(batch_ops));
            mlipper::likelihood::partials::UpdatePartialsUpwardLevelKernel<<<grid, cfg.block, 0, stream>>>(
                D, d_ops + start + batch_start, batch_ops);
            CUDA_CHECK(cudaGetLastError());
        }
    }
}

static void launch_downward_clv_update(
    const DeviceTree& D,
    NodeOpInfo* d_ops,
    int num_ops,
    int num_sms,
    cudaStream_t stream)
{
    if (num_ops <= 0 || !d_ops) return;

    static thread_local LaunchConfig cfg_storage;
    const LaunchConfig& cfg = initialized_launch_config(
        cfg_storage,
        mlipper::likelihood::partials::UpdatePartialsDownwardKernel);
    const int grid = capped_site_grid(D.sites, cfg, num_sms);

    mlipper::likelihood::partials::UpdatePartialsDownwardKernel<<<grid, cfg.block, 0, stream>>>(D, d_ops, num_ops);
    CUDA_CHECK(cudaGetLastError());
}

static void launch_downward_clv_update_levelized(
    const DeviceTree& D,
    NodeOpInfo* d_ops,
    const std::vector<int>& level_offsets,
    int num_sms,
    cudaStream_t stream)
{
    if (!d_ops || D.sites == 0 || level_offsets.size() < 2) return;

    static thread_local LaunchConfig cfg_storage;
    const LaunchConfig& cfg = initialized_launch_config(
        cfg_storage,
        mlipper::likelihood::partials::UpdatePartialsDownwardLevelKernel);
    unsigned int grid_x = static_cast<unsigned int>(
        std::max(1, capped_site_grid(D.sites, cfg, num_sms)));

    constexpr int kMaxGridY = 65535;
    for (size_t level = 1; level + 1 < level_offsets.size(); ++level) {
        const int level_start = level_offsets[level];
        const int level_count = level_offsets[level + 1] - level_start;
        if (level_count <= 0) continue;

        for (int batch_start = 0; batch_start < level_count; batch_start += kMaxGridY) {
            const int batch_ops = std::min(level_count - batch_start, kMaxGridY);
            dim3 grid(grid_x, static_cast<unsigned int>(batch_ops));
            mlipper::likelihood::partials::UpdatePartialsDownwardLevelKernel<<<grid, cfg.block, 0, stream>>>(
                D,
                d_ops + level_start + batch_start,
                batch_ops);
            CUDA_CHECK(cudaGetLastError());
        }
    }
}

static int get_device_sm_count()
{
    // CUDA device selection is thread-local, so cache device properties at the
    // same scope to avoid cross-thread races in multi-GPU workflows.
    static thread_local int cached_device = -1;
    static thread_local int sm_count = 0;
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

    cuda_malloc_bytes(
        placement_ops.d_ops,
        mlipper::util::checked_allocation_bytes<NodeOpInfo>(
            static_cast<size_t>(required_ops), "placement operations"));
    placement_ops.capacity = required_ops;
}

static void upload_ops_to_device(
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    cudaStream_t stream)
{
    if (host_ops.size() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error(
            "placement operation count exceeds integer indexing limits");
    }
    const int num_ops = static_cast<int>(host_ops.size());
    if (num_ops <= 0) {
        placement_ops.num_ops = 0;
        return;
    }

    ensure_placement_op_capacity(placement_ops, num_ops, stream);
    copy_host_to_device_async(
        placement_ops.d_ops,
        host_ops.data(),
        mlipper::util::checked_allocation_bytes<NodeOpInfo>(
            static_cast<size_t>(num_ops), "uploaded placement operations"),
        stream);
    placement_ops.num_ops = num_ops;
}

static void run_upward_clv_stage(
    const DeviceTree& D,
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    int num_sms,
    cudaStream_t stream)
{
    upload_ops_to_device(placement_ops, host_ops, stream);
    launch_upward_clv_update(
        D, placement_ops.d_ops, placement_ops.num_ops, num_sms, stream);
}

static void run_upward_clv_stage_levelized(
    const DeviceTree& D,
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    const std::vector<int>& level_offsets,
    int num_sms,
    cudaStream_t stream)
{
    upload_ops_to_device(placement_ops, host_ops, stream);
    launch_upward_clv_update_levelized(
        D, placement_ops.d_ops, level_offsets, num_sms, stream);
}

static void run_downward_clv_stage(
    const DeviceTree& D,
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    int num_sms,
    cudaStream_t stream)
{
    upload_ops_to_device(placement_ops, host_ops, stream);
    launch_downward_clv_update(
        D, placement_ops.d_ops, placement_ops.num_ops, num_sms, stream);
}

static void run_downward_clv_stage_levelized(
    const DeviceTree& D,
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    const std::vector<int>& level_offsets,
    int num_sms,
    cudaStream_t stream)
{
    upload_ops_to_device(placement_ops, host_ops, stream);
    launch_downward_clv_update_levelized(
        D, placement_ops.d_ops, level_offsets, num_sms, stream);
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
    mlipper::likelihood::partials::InitializeTipPartialsKernel<<<grid, block, 0, stream>>>(D);
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

void UpdateTreeClvs(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    cudaStream_t stream)
{
    const int sm_count = get_device_sm_count();
    rebuild_node_to_tip_map(T, H, placement_ops.node_to_tip);
    build_upward_ops_host_levelized(
        T, placement_ops.node_to_tip, placement_ops.upward_ops_host,
        placement_ops.upward_level_offsets_host);
    run_upward_clv_stage_levelized(
        D,
        placement_ops,
        placement_ops.upward_ops_host,
        placement_ops.upward_level_offsets_host,
        sm_count,
        stream);

    build_downward_ops_host_levelized(
        T,
        placement_ops.node_to_tip,
        placement_ops.downward_ops_host,
        placement_ops.downward_level_offsets_host);
    run_downward_clv_stage_levelized(
        D,
        placement_ops,
        placement_ops.downward_ops_host,
        placement_ops.downward_level_offsets_host,
        sm_count,
        stream);
}

void UpdateTreeClvsPreservingTipClvs(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    cudaStream_t stream)
{
    const int sm_count = get_device_sm_count();
    rebuild_node_to_tip_map(T, H, placement_ops.node_to_tip);
    build_upward_ops_host_levelized(
        T,
        placement_ops.node_to_tip,
        placement_ops.upward_ops_host,
        placement_ops.upward_level_offsets_host);
    upload_ops_to_device(
        placement_ops,
        placement_ops.upward_ops_host,
        stream);
    launch_upward_clv_update_levelized(
        D,
        placement_ops.d_ops,
        placement_ops.upward_level_offsets_host,
        sm_count,
        stream,
        false);

    build_downward_ops_host_levelized(
        T,
        placement_ops.node_to_tip,
        placement_ops.downward_ops_host,
        placement_ops.downward_level_offsets_host);
    run_downward_clv_stage_levelized(
        D,
        placement_ops,
        placement_ops.downward_ops_host,
        placement_ops.downward_level_offsets_host,
        sm_count,
        stream);
}

void UpdateTreeClvsUpwardOnly(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    cudaStream_t stream)
{
    const int sm_count = get_device_sm_count();
    rebuild_node_to_tip_map(T, H, placement_ops.node_to_tip);
    build_upward_ops_host_levelized(
        T, placement_ops.node_to_tip, placement_ops.upward_ops_host,
        placement_ops.upward_level_offsets_host);
    upload_ops_to_device(placement_ops, placement_ops.upward_ops_host, stream);
    launch_upward_clv_update_levelized(
        D, placement_ops.d_ops, placement_ops.upward_level_offsets_host,
        sm_count, stream);
}

void UpdateTreeClvsUpwardOnlyPrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_upward_ops,
    cudaStream_t stream)
{
    if (!prepared_upward_ops.d_ops || prepared_upward_ops.num_ops <= 0) {
        throw std::runtime_error(
            "UpdateTreeClvsUpwardOnlyPrepared requires uploaded upward operations");
    }
    launch_upward_clv_update_levelized(
        D, prepared_upward_ops.d_ops,
        prepared_upward_ops.upward_level_offsets_host,
        get_device_sm_count(), stream);
}

void PrepareTreeClvOperations(
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& prepared_upward_ops,
    PlacementOpBuffer& prepared_downward_ops,
    cudaStream_t stream)
{
    rebuild_node_to_tip_map(T, H, prepared_upward_ops.node_to_tip);
    build_upward_ops_host_levelized(
        T,
        prepared_upward_ops.node_to_tip,
        prepared_upward_ops.upward_ops_host,
        prepared_upward_ops.upward_level_offsets_host);
    build_downward_ops_host_levelized(
        T,
        prepared_upward_ops.node_to_tip,
        prepared_downward_ops.downward_ops_host,
        prepared_downward_ops.downward_level_offsets_host);
    upload_ops_to_device(
        prepared_upward_ops,
        prepared_upward_ops.upward_ops_host,
        stream);
    upload_ops_to_device(
        prepared_downward_ops,
        prepared_downward_ops.downward_ops_host,
        stream);
}

void UpdateTreeClvsPrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_upward_ops,
    PlacementOpBuffer& prepared_downward_ops,
    cudaStream_t stream)
{
    if (!prepared_upward_ops.d_ops ||
        prepared_upward_ops.num_ops <= 0 ||
        !prepared_downward_ops.d_ops ||
        prepared_downward_ops.num_ops <= 0) {
        throw std::runtime_error(
            "UpdateTreeClvsPrepared requires uploaded upward and downward operations");
    }
    const int sm_count = get_device_sm_count();
    launch_upward_clv_update_levelized(
        D,
        prepared_upward_ops.d_ops,
        prepared_upward_ops.upward_level_offsets_host,
        sm_count,
        stream);
    launch_downward_clv_update_levelized(
        D,
        prepared_downward_ops.d_ops,
        prepared_downward_ops.downward_level_offsets_host,
        sm_count,
        stream);
}

void UpdateTreeClvsDownwardOnlyPrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_downward_ops,
    cudaStream_t stream)
{
    if (!prepared_downward_ops.d_ops ||
        prepared_downward_ops.num_ops <= 0) {
        throw std::runtime_error(
            "UpdateTreeClvsDownwardOnlyPrepared requires uploaded downward operations");
    }
    launch_downward_clv_update_levelized(
        D,
        prepared_downward_ops.d_ops,
        prepared_downward_ops.downward_level_offsets_host,
        get_device_sm_count(),
        stream);
}

void BuildSingleTreeEdgeOutsideWarpSitePrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_downward_ops,
    int operation_index,
    cudaStream_t stream)
{
    if (!prepared_downward_ops.d_ops || operation_index < 0 ||
        operation_index >= prepared_downward_ops.num_ops) {
        throw std::runtime_error(
            "BuildSingleTreeEdgeOutsideWarpSitePrepared: invalid operation index");
    }
    if (!supports_dna_g4_fast_path(D)) {
        constexpr int kGenericBlockSize = 256;
        const unsigned int grid_x = static_cast<unsigned int>(std::max<size_t>(
            1, (D.sites + kGenericBlockSize - 1) / kGenericBlockSize));
        mlipper::likelihood::partials::UpdatePartialsDownwardLevelKernel
            <<<dim3(grid_x, 1), kGenericBlockSize, 0, stream>>>(
                D, prepared_downward_ops.d_ops + operation_index, 1);
        CUDA_CHECK(cudaGetLastError());
        return;
    }
    constexpr int kBlockSize = 256;
    constexpr int kComponentsPerSite = 16;
    const unsigned int grid_x = static_cast<unsigned int>(std::max<size_t>(
        1, (D.sites * kComponentsPerSite + kBlockSize - 1) / kBlockSize));
    mlipper::likelihood::partials::BuildTreeEdgeOutsideWarpSiteKernel
        <<<dim3(grid_x, 1), kBlockSize, 64 * sizeof(fp_t), stream>>>(
            D, prepared_downward_ops.d_ops + operation_index);
    CUDA_CHECK(cudaGetLastError());
}

void UpdateSingleTreeClvUpwardWarpSitePrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_upward_ops,
    int operation_index,
    cudaStream_t stream)
{
    if (!prepared_upward_ops.d_ops || operation_index < 0 ||
        operation_index >= prepared_upward_ops.num_ops) {
        throw std::runtime_error(
            "UpdateSingleTreeClvUpwardWarpSitePrepared: invalid operation index");
    }
    if (!supports_dna_g4_fast_path(D)) {
        constexpr int kGenericBlockSize = 256;
        const unsigned int grid_x = static_cast<unsigned int>(std::max<size_t>(
            1, (D.sites + kGenericBlockSize - 1) / kGenericBlockSize));
        mlipper::likelihood::partials::UpdatePartialsUpwardLevelKernel
            <<<dim3(grid_x, 1), kGenericBlockSize, 0, stream>>>(
                D, prepared_upward_ops.d_ops + operation_index, 1);
        CUDA_CHECK(cudaGetLastError());
        return;
    }
    constexpr int kBlockSize = 256;
    constexpr int kComponentsPerSite = 16;
    const unsigned int grid_x = static_cast<unsigned int>(std::max<size_t>(
        1, (D.sites * kComponentsPerSite + kBlockSize - 1) / kBlockSize));
    mlipper::likelihood::partials::UpdateTreeUpwardWarpSiteKernel
        <<<dim3(grid_x, 1), kBlockSize, 128 * sizeof(fp_t), stream>>>(
            D, prepared_upward_ops.d_ops + operation_index);
    CUDA_CHECK(cudaGetLastError());
}

void RefreshSingleTreeChildDownPrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_downward_ops,
    int operation_index,
    int target_id,
    cudaStream_t stream)
{
    if (!prepared_downward_ops.d_ops || operation_index < 0 ||
        operation_index >= prepared_downward_ops.num_ops) {
        throw std::runtime_error(
            "RefreshSingleTreeChildDownPrepared: invalid operation index");
    }
    if (!supports_dna_g4_fast_path(D)) {
        constexpr int kGenericBlockSize = 256;
        const unsigned int grid_x = static_cast<unsigned int>(std::max<size_t>(
            1, (D.sites + kGenericBlockSize - 1) / kGenericBlockSize));
        mlipper::likelihood::partials::UpdatePartialsDownwardLevelKernel
            <<<dim3(grid_x, 1), kGenericBlockSize, 0, stream>>>(
                D, prepared_downward_ops.d_ops + operation_index, 1);
        CUDA_CHECK(cudaGetLastError());
        return;
    }
    constexpr int kBlockSize = 256;
    constexpr int kComponentsPerSite = 16;
    const unsigned int grid_x = static_cast<unsigned int>(std::max<size_t>(
        1, (D.sites * kComponentsPerSite + kBlockSize - 1) / kBlockSize));
    mlipper::likelihood::partials::RefreshTreeChildDownWarpSiteKernel
        <<<dim3(grid_x, 1), kBlockSize, 64 * sizeof(fp_t), stream>>>(
            D, target_id);
    CUDA_CHECK(cudaGetLastError());
}

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

    run_downward_clv_stage(
        D,
        placement_ops,
        required_downward_ops,
        sm_count,
        stream);
}

static void UpdateTreeClvsAfterInsertion(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    int upward_start_node,
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

    build_downward_ops_host_levelized(
        T,
        placement_ops.node_to_tip,
        placement_ops.downward_ops_host,
        placement_ops.downward_level_offsets_host);
    run_downward_clv_stage_levelized(
        D,
        placement_ops,
        placement_ops.downward_ops_host,
        placement_ops.downward_level_offsets_host,
        sm_count,
        stream);
}

static InsertResult insert_query_with_intermediate(
    TreeBuildResult& T,
    const std::string& raw_name,
    int target_id,
    double pendant,
    double proximal)
{
    InsertResult out{};
    if (target_id < 0 ||
        target_id >= static_cast<int>(T.nodes.size())) {
        throw std::runtime_error("insert_query_with_intermediate: invalid target_id.");
    }
    if (T.nodes.size() >
        static_cast<size_t>(std::numeric_limits<int>::max() - 2)) {
        throw std::length_error(
            "insert_query_with_intermediate: tree exceeds integer indexing limits");
    }
    int parent_id = T.nodes[target_id].parent;
    const double total = T.nodes[target_id].branch_length_to_parent;
    double proximal_len = 0.0;
    double distal_len = 0.0;
    normalize_split_branch_lengths(total, proximal, OPT_BRANCH_LEN_MIN, proximal_len, distal_len);
    double pendant_len = sanitize_branch_length(pendant);

    const int new_internal_id = static_cast<int>(T.nodes.size());
    const int new_tip_id = new_internal_id + 1;

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

static void append_query_tip_to_host_packing(
    HostPacking& H,
    const PlacementQueryBatch& Q,
    int query_idx,
    int new_tip_node_id,
    size_t sites)
{
    if (query_idx < 0 || static_cast<size_t>(query_idx) >= Q.count) {
        throw std::runtime_error("append_query_tip_to_host_packing: query_idx out of range.");
    }
    const size_t needed = mlipper::util::checked_mul_size(
        static_cast<size_t>(query_idx) + 1, sites,
        "committed query character prefix");
    if (Q.query_chars.size() < needed) {
        throw std::runtime_error("append_query_tip_to_host_packing: query_chars buffer too small.");
    }

    H.tip_node_ids.push_back(new_tip_node_id);
    const size_t old = H.tipchars.size();
    H.tipchars.resize(mlipper::util::checked_add_size(
        old, sites, "committed tip character matrix"));
    std::memcpy(
        H.tipchars.data() + old,
        Q.query_chars.data() + static_cast<size_t>(query_idx) * sites,
        sites * sizeof(uint8_t));
}

static void update_insertion_device(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
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

    copy_host_to_device_async(
        D.d_tipchars + (size_t)old_tips * D.sites,
        H.tipchars.data() + (size_t)old_tips * D.sites,
        sizeof(uint8_t) * D.sites,
        stream);
    copy_host_to_device_async(
        D.d_tip_node_ids + (size_t)old_tips,
        H.tip_node_ids.data() + (size_t)old_tips,
        sizeof(int),
        stream);

    const size_t stride = D.pmat_per_node_elems();
    const fp_t* mid_prox_src = H.pmats_mid_prox.empty()
        ? H.pmats_mid.data()
        : H.pmats_mid_prox.data();
    const fp_t* mid_dist_src = H.pmats_mid_dist.empty()
        ? H.pmats_mid.data()
        : H.pmats_mid_dist.data();

    auto sync_node = [&](int nid) {
        if (nid < 0 || nid >= newN) {
            throw std::runtime_error("update_insertion_device: node id out of range.");
        }
        const size_t off = (size_t)nid * stride;
        copy_host_to_device_async(
            D.d_blen + (size_t)nid,
            &H.blen[(size_t)nid],
            sizeof(fp_t),
            stream);

        auto copy_pmat_slice = [&](fp_t* dst, const fp_t* src) {
            copy_host_to_device_async(
                dst + off, src + off, sizeof(fp_t) * stride, stream);
        };
        copy_pmat_slice(D.d_pmat, H.pmats.data());
        copy_pmat_slice(D.d_pmat_mid, H.pmats_mid.data());
        copy_pmat_slice(D.d_pmat_mid_prox, mid_prox_src);
        copy_pmat_slice(D.d_pmat_mid_dist, mid_dist_src);
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
    const RawPlacementResult& placement,
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
    const double total_branch_length = tree.nodes[placement.target_id].branch_length_to_parent;
    const double commit_proximal_length = total_branch_length - placement.distal_length;

    // Topology is authoritative. Rebuild each derived representation in order
    // before any subsequent query is allowed to observe the inserted tip.
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

    populate_host_topology(tree, host);
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
        placement.target_id,
        insert_result.internal_id,
        insert_result.tip_id,
        stream);

    UpdateTreeClvsAfterInsertion(
        D,
        tree,
        host,
        placement_ops,
        insert_result.internal_id,
        stream);
}

class MainPlacementScratch {
public:
    explicit MainPlacementScratch(cudaStream_t stream) : stream_(stream) {}

    ~MainPlacementScratch() noexcept
    {
        if (device_id_ >= 0) cudaSetDevice(device_id_);
        // DeviceBuffer destruction uses cudaFree. Synchronize the stream first
        // so no placement kernel still references these borrowed scratch views.
        if (sumtable_ || likelihoods_ || query_pmat_) {
            cudaStreamSynchronize(stream_);
        }
    }

    MainPlacementScratch(const MainPlacementScratch&) = delete;
    MainPlacementScratch& operator=(const MainPlacementScratch&) = delete;

    void ensureCapacity(const DeviceTree& D, size_t required_ops)
    {
        if (required_ops == 0) return;
        ensure_device_tree_current_device(D, "MainPlacementScratch::ensureCapacity");

        const size_t sumtable_stride = D.per_node_elems();
        const size_t query_pmat_stride = D.pmat_per_node_elems();
        if (sumtable_stride == 0 || query_pmat_stride == 0) {
            throw std::runtime_error(
                "MainPlacementScratch::ensureCapacity: invalid scratch stride.");
        }

        constexpr size_t kCapacityChunkOps = 256;
        const size_t adjusted_ops = mlipper::util::checked_add_size(
            required_ops,
            kCapacityChunkOps - 1,
            "main placement scratch rounding");
        const size_t target_ops = mlipper::util::checked_mul_size(
            adjusted_ops / kCapacityChunkOps,
            kCapacityChunkOps,
            "main placement scratch rounding");
        sumtable_.ensureCapacity(mlipper::util::checked_mul_size(
            sumtable_stride, target_ops, "main placement sumtable"));
        likelihoods_.ensureCapacity(target_ops);
        query_pmat_.ensureCapacity(mlipper::util::checked_mul_size(
            query_pmat_stride, target_ops, "main placement query matrices"));

        device_id_ = D.device_id;
        view_.d_sumtable = sumtable_.get();
        view_.d_likelihoods = likelihoods_.get();
        view_.d_query_pmat = query_pmat_.get();
        view_.sumtable_capacity_ops = target_ops;
        view_.likelihood_capacity_ops = target_ops;
        view_.query_pmat_capacity_ops = target_ops;
    }

    const PlacementScratchOverride& view() const noexcept { return view_; }

private:
    cudaStream_t stream_ = nullptr;
    int device_id_ = -1;
    mlipper::gpu::DeviceBuffer<fp_t> sumtable_;
    mlipper::gpu::DeviceBuffer<fp_t> likelihoods_;
    mlipper::gpu::DeviceBuffer<fp_t> query_pmat_;
    PlacementScratchOverride view_{};
};

static RawPlacementResult evaluate_single_placement_query(
    DeviceTree& D,
    const PlacementOpBuffer& placement_ops,
    int query_idx,
    int smoothing,
    cudaStream_t stream,
    const PlacementScratchOverride* scratch_override)
{
    const size_t node_bytes =
        mlipper::util::checked_allocation_bytes<fp_t>(
            static_cast<size_t>(D.N), "placement branch-length scratch");
    zero_length_scratch_async(D, node_bytes, stream);

    build_query_clv(D, query_idx, stream);
    CUDA_CHECK(cudaGetLastError());

    DeviceTree query_view = make_query_view(D, query_idx);
    return EvaluatePlacementCandidates(
        query_view,
        placement_ops.d_ops,
        placement_ops.num_ops,
        smoothing,
        stream,
        true,
        placement_ops.tuning,
        scratch_override);
}

void EvaluatePlacementQueries(
    DeviceTree& D,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    PlacementCommitContext& commit_ctx,
    std::vector<RawPlacementResult>* placement_results_out,
    int smoothing,
    bool commit_to_tree,
    cudaStream_t stream)
{
    const int query_count = D.placement_queries;
    if (placement_results_out) {
        placement_results_out->clear();
        placement_results_out->reserve(static_cast<size_t>(query_count));
    }

    if (!commit_ctx.placement_ops) {
        throw std::runtime_error("EvaluatePlacementQueries: placement ops are null.");
    }
    if (commit_to_tree) {
        validate_commit_context(commit_ctx);
    }

    const PlacementOpBuffer& placement_ops = *commit_ctx.placement_ops;
    MainPlacementScratch placement_scratch(stream);
    for (int query_idx = 0; query_idx < query_count; ++query_idx) {
        placement_scratch.ensureCapacity(
            D,
            static_cast<size_t>(placement_ops.num_ops));
        const PlacementScratchOverride& scratch_view = placement_scratch.view();
        RawPlacementResult placement = evaluate_single_placement_query(
            D,
            placement_ops,
            query_idx,
            smoothing,
            stream,
            &scratch_view);
        if (placement_results_out) {
            placement_results_out->push_back(placement);
        }
        if (commit_to_tree) {
            // Commit refreshes topology, PMATs, traversal operations, and CLVs;
            // the following iteration therefore scores against the new tree.
            commit_placement_result(
                D,
                er,
                rate_multipliers,
                commit_ctx,
                placement,
                query_idx,
                stream);
        }
    }
}
