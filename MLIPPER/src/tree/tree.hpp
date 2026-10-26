#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cuda_runtime.h>
#include <libpll/pll.h>

#include "pmatrix/pmat.h"
#include "util/checked_size.hpp"
#include "util/precision.hpp"

struct RawPlacementResult;
struct NodeOpInfo;

struct PlacementTuningConfig {
    static constexpr int kDefaultFullOptPasses = 4;
    static constexpr int kDefaultExportPlacementTopK = 5;
#if !defined(MLIPPER_USE_DOUBLE)
    static constexpr int kDefaultDoubleRerankUlpFactor = 4;
#endif

    int full_opt_passes = kDefaultFullOptPasses;
    // A positive accumulated-LWR threshold supersedes this fixed limit.
    int export_placement_topk = kDefaultExportPlacementTopK;
    double export_accumulated_lwr_threshold = 0.0;
#if !defined(MLIPPER_USE_DOUBLE)
    bool enable_double_rerank = true;
    int double_rerank_ulp_factor = kDefaultDoubleRerankUlpFactor;
    double double_rerank_gap_top2 = 0.0;
#endif
    double pendant_branch_min = 1.0e-4;
    double split_branch_min = 1.0e-4;
};

struct PlacementOpBuffer {
    // d_ops is owned by this aggregate by convention and released through
    // free_placement_op_buffer(); the host vectors define reusable schedules.
    NodeOpInfo* d_ops = nullptr;
    int num_ops = 0, capacity = 0;
    std::vector<int> node_to_tip;
    std::vector<NodeOpInfo> upward_ops_host;
    // Exclusive offsets for upward ops grouped by height from the tips.
    std::vector<int> upward_level_offsets_host;
    std::vector<NodeOpInfo> downward_ops_host;
    // Exclusive offsets for downward ops grouped by child depth.
    std::vector<int> downward_level_offsets_host;
    PlacementTuningConfig tuning;
};

namespace parse {
struct ModelConfig;
}

class CudaRuntimeError : public std::runtime_error {
public:
    CudaRuntimeError(cudaError_t code, std::string message)
        : std::runtime_error(std::move(message)), code_(code)
    {
    }

    cudaError_t code() const noexcept { return code_; }

private:
    cudaError_t code_;
};

#ifndef CUDA_CHECK
#define CUDA_CHECK(expr) do { \
    cudaError_t _e = (expr);  \
    if (_e != cudaSuccess) {  \
      throw CudaRuntimeError(_e, std::string("[CUDA] ") + cudaGetErrorString(_e) + \
          " at " + __FILE__ + ":" + std::to_string(__LINE__) + \
          " while calling " + #expr); \
    }                         \
} while(0)
#endif

inline uint8_t encode_state_DNA5(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't':
        case 'U': case 'u': return 3;
        case '-':           return 4; // gap
        default:            return 4; // Unknown or ambiguous state.
    }
}

inline uint8_t encode_state_DNA4_mask(char c) {
    switch (c) {
        case 'A': case 'a': return 1u << 0;
        case 'C': case 'c': return 1u << 1;
        case 'G': case 'g': return 1u << 2;
        case 'T': case 't':
        case 'U': case 'u': return 1u << 3;
        case 'R': case 'r': return (1u << 0) | (1u << 2); // A/G
        case 'Y': case 'y': return (1u << 1) | (1u << 3); // C/T
        case 'S': case 's': return (1u << 1) | (1u << 2); // C/G
        case 'W': case 'w': return (1u << 0) | (1u << 3); // A/T
        case 'K': case 'k': return (1u << 2) | (1u << 3); // G/T
        case 'M': case 'm': return (1u << 0) | (1u << 1); // A/C
        case 'B': case 'b': return (1u << 1) | (1u << 2) | (1u << 3); // C/G/T
        case 'D': case 'd': return (1u << 0) | (1u << 2) | (1u << 3); // A/G/T
        case 'H': case 'h': return (1u << 0) | (1u << 1) | (1u << 3); // A/C/T
        case 'V': case 'v': return (1u << 0) | (1u << 1) | (1u << 2); // A/C/G
        case 'N': case 'n':
        case '-':
        case '.':
        default:            return (1u << 0) | (1u << 1) | (1u << 2) | (1u << 3); // unknown/gap
    }
}

namespace mlipper {

struct SequenceRecord {
    std::string name;
    std::string sequence;
};

} // namespace mlipper


// GPU layout assumes every node has the same CLV size (sites * rate_cats * states).
struct TreeNode {
    int   id = -1;
    // Stable across topology edits; unlike id, this is not an array index.
    // Local refinement uses it to find committed attachment and query nodes.
    int   stable_node_label = -1;
    bool  is_tip = false;
    int   left = -1, right = -1, parent = -1;
    fp_t branch_length_to_parent = fp_t(0); // Number after the colon in Newick
    std::string name;   // only for tips
};

// Authoritative mutable CPU topology. All array-index relationships use node
// IDs; HostPacking and DeviceTree are derived representations of this state.
struct TreeBuildResult {
    std::vector<TreeNode> nodes; // 0..N-1
    int root_id = -1;
    std::vector<int> postorder; // children -> parent
    std::vector<int> preorder;  // parent -> children
    std::unordered_map<std::string,int> tip_node_by_name;
};

// Parse a rooted Newick tree and assign dense node IDs. Tip names must match
// msa_tip_names exactly; the returned preorder/postorder and parent/child links
// are mutually consistent and branch lengths are normalized for MLIPPER use.
TreeBuildResult build_tree_from_newick_with_pll(
    const std::vector<std::string>& msa_tip_names,
    const std::string& newick_text);

std::vector<double> build_gtr_q_matrix(
    int states,
    const parse::ModelConfig& model,
    const std::vector<double>& pi);

// Trivially copyable CUDA descriptor passed to kernels. Pointer fields do not
// express ownership; OwnedDeviceTree releases the allocations it owns.
enum class DownwardPmatIndexing : std::uint8_t {
    Columns,
    Rows,
};

struct DeviceTree {
    int     N = 0, tips = 0, inners = 0, placement_queries = 0;
    int     capacity_N = 0, capacity_tips = 0, query_capacity = 0;
    int     root_id = -1, device_id = -1;
    size_t  sites = 0;
    int     states = 0, rate_cats = 0;
    unsigned int log2_stride = 0;
    bool    per_rate_scaling = false;

    fp_t   *d_lambdas = nullptr, *d_V = nullptr, *d_Vinv = nullptr;
    fp_t   *d_frequencies = nullptr, *d_rate_weights = nullptr;

    fp_t   *d_blen = nullptr, *d_new_pendant_length = nullptr, *d_new_proximal_length = nullptr;
    fp_t   *d_prev_pendant_length = nullptr, *d_prev_proximal_length = nullptr;

    uint8_t *d_tipchars = nullptr; // [tips * sites], DNA4 bitmask when states==4 else DNA5 code
    int     *d_tip_node_ids = nullptr; // [tips], maps tip index -> node id for tip-CLV initialization

    fp_t    *d_clv_up = nullptr, *d_clv_down = nullptr;
    fp_t    *d_edge_midpoint_clv = nullptr, *d_edge_outside_clv = nullptr;
    unsigned *d_pattern_weights_u = nullptr;

    unsigned *d_scaler_storage = nullptr;
    unsigned *d_site_scaler_up = nullptr, *d_site_scaler_down = nullptr;
    unsigned *d_edge_midpoint_scaler = nullptr, *d_edge_outside_scaler = nullptr;

    fp_t* d_pmat = nullptr, *d_pmat_mid = nullptr;
    fp_t* d_pmat_mid_prox = nullptr, *d_pmat_mid_dist = nullptr;
    unsigned int* d_tipmap = nullptr;       // decode table for tip/query chars -> state bitmask

    uint8_t* d_query_chars = nullptr; // [num_queries * sites], same encoding contract as d_tipchars
    fp_t   *d_query_clv = nullptr, *d_query_pmat = nullptr;
    DownwardPmatIndexing downward_pmat_indexing =
        DownwardPmatIndexing::Columns;

    size_t per_node_elems() const {
        return mlipper::util::checked_product(
            "DeviceTree CLV stride", sites, static_cast<size_t>(rate_cats),
            static_cast<size_t>(states));
    }
    size_t scaler_elems() const {
        return per_rate_scaling
            ? mlipper::util::checked_mul_size(
                  sites, static_cast<size_t>(rate_cats), "DeviceTree scaler stride")
            : sites;
    }
    size_t scaler_pool_elems() const {
        return mlipper::util::checked_mul_size(
            static_cast<size_t>(capacity_N), scaler_elems(),
            "DeviceTree scaler pool");
    }
    size_t scaler_storage_elems() const {
        return mlipper::util::checked_mul_size(
            scaler_pool_elems(), 4, "DeviceTree scaler storage");
    }
    size_t pmat_per_node_elems() const {
        const size_t state_count = static_cast<size_t>(states);
        return mlipper::util::checked_product(
            "DeviceTree transition-matrix stride",
            static_cast<size_t>(rate_cats), state_count, state_count);
    }
};

static_assert(
    std::is_trivially_copyable_v<DeviceTree>,
    "DeviceTree must remain a trivially-copyable CUDA view");

// Owns the allocations referenced by its DeviceTree base. DeviceTree remains
// a trivially-copyable, non-owning kernel descriptor; only this type releases
// the pointers. Inheritance keeps existing kernel call sites as plain
// DeviceTree references without adding conversion machinery.
struct OwnedDeviceTree : DeviceTree {
    OwnedDeviceTree() = default;
    ~OwnedDeviceTree() noexcept;

    OwnedDeviceTree(const OwnedDeviceTree&) = delete;
    OwnedDeviceTree& operator=(const OwnedDeviceTree&) = delete;
    OwnedDeviceTree(OwnedDeviceTree&& other) noexcept;
    OwnedDeviceTree& operator=(OwnedDeviceTree&& other) noexcept;

    void reset() noexcept;
};

static_assert(!std::is_copy_constructible_v<OwnedDeviceTree>);
static_assert(!std::is_copy_assignable_v<OwnedDeviceTree>);
static_assert(std::is_nothrow_move_constructible_v<OwnedDeviceTree>);
static_assert(std::is_nothrow_move_assignable_v<OwnedDeviceTree>);

inline void ensure_device_tree_current_device(
    const DeviceTree& D,
    const char* context)
{
    if (D.device_id < 0) {
        return;
    }
    int current_device = -1;
    CUDA_CHECK(cudaGetDevice(&current_device));
    if (current_device == D.device_id) {
        return;
    }
    const cudaError_t err = cudaSetDevice(D.device_id);
    if (err != cudaSuccess) {
        throw std::runtime_error(
            std::string("[CUDA] failed to switch to DeviceTree device ") +
            std::to_string(D.device_id) + " from device " +
            std::to_string(current_device) + " in " + context + ": " +
            cudaGetErrorString(err));
    }
}

// Contiguous host-side transfer representation derived from TreeBuildResult.
// Topology edits must update TreeBuildResult first and then rebuild this data.
struct HostPacking {
    std::vector<int>     postorder, preorder, parent, left, right;
    std::vector<uint8_t> is_tip;
    std::vector<fp_t>    blen;

    std::vector<int>     tip_node_ids;      // size = tips
    std::vector<uint8_t> tipchars;          // size = tips * sites

    std::vector<unsigned> pattern_weights;  // size = sites
    std::vector<fp_t>     pmats;
    std::vector<fp_t>     pmats_mid;        // half-branch pmats for midpoint calculations
    std::vector<fp_t>     pmats_mid_prox;   // proximal branch pmats for midpoint calculations
    std::vector<fp_t>     pmats_mid_dist;   // distal branch pmats for midpoint calculations
};

struct PlacementQueryBatch {
    size_t count = 0;
    // Query-major encoded rows: [count][sites].
    std::vector<uint8_t> query_chars;

    bool empty() const { return count == 0; }
    size_t size() const { return count; }
};

// Non-owning references required to commit a placement across every mutable
// CPU/host representation. They may be null only for scoring without commit.
struct PlacementCommitContext {
    TreeBuildResult* tree = nullptr;
    HostPacking* host = nullptr;
    PlacementQueryBatch* queries = nullptr;
    PlacementOpBuffer* placement_ops = nullptr;
    const std::vector<std::string>* query_names = nullptr;
    std::vector<std::string>* inserted_query_names = nullptr;
};

// Maps one full-tree NNI target into a compact scoring slot. When
// direct_edge_outside_src is valid its cached outside message is copied;
// otherwise parent_down_src and sibling_up_src are combined through the sibling
// and second PMATs to reconstruct the same message in dst_target.
struct DirectNNIContextOp {
    int target_src = -1;
    int direct_edge_outside_src = -1;
    int parent_down_src = -1;
    int sibling_up_src = -1;
    int second_pmat_src = -1;
    int dst_target = -1;
};

struct OwnedSubtreeWorkspace {
    OwnedDeviceTree dev{};
};

struct SubtreeWorkspaceLoadConfig {
    const TreeBuildResult& tree;
    const HostPacking& host;
    const EigResult& eig;
    const std::vector<double>& rate_weights;
    const std::vector<double>& rate_multipliers;
    const std::vector<double>& pi;
    size_t sites = 0;
    int states = 0;
    int rate_cats = 0;
    bool per_rate_scaling = false;
    const PlacementQueryBatch* queries = nullptr;
    bool commit_to_tree = false;
    // Explicit capacities. Callers must reserve at least the loaded query count.
    int query_capacity = 0;
    // Number of future insertions for which node/tip/CLV storage is reserved.
    int insert_capacity = 0;
};

HostPacking pack_host_arrays_from_tree_and_msa(
    const TreeBuildResult& T,
    const std::vector<std::string>& msa_tip_names,
    const std::vector<std::string>& msa_rows,
    size_t sites,
    int states
);
void fill_pmats_in_host_packing(
    const TreeBuildResult&       T,
    HostPacking&                 H,
    const EigResult&             er,
    const std::vector<double>&   rate_multipliers,   // len = rate_cats (per-category rate multipliers)
    int states,
    int rate_cats,
    const int* changed_nodes = nullptr,
    int num_changed_nodes = 0,
    bool include_midpoint_pmats = true
);

// Allocate and populate an owning device tree on the caller's current CUDA
// device. query storage is sized for the loaded batch; insert_capacity reserves
// additional node/tip slots for commit mode. On failure device_tree remains
// empty. Kernel-facing DeviceTree copies are non-owning views of these buffers.
void allocate_device_tree_on_current_gpu(
    OwnedDeviceTree& device_tree,
    const TreeBuildResult& tree,
    const HostPacking& host,
    const EigResult& eig,
    const std::vector<double>& rate_weights,
    const std::vector<double>& rate_multipliers,
    const std::vector<double>& frequencies,
    size_t sites,
    int states,
    int rate_cats,
    bool per_rate_scaling,
    const PlacementQueryBatch* queries,
    bool commit_to_tree,
    int insert_capacity,
    bool allocate_directional_clvs = true);

DeviceTree make_query_view(const DeviceTree& D, int query_idx);
void copy_unscaled_up_clv_to_query_slot(
    const DeviceTree& src,
    int src_node_id,
    DeviceTree& dst,
    int dst_query_idx,
    cudaStream_t stream = 0);
void copy_upward_state(const DeviceTree& src, DeviceTree& dst, cudaStream_t stream = 0);
void copy_selected_upward_state(
    const DeviceTree& src,
    DeviceTree& dst,
    const int* d_node_ids,
    int node_count,
    cudaStream_t stream = 0);
void build_direct_nni_target_contexts(
    const DeviceTree& src,
    DeviceTree& dst,
    const std::vector<DirectNNIContextOp>& ops,
    cudaStream_t stream = 0);

void reload_device_tree_live_data(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const PlacementQueryBatch* queries = nullptr,
    cudaStream_t stream = 0);
void reload_device_tree_live_data_preserving_clvs(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const PlacementQueryBatch* queries = nullptr,
    cudaStream_t stream = 0);
void reload_device_tree_live_data_local_spr(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const HostPacking& base_H,
    int current_main_pmat_node,
    int& previous_main_pmat_node,
    const PlacementQueryBatch* queries = nullptr,
    cudaStream_t stream = 0);
bool subtree_workspace_requires_rebuild(
    const OwnedSubtreeWorkspace& workspace,
    const TreeBuildResult& tree,
    const HostPacking& host,
    size_t sites,
    int states,
    int rate_cats,
    bool per_rate_scaling);
// Ensures storage for exactly required_queries active slots. Reallocation does
// not preserve existing query characters or CLVs; callers must reload them.
void ensure_device_tree_query_capacity(
    DeviceTree& D,
    int required_queries,
    const char* context);
void release_subtree_workspace(OwnedSubtreeWorkspace& workspace) noexcept;
void load_subtree_workspace(
    OwnedSubtreeWorkspace& workspace,
    const SubtreeWorkspaceLoadConfig& config,
    cudaStream_t stream = 0,
    const char* context = "load_subtree_workspace");

struct TopologyRefinementState {
    // Non-owning device view. GPU buffers remain owned by MlipperSession or a
    // session workspace; this aggregate must never free them.
    DeviceTree device;
    TreeBuildResult tree;
    HostPacking host_packing;
    EigResult eig;
    PlacementQueryBatch queries;
};

void launch_init_tip_clv(const DeviceTree& D, cudaStream_t stream = 0);

void free_placement_op_buffer(PlacementOpBuffer& placement_ops, cudaStream_t stream = 0);

void UploadPlacementOps(
    PlacementOpBuffer& placement_ops, const std::vector<NodeOpInfo>& host_ops, cudaStream_t stream = 0);

// Rebuild the full tree CLV state after an initial upload or global topology change.
void UpdateTreeClvs(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    cudaStream_t stream = 0);

// Rebuild internal upward/downward CLVs while preserving the CLVs currently
// installed in tip slots. This supports virtual D&C boundary tips carrying
// precomputed directional messages instead of sequence characters.
void UpdateTreeClvsPreservingTipClvs(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    cudaStream_t stream = 0);

// Root-likelihood/model-optimization path: rebuild only postorder/upward CLVs.
// Reuses the placement-op buffer and does not allocate or touch down/mid CLVs.
void UpdateTreeClvsUpwardOnly(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    cudaStream_t stream = 0);

void UpdateTreeClvsUpwardOnlyPrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_upward_ops,
    cudaStream_t stream = 0);

void PrepareTreeClvOperations(
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& prepared_upward_ops,
    PlacementOpBuffer& prepared_downward_ops,
    cudaStream_t stream = 0);

void UpdateTreeClvsPrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_upward_ops,
    PlacementOpBuffer& prepared_downward_ops,
    cudaStream_t stream = 0);

// Complete a derivative pass when the current candidate's upward CLVs are
// already resident in D.
void UpdateTreeClvsDownwardOnlyPrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_downward_ops,
    cudaStream_t stream = 0);

// Sequential branch-coordinate traversal primitives. Each call refreshes one
// directional message across all sites without rebuilding the full tree.
void BuildSingleTreeEdgeOutsideWarpSitePrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_downward_ops,
    int operation_index,
    cudaStream_t stream = 0);
void UpdateSingleTreeClvUpwardWarpSitePrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_upward_ops,
    int operation_index,
    cudaStream_t stream = 0);
void RefreshSingleTreeChildDownPrepared(
    DeviceTree& D,
    PlacementOpBuffer& prepared_downward_ops,
    int operation_index,
    int target_id,
    cudaStream_t stream = 0);

// Recompute only the CLVs affected by a prune/local-regraft scoring pass.
// upward_start_node < 0 skips the upward path refresh and only applies the
// supplied downward ops.
void UpdateTreeClvsAfterPrune(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
    int upward_start_node,
    const std::vector<NodeOpInfo>& required_downward_ops,
    cudaStream_t stream = 0);

// Evaluates loaded queries in order. With commit_to_tree=false, D and the CPU
// tree remain unchanged. With commit_to_tree=true, each accepted placement is
// inserted and all dependent host/device traversal state is refreshed before
// scoring the next query.
void EvaluatePlacementQueries(
    DeviceTree& D,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    PlacementCommitContext& commit_ctx,
    std::vector<struct RawPlacementResult>* placement_results_out = nullptr,
    int smoothing = 1,
    bool commit_to_tree = true,
    cudaStream_t stream = 0);

std::vector<mlipper::SequenceRecord> build_placement_query(
    const std::vector<std::string>& msa_tip_names,
    const std::vector<std::string>& msa_rows);
