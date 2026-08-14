#pragma once
#include <cuda_runtime.h>
#include <utility>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <cstddef>
#include <cstdint>
#include <libpll/pll.h>

#include "pmatrix/pmat.h"
#include "util/precision.hpp"

struct PlacementResult;
struct NodeOpInfo;

struct CommitTimingStats {
    double initial_upward_host_ms = 0.0;
    double initial_downward_host_ms = 0.0;
    double initial_upward_stage_ms = 0.0;
    double initial_downward_stage_ms = 0.0;

    double query_reset_stage_ms = 0.0;
    double query_build_clv_stage_ms = 0.0;
    double query_kernel_total_ms = 0.0;

    double insertion_pre_clv_ms = 0.0;
    double insertion_upward_host_ms = 0.0;
    double insertion_downward_host_ms = 0.0;
    double insertion_upward_stage_ms = 0.0;
    double insertion_downward_stage_ms = 0.0;

    long long initial_upward_ops = 0;
    long long initial_downward_ops = 0;
    long long insertion_upward_ops = 0;
    long long insertion_downward_ops = 0;

    int initial_updates = 0;
    int query_evals = 0;
    int insertion_updates = 0;
};

struct PlacementOpBuffer {
    NodeOpInfo* d_ops = nullptr;
    int num_ops = 0;
    int capacity = 0;
    std::vector<int> node_to_tip;
    std::vector<NodeOpInfo> upward_ops_host;
    std::vector<NodeOpInfo> downward_ops_host;
    // Exclusive offsets for downward ops grouped by child depth.
    std::vector<int> downward_level_offsets_host;
    bool profile_commit_timing = false;
    CommitTimingStats timing;
};

namespace parse {
struct ModelConfig;
}

// ===== CUDA error helper =====
#ifndef CUDA_CHECK
#define CUDA_CHECK(expr) do { \
    cudaError_t _e = (expr);  \
    if (_e != cudaSuccess) {  \
      throw std::runtime_error(std::string("[CUDA] ") + cudaGetErrorString(_e)); \
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
        default:            return 4; // Treat as gap for now; switch to bitmask if you need IUPAC
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

struct NewPlacementQuery {
    std::pair<int,int> node_id_pair{-1, -1}; // node where to insert {first = parent, second = child}
    fp_t pendant = fp_t(0);                  // branch length of new insertion
    fp_t distal = fp_t(0);                   // branch length to new insertion
    std::string msa;                         // sequence
    std::string msa_name;                    // sequence name
};


// GPU layout assumes every node has the same CLV size (sites * rate_cats * states).
struct TreeNode {
    int   id = -1;
    bool  is_tip = false;
    int   left = -1;    // child id
    int   right = -1;   // child id
    int   parent = -1;  // parent id (root = -1)
    fp_t branch_length_to_parent = fp_t(0); // Number after the colon in Newick
    fp_t branch_length_to_insert = fp_t(0); // for insertion operations
    std::string name;   // only for tips
    // GPU offsets: use these directly with scaler pool
    size_t scaler_offset = 0;  // elements-based offset (per-site or per-rate)
};

struct TreeBuildResult {
    std::vector<TreeNode> nodes; // 0..N-1
    int root_id = -1;
    std::vector<int> postorder; // children -> parent
    std::vector<int> preorder;  // parent -> children
    std::unordered_map<std::string,int> tip_node_by_name;
};

// ----- Tree/model construction -----

TreeBuildResult build_tree_from_newick_with_pll(
    const std::vector<std::string>& msa_tip_names,
    const std::string& newick_text,
    size_t sites,
    int states,
    int rate_cats,
    bool per_rate_scaling);

std::vector<double> build_mixture_weights(const parse::ModelConfig& model, int rate_cats);
std::vector<double> build_gamma_rate_categories(double alpha, int rate_cats);
std::vector<double> build_gtr_q_matrix(
    int states,
    const parse::ModelConfig& model,
    const std::vector<double>& pi);

// ----- Device-side tree storage -----

struct DeviceTree {
    int     N = 0; // nodes = tips + inners
    int     tips = 0;
    int     inners = 0;
    int     placement_queries = 0;
    int     capacity_N = 0;
    int     capacity_tips = 0;
    int     query_capacity = 0;
    int     root_id = -1;
    size_t  sites = 0;
    int     states = 0;
    int     rate_cats = 0;
    unsigned int log2_stride = 0;
    bool    per_rate_scaling = false;

    fp_t   *d_lambdas = nullptr;      // [rate_cats * states]
    fp_t   *d_V = nullptr;            // [states * states] row-major
    fp_t   *d_Vinv = nullptr;         // [states * states] row-major
    fp_t   *d_rate_w = nullptr;       // [rate_cats]
    fp_t   *d_frequencies = nullptr;  // [states]
    fp_t   *d_rate_weights = nullptr; // [rate_cats]

    fp_t   *d_blen = nullptr;                 // [N]
    fp_t   *d_new_pendant_length = nullptr;  // [N]
    fp_t   *d_new_proximal_length = nullptr; // [N]
    fp_t   *d_prev_pendant_length = nullptr; // [N]
    fp_t   *d_prev_proximal_length = nullptr; // [N] previous proximal branch lengths (smoothing rollback)

    uint8_t *d_tipchars = nullptr; // [tips * sites], DNA4 bitmask when states==4 else DNA5 code
    int     *d_tip_node_ids = nullptr; // [tips], maps tip index -> node id for tip-CLV initialization

    fp_t    *d_clv_up = nullptr;
    fp_t    *d_clv_down = nullptr;
    size_t   clv_down_offset_elems = 0;
    fp_t    *d_clv_mid = nullptr;
    fp_t    *d_clv_mid_base = nullptr;
    double  *d_root_loglik_total = nullptr;
    fp_t    *d_sumtable = nullptr;
    fp_t    *d_likelihoods = nullptr;
    size_t   sumtable_capacity_ops = 0;
    size_t   likelihood_capacity_ops = 0;
    unsigned *d_pattern_weights_u = nullptr;

    // Root-likelihood compatibility view. Do not alias this to one of the
    // branch-local scaler pools globally; set it explicitly only for the
    // root-specific path that still expects a flat site-scaler pointer.
    unsigned *d_site_scaler = nullptr;
    unsigned *d_site_scaler_storage = nullptr;
    unsigned *d_site_scaler_up = nullptr;
    unsigned *d_site_scaler_down = nullptr;
    unsigned *d_site_scaler_mid = nullptr;
    unsigned *d_site_scaler_mid_base = nullptr;

    fp_t* d_pmat = nullptr;
    fp_t* d_pmat_mid = nullptr;             // half-branch pmats for midpoint calculations
    fp_t* d_pmat_mid_prox = nullptr;        // proximal branch pmats for midpoint calculations
    fp_t* d_pmat_mid_dist = nullptr;        // distal branch pmats for midpoint calculations
    unsigned int* d_tipmap = nullptr;       // decode table for tip/query chars -> state bitmask

    uint8_t* d_query_chars = nullptr; // [num_queries * sites], same encoding contract as d_tipchars
    fp_t   *d_query_clv = nullptr;    // [sites * rate_cats * states]
    fp_t   *d_query_pmat = nullptr;

    size_t per_node_elems() const {
        return sites * static_cast<size_t>(rate_cats) * static_cast<size_t>(states);
    }
    size_t clv_pool_elems() const {
        return static_cast<size_t>(N) * per_node_elems();
    }
    size_t scaler_elems() const {
        return per_rate_scaling ? (sites * static_cast<size_t>(rate_cats)) : sites;
    }
    size_t scaler_pool_elems() const {
        return static_cast<size_t>(capacity_N) * scaler_elems();
    }
    size_t scaler_storage_elems() const {
        return scaler_pool_elems() * 4;
    }
    size_t pmat_per_node_elems() const {
        const size_t state_count = static_cast<size_t>(states);
        return static_cast<size_t>(rate_cats) * state_count * state_count;
    }
};

// ----- Host-side staging -----

struct HostPacking {
    std::vector<int>     postorder, preorder, parent, left, right;
    std::vector<uint8_t> is_tip;
    std::vector<fp_t>    blen;

    std::vector<int>     tip_node_ids;      // size = tips
    std::vector<uint8_t> tipchars;          // size = tips * sites

    std::vector<unsigned> site_scaler;      // size = sites or sites*rate
    std::vector<unsigned> pattern_weights;  // size = sites
    std::vector<fp_t>     pmats;
    std::vector<fp_t>     pmats_mid;        // half-branch pmats for midpoint calculations
    std::vector<fp_t>     pmats_mid_prox;   // proximal branch pmats for midpoint calculations
    std::vector<fp_t>     pmats_mid_dist;   // distal branch pmats for midpoint calculations
};

struct PlacementQueryBatch {
    size_t count = 0;
    std::vector<fp_t> branch_lengths;
    std::vector<uint8_t> query_chars;
    std::vector<fp_t> query_pmats;

    bool empty() const { return count == 0; }
    size_t size() const { return count; }
};

struct PlacementCommitContext {
    TreeBuildResult* tree = nullptr;
    HostPacking* host = nullptr;
    PlacementQueryBatch* queries = nullptr;
    PlacementOpBuffer* placement_ops = nullptr;
    const std::vector<std::string>* query_names = nullptr;
    std::vector<std::string>* inserted_query_names = nullptr;
};

struct DeviceTreeReloadTimingStats {
    double branch_copy_ms = 0.0;
    double branch_reset_ms = 0.0;
    double tipchar_copy_ms = 0.0;
    double clv_reset_ms = 0.0;
    double root_seed_ms = 0.0;
    double pmat_copy_ms = 0.0;
    double query_copy_ms = 0.0;
    double pattern_copy_ms = 0.0;
};

// ----- Host packing and upload -----
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
    int num_changed_nodes = 0
);

DeviceTree make_query_view(const DeviceTree& D, int query_idx);
void build_query_clv(
    const DeviceTree& D,
    int query_idx,
    cudaStream_t stream = 0);
void copy_unscaled_up_clv_to_query_slot(
    const DeviceTree& src,
    int src_node_id,
    DeviceTree& dst,
    int dst_query_idx,
    cudaStream_t stream = 0);
void copy_upward_state(
    const DeviceTree& src,
    DeviceTree& dst,
    cudaStream_t stream = 0);

DeviceTree upload_to_gpu(
    const TreeBuildResult& T,
    const HostPacking& H,
    const EigResult& er,
    const std::vector<double>& rate_weights,
    const std::vector<double>& rate_multipliers,
    const std::vector<double>& pi,
    size_t sites, int states, int rate_cats, bool per_rate_scaling,
    const PlacementQueryBatch* queries = nullptr,
    bool commit_to_tree = false
);
void reload_device_tree_live_data(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const PlacementQueryBatch* queries = nullptr,
    cudaStream_t stream = 0,
    DeviceTreeReloadTimingStats* timing = nullptr);
void reload_device_tree_live_data_local_spr(
    DeviceTree& D,
    const TreeBuildResult& T,
    const HostPacking& H,
    const HostPacking& base_H,
    int current_main_pmat_node,
    int& previous_main_pmat_node,
    const PlacementQueryBatch* queries = nullptr,
    cudaStream_t stream = 0,
    DeviceTreeReloadTimingStats* timing = nullptr);

struct BuildToGpuResult {
    DeviceTree dev;
    TreeBuildResult tree;
    HostPacking hostPack;
    EigResult eig;
    PlacementQueryBatch queries;
};

BuildToGpuResult BuildAllToGPU(
    const std::vector<std::string>& msa_tip_names,
    const std::vector<std::string>& msa_rows,
    const std::string& newick_text,
    const std::vector<double>& Q_rowmajor,   // size = states*states
    const std::vector<double>& pi,           // size = states
    const std::vector<double>& rate_multipliers,   // size = rate_cats
    const std::vector<double>& rate_weights, // size = rate_cats
    const std::vector<unsigned>& pattern_weights,
    size_t sites, int states, int rate_cats, bool per_rate_scaling,
    const std::vector<NewPlacementQuery>& placement_queries,
    bool commit_to_tree = false);
BuildToGpuResult BuildAllToGPU(
    const std::vector<std::string>& msa_tip_names,
    const std::vector<std::string>& msa_rows,
    const TreeBuildResult& tree,
    const std::vector<double>& Q_rowmajor,   // size = states*states
    const std::vector<double>& pi,           // size = states
    const std::vector<double>& rate_multipliers,   // size = rate_cats
    const std::vector<double>& rate_weights, // size = rate_cats
    const std::vector<unsigned>& pattern_weights,
    size_t sites, int states, int rate_cats, bool per_rate_scaling,
    const std::vector<NewPlacementQuery>& placement_queries,
    bool commit_to_tree = false);

// ----- Device lifecycle and evaluation -----

void free_device_tree(DeviceTree& D);

void launch_init_tip_clv(const DeviceTree& D, cudaStream_t stream = 0);

double eval_root_loglikelihood(
    const DeviceTree& D,
    int root_id,
    const std::vector<double>& pi,
    const std::vector<double>& rate_weights,
    cudaStream_t stream
);

void LaunchPreorderDownwardClv(
    const DeviceTree&      D,
    const TreeBuildResult& T,
    const HostPacking&     H,
    cudaStream_t           stream = 0);

void free_placement_op_buffer(
    PlacementOpBuffer& placement_ops,
    cudaStream_t stream = 0);

void DownloadClvDump(
    const DeviceTree& D,
    std::vector<fp_t>& clv_up,
    std::vector<unsigned>& scaler_up,
    cudaStream_t stream = 0);

void UploadPlacementOps(
    PlacementOpBuffer& placement_ops,
    const std::vector<NodeOpInfo>& host_ops,
    cudaStream_t stream = 0);

// Rebuild the full tree CLV state after an initial upload or global topology change.
void UpdateTreeClvs(
    DeviceTree& D,
    TreeBuildResult& T,
    HostPacking& H,
    PlacementOpBuffer& placement_ops,
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

void EvaluatePlacementQueries(
    DeviceTree& D,
    const EigResult& er,
    const std::vector<double>& rate_multipliers,
    PlacementCommitContext& commit_ctx,
    std::vector<struct PlacementResult>* placement_results_out = nullptr,
    int smoothing = 1,
    bool commit_to_tree = true,
    cudaStream_t stream = 0);

// ----- Placement query input helpers -----

std::vector<NewPlacementQuery> build_placement_query(const std::string& msa_path);
std::vector<NewPlacementQuery> build_placement_query(
    const std::vector<std::string>& msa_tip_names,
    const std::vector<std::string>& msa_rows);
