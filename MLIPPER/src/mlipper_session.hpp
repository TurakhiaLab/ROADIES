#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <cuda_runtime.h>

#include "gpu/gpu_admission.hpp"
#include "io/parse_file.hpp"
#include "optimize/model_optimization_backend.hpp"
#include "placement/placement.cuh"
#include "util/mlipper_util.h"
#include "tree/divide_and_conquer.hpp"
#include "tree/tree.hpp"

struct LocalSPRSessionWorkspaceSet;

namespace mlipper {

enum class TopologyMoveType;
namespace model_optimization {
class ModelOptimizationBackend;
struct ModelOptimizationContext;
}

struct PlacementResult {
    // Stable labels of the original target-edge endpoints before insertion.
    int parent_node_label = -1;
    int child_node_label = -1;

    // Set after commit. The attachment node is the new degree-three internal
    // node connecting the parent side, child side, and inserted query tip.
    int attachment_node_label = -1;
    int query_tip_label = -1;

    double parent_side_branch_length = 0.0;
    double child_side_branch_length = 0.0;
    double pendant_branch_length = 0.0;
    double log_likelihood = 0.0;
    std::string query_name;
};

struct PlacementBatchResult {
    std::vector<PlacementResult> placements;
    std::vector<::RawPlacementResult> ranked_placements;
    std::vector<std::string> query_names;
};

struct MlipperPlacementParams {
    bool commit_to_tree = false;
    bool local_spr = false;
    // Zero uses the internal fixed-size ranking. A value in (0, 1] collects
    // enough candidates to cover this accumulated LWR for jplace output.
    double accumulated_lwr_threshold = 0.0;
    // Maximum insertions before the caller rebuilds the resident GPU tree.
    // Zero reserves space for every loaded query.
    int insertion_capacity = 0;
};

struct MlipperDivideAndConquerParams {
    int tip_budget = 1000;
};

struct MlipperLocalSPRParams {
    int local_spr_radius = 4;
    int local_spr_cluster_threshold = 3;
    int local_spr_rounds = 1;
};

struct DivideAndConquerNNIOptions {
    int core_edges = 80;
    int max_sweeps = 1;
};

struct DivideAndConquerNNIResult {
    int sweeps = 0;
    int sectors = 0;
    int accepted_moves = 0;
    size_t covered_edges = 0;
    size_t final_internal_edges = 0;
    bool converged = false;
};

struct DivideAndConquerFinalOptimizationOptions {
    int max_rounds = 10;
    int model_parameter_rounds = 100;
    int core_edges = 80;
    int branch_length_sweeps = 8;
    int branch_length_newton_iterations = 30;
    double likelihood_tolerance = 1.0e-6;
};

struct DivideAndConquerFinalOptimizationResult {
    FinalModelOptimizationResult optimization;
    size_t branch_sectors = 0;
    size_t covered_branches = 0;
};

// Owns the mutable state for one MLIPPER run. The expected lifecycle is:
//
//   load tree/alignment/model -> initializeCPU -> initialize a GPU mode
//     -> place/refine/optimize -> serialize results
//
// Topology changes invalidate derived host packing, traversal operations, and
// GPU likelihood buffers. Workflow methods keep those representations in sync;
// callers should not attempt to mix initialization modes manually.
class MlipperSession {
public:
    MlipperSession();
    ~MlipperSession();

    MlipperSession(const MlipperSession&) = delete;
    MlipperSession& operator=(const MlipperSession&) = delete;
    MlipperSession(MlipperSession&&) = delete;
    MlipperSession& operator=(MlipperSession&&) = delete;

    // Input stage. Replacing any input invalidates derived CPU state;
    // initializeCPU() and the workflow's GPU initializer must run again before
    // execution. The GPU initializer releases any allocation from the old input.
    std::string loadBackboneTree(const std::string& newick_text);
    const parse::Alignment& loadAlignment(const parse::Alignment& tree_alignment);
    const parse::Alignment& loadAlignment(
        const parse::Alignment& tree_alignment,
        const parse::Alignment& query_alignment);
    void setPatternWeights(const std::vector<unsigned>& pattern_weights);
    const parse::ModelConfig& loadModel(
        const parse::ModelConfig& model_config,
        bool model_uses_empirical_freqs = false);

    // Initialization stage. Exactly one workflow-specific GPU initializer is
    // selected after the shared CPU representation has been constructed.
    void initializeCPU();
    void initializeGPU(
        const MlipperPlacementParams& params,
        const MlipperGpuConfig& gpu_config = {});
    void initializeDivideAndConquerGPU(
        const MlipperDivideAndConquerParams& params,
        const MlipperGpuConfig& gpu_config = {});
    // Initialize on a device reserved before the session's GPU setup. The
    // reservation is passed by value and its ownership moves into this session,
    // keeping the same admission entry alive for the session lifetime. D&C uses
    // this overload after DIPPER has built its starting tree on that device.
    void initializeDivideAndConquerGPU(
        const MlipperDivideAndConquerParams& params,
        gpu::DeviceReservation reservation);

    // Execution stage. These entry points preserve synchronization among the
    // mutable CPU topology and all derived host/device representations.
    PlacementBatchResult findBestLoadedPlacements();
    void runSmallTipBatches(
        const MlipperPlacementParams& params,
        const MlipperLocalSPRParams& local_spr_params,
        const MlipperGpuConfig& gpu_config = {});
    DivideAndConquerNNIResult runDivideAndConquerNNI(
        const DivideAndConquerNNIOptions& options = {});
    DivideAndConquerFinalOptimizationResult
    runDivideAndConquerFinalOptimization(
        const DivideAndConquerFinalOptimizationOptions& options = {});
    FinalModelOptimizationResult runFinalModelOptimization(
        const FinalModelOptimizationOptions& options = {});

    // Output stage. Serialization observes session state but does not change
    // topology, likelihood caches, or GPU ownership.
    void writeTree(
        const std::string& output_path,
        double collapse_internal_epsilon = -1.0);
    void writeJplace(
        const std::string& output_path,
        const std::string& invocation,
        const PlacementBatchResult& placements) const;

private:
    enum class GpuMode {
        None,
        ResidentTree,
        SiteBatchedModel,
        DivideAndConquer,
    };

    // Internal workflow operations. Public callers use the workflow-level
    // small-tip, D&C, and final-optimization entry points above.
    struct TopologyRefinementParams {
        int radius = 0;
        int cluster_threshold = 0;
        int topk_per_unit = 0;
        int rounds = 0;
    };
    struct NNIOwnedSplit {
        size_t side_size = 0;
        std::uint64_t hash_a = 0;
        std::uint64_t hash_b = 0;
    };

    void releasePlacementWorkspace();
    void releaseDivideAndConquerSubtreeCache();
    void releaseOwnedDeviceTreeMemory();
    size_t loadedQueryCount() const noexcept;
    PlacementBatchResult findBestPlacements(
        const std::vector<SequenceRecord>& queries);
    std::vector<PlacementResult> placeAndCommitQueries(
        const std::vector<SequenceRecord>& queries);
    void runLocalSPR(
        const MlipperLocalSPRParams& params,
        const std::vector<PlacementResult>& recent_committed_placements);
    void runTopologyRefinement(
        const TopologyRefinementParams& params,
        const std::vector<divide_and_conquer::TreeEdgeEndpoints>& anchors,
        TopologyMoveType move_type,
        const std::vector<NNIOwnedSplit>& owned_nni_splits);

    std::string backbone_tree_newick_;
    parse::Alignment tree_alignment_;
    parse::Alignment query_alignment_;
    parse::ModelConfig model_config_;
    bool has_backbone_tree_ = false;
    bool has_tree_alignment_ = false;
    bool has_query_alignment_ = false;
    bool has_model_config_ = false;
    bool model_uses_empirical_freqs_ = false;
    bool has_initialized_cpu_state_ = false;
    size_t model_site_batch_size_ = 0;
    bool has_custom_pattern_weights_ = false;
    std::vector<unsigned> loaded_pattern_weights_;

    // The reservation accounts for projected cross-process use; the owned tree
    // and workspaces below hold the actual CUDA allocations.
    gpu::DeviceReservation gpu_reservation_;
    int active_gpu_device_ = -1;
    OwnedDeviceTree owned_device_tree_;
    PlacementOpBuffer placement_ops_;
    GpuMode gpu_mode_ = GpuMode::None;
    int divide_and_conquer_tip_budget_ = 1000;
    OwnedSubtreeWorkspace subtree_workspace_;
    PlacementOpBuffer subtree_placement_ops_;
    std::unique_ptr<LocalSPRSessionWorkspaceSet>
        topology_refinement_session_workspaces_;
    int last_topology_refinement_accepted_move_count_ = 0;
    PlacementQueryBatch pending_queries_;
    std::vector<std::string> pending_query_names_;
    TreeBuildResult cpu_tree_;
    HostPacking cpu_host_packing_;
    EigResult cpu_eig_;
    std::vector<double> cpu_rate_weights_;
    std::vector<double> cpu_rate_multipliers_;
    std::vector<double> cpu_pi_;
    model_optimization::RootSiteBatchWorkspace model_site_batch_workspace_;
    model_optimization::GlobalModelOptimizer global_model_optimizer_;

    bool usesModelSiteBatching() const noexcept;
    model_optimization::ModelOptimizationContext
        modelOptimizationContext();
    void enterGlobalParameterOptimizationMode(
        bool require_single_site_batch,
        GlobalOptimizationBackend backend);
    mlipper::BranchLengthOptimizationResult
    optimizeSectorPartitionedBranchLengths(
        const DivideAndConquerFinalOptimizationOptions& options,
        size_t& sector_count,
        size_t& covered_branch_count);

    int next_stable_node_label_ = 0;

    void clearInitializedCpuState();
    void validateLoadedInputsAgainstDeviceTree(
        const DeviceTree& device_tree,
        const char* context) const;
    void initializeGPUWorkspace(
        const MlipperPlacementParams& params,
        const MlipperGpuConfig& gpu_config,
        const PlacementQueryBatch* queries,
        bool prepare_site_batched_model);
    void initializeSiteBatchedModelOptimization(
        const MlipperGpuConfig& gpu_config);
    cudaStream_t currentStream() const noexcept;
    void clearPendingPlacement();
    void releasePlacementOps() noexcept;
    void releaseSubtreeWorkspace() noexcept;
    LocalSPRSessionWorkspaceSet& ensureTopologyRefinementWorkspaces();
    void releaseTopologyRefinementWorkspaces() noexcept;
    void runSectorNNI(const std::vector<NNIOwnedSplit>& owned_splits);
    void loadSubtreeWorkspace(
        TreeBuildResult& tree,
        HostPacking& host_packing,
        const PlacementQueryBatch& queries,
        const char* context,
        bool commit_to_tree = false,
        int insertion_capacity = 0);
    void ensurePlacementReady(const char* context) const;
    SequenceRecord buildLoadedQueryRecord(int query_index) const;
    void installStreamingDirectionalBoundaryMessages(
        const std::vector<divide_and_conquer::DirectionalBoundaryPort>& ports);
    HostPacking buildCurrentFullTreeHostPacking() const;
    void rebuildCpuTreeFromCurrentTopology();
    PlacementResult makePlacementResultForTarget(
        int target_id,
        double loglikelihood,
        double distal_length,
        double pendant_length,
        const std::string& tip_name) const;
    int prepareGpuDevice(
        const MlipperGpuConfig& gpu_config,
        int estimated_process_memory_mb,
        gpu::DeviceReservation* transferred_reservation = nullptr);
    void initializeDivideAndConquerGPUImpl(
        const MlipperDivideAndConquerParams& params,
        const MlipperGpuConfig& gpu_config,
        gpu::DeviceReservation* transferred_reservation);
    void releaseGpuReservation() noexcept;
    void ensureGpuResources();

};

} // namespace mlipper
