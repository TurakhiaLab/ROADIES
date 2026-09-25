#pragma once

#include <string>
#include <utility>
#include <vector>

#include "optimize/model_parameter_optimizer.hpp"
#include "gpu/device_buffer.hpp"
#include "io/parse_file.hpp"
#include "placement/placement.cuh"

namespace mlipper {

namespace model_optimization {

// The two workspace kinds have different dominant per-site buffers and are
// budgeted separately when an alignment must be processed in site batches.
enum class OptimizationWorkspaceKind {
    GlobalParameters,
    AllBranchGradient,
};

// Policy for choosing a site batch from current free VRAM. target_workspace is
// an upper preference, while free_memory_fraction and safety_factor retain
// headroom for resident tree state and allocator overhead.
struct SiteBatchSelectionConfig {
    size_t min_sites = 128;
    size_t alignment = 128;
    size_t target_workspace_bytes = size_t{6} << 30;
    double free_memory_fraction = 0.15;
    double safety_factor = 1.15;
};

struct SiteBatchSelection {
    size_t batch_sites = 0;
    size_t batch_count = 0;
    size_t estimated_bytes_per_site = 0;
    size_t estimated_workspace_bytes = 0;
    size_t memory_budget_bytes = 0;
};

// Return a deterministic aligned batch size and its estimated memory budget.
// free_gpu_bytes == 0 queries the active CUDA device.
SiteBatchSelection chooseOptimizationSiteBatchSize(
    size_t nodes,
    size_t tips,
    size_t total_sites,
    size_t rate_categories,
    size_t states,
    bool per_rate_scaling,
    OptimizationWorkspaceKind workspace_kind,
    size_t free_gpu_bytes = 0,
    const SiteBatchSelectionConfig& config = {});

// Owns a reusable, upward-only DeviceTree slice. Model arrays and transition
// matrices remain aliases into resident_model and therefore must outlive evaluate().
class RootSiteBatchWorkspace {
public:
    ~RootSiteBatchWorkspace();
    RootSiteBatchWorkspace(const RootSiteBatchWorkspace&) = delete;
    RootSiteBatchWorkspace& operator=(const RootSiteBatchWorkspace&) = delete;
    RootSiteBatchWorkspace() = default;
    void release();
    double evaluate(const DeviceTree& resident_model, TreeBuildResult& tree,
        HostPacking& host, size_t batch_sites, cudaStream_t stream = 0);
    double optimizeSequentialBranchNewtonSweeps(
        const DeviceTree& resident_model, TreeBuildResult& tree,
        HostPacking& host, size_t batch_sites, int sweeps,
        int newton_iterations, cudaStream_t stream = 0);

private:
    void ensureCapacity(const DeviceTree&, size_t);
    void ensureUpwardScheduleCurrent(
        const TreeBuildResult&, const HostPacking&, cudaStream_t);
    void ensureFullTreeOperationsCurrent(
        TreeBuildResult&, HostPacking&, cudaStream_t);
    DeviceTree& prepareUpwardBatch(
        const DeviceTree&, size_t begin, size_t count, cudaStream_t);
    DeviceTree batch_tree_{};
    PlacementOpBuffer prepared_upward_ops_{};
    PlacementOpBuffer prepared_downward_ops_{};
    bool upward_ops_prepared_ = false;
    bool full_tree_ops_prepared_ = false;
    std::vector<int> topology_schedule_signature_;
    size_t site_capacity_ = 0;
    std::vector<fp_t> root_down_seed_;
    gpu::DeviceBuffer<double> edge_gradient_;
    gpu::DeviceBuffer<double> edge_hessian_;
    gpu::DeviceBuffer<fp_t> branch_sumtable_;
    gpu::DeviceBuffer<double> branch_newton_state_;
    gpu::DeviceBuffer<int> branch_newton_failure_;
};

struct ModelOptimizationContext {
    // All references are non-owning and belong to the workflow/session.
    // ModelOptimizationBackend mutates them as one logical model state.
    bool& has_initialized_cpu_state;
    DeviceTree& device_tree;
    TreeBuildResult& cpu_tree;
    HostPacking& host_packing;
    PlacementOpBuffer& placement_ops;
    parse::ModelConfig& model_config;
    bool& model_uses_empirical_freqs;
    EigResult& eig;
    std::vector<double>& rate_multipliers;
    std::vector<double>& frequencies;
    RootSiteBatchWorkspace& site_batch_workspace;
    size_t& site_batch_size;
    std::string& tree_newick;
    bool site_batched = false;
    cudaStream_t stream = nullptr;
};

// Adapts the numerical model optimizer to the tree, GPU, and likelihood state
// supplied by a workflow owner. It intentionally has no MlipperSession
// dependency so it can be reused and tested independently.
class ModelOptimizationBackend {
public:
    explicit ModelOptimizationBackend(ModelOptimizationContext context)
        : context_(std::move(context)) {}

    double evaluate();
    void installSubstitutionModel(
        const std::vector<double>& frequencies,
        const std::vector<double>& rates,
        bool include_midpoint_pmats = true,
        bool recompute_resident_clvs = true,
        bool gpu_pmats_only = false);
    void installGammaAlpha(
        double alpha,
        bool include_midpoint_pmats = true,
        bool gpu_pmats_only = false);
    mlipper::BranchLengthOptimizationResult optimizeBranchLengths(
        int sweeps, int newton_iterations);
    void refreshBranchLengthTransitionMatrices(
        bool include_midpoint_pmats = true,
        bool recompute_resident_clvs = true);
    GlobalModelOptimizerCallbacks callbacks();

private:
    bool usesSiteBatching() const noexcept;
    ModelOptimizationContext context_;
};

} // namespace model_optimization
} // namespace mlipper
