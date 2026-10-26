#include "optimize/model_optimization_backend.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "likelihood/root_loglikelihood.cuh"
#include "io/tree_newick.hpp"
#include "util/model_utils.hpp"
#include "optimize/optimization_types.hpp"
#include "placement/derivative.cuh"
#include "pmatrix/pmat_gpu.cuh"
#include "tree/tree_topology_utils.hpp"
#include "util/checked_size.hpp"
#include "util/mlipper_util.h"

namespace {

void copy_host_to_device(void* destination, const void* source, size_t bytes)
{
    CUDA_CHECK(cudaMemcpy(
        destination, source, bytes, cudaMemcpyHostToDevice));
}

void upload_transition_matrices(
    DeviceTree& device_tree, const HostPacking& host, bool include_midpoints)
{
    copy_host_to_device(
        device_tree.d_pmat, host.pmats.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            host.pmats.size(), "primary PMAT upload"));
    if (!include_midpoints) return;
    copy_host_to_device(
        device_tree.d_pmat_mid, host.pmats_mid.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            host.pmats_mid.size(), "midpoint PMAT upload"));
    const std::vector<fp_t>& proximal = host.pmats_mid_prox.empty()
        ? host.pmats_mid : host.pmats_mid_prox;
    const std::vector<fp_t>& distal = host.pmats_mid_dist.empty()
        ? host.pmats_mid : host.pmats_mid_dist;
    copy_host_to_device(
        device_tree.d_pmat_mid_prox, proximal.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            proximal.size(), "proximal midpoint PMAT upload"));
    copy_host_to_device(
        device_tree.d_pmat_mid_dist, distal.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            distal.size(), "distal midpoint PMAT upload"));
}

std::vector<int> topology_schedule_signature(
    const TreeBuildResult& tree, const HostPacking& host)
{
    std::vector<int> signature;
    if (tree.nodes.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
        tree.postorder.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
        tree.preorder.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
        host.tip_node_ids.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("Topology schedule exceeds integer indexing.");
    }
    size_t signature_capacity = mlipper::util::checked_product(
        "topology schedule signature", size_t{5}, tree.nodes.size());
    signature_capacity = mlipper::util::checked_add_size(
        signature_capacity, tree.postorder.size(), "topology schedule signature");
    signature_capacity = mlipper::util::checked_add_size(
        signature_capacity, tree.preorder.size(), "topology schedule signature");
    signature_capacity = mlipper::util::checked_add_size(
        signature_capacity, host.tip_node_ids.size(), "topology schedule signature");
    signature.reserve(mlipper::util::checked_add_size(
        signature_capacity, size_t{7}, "topology schedule signature"));
    signature.push_back(static_cast<int>(tree.nodes.size()));
    signature.push_back(tree.root_id);
    for (const TreeNode& node : tree.nodes) {
        signature.push_back(node.id);
        signature.push_back(node.parent);
        signature.push_back(node.left);
        signature.push_back(node.right);
        signature.push_back(node.is_tip ? 1 : 0);
    }
    signature.push_back(static_cast<int>(tree.postorder.size()));
    signature.insert(
        signature.end(), tree.postorder.begin(), tree.postorder.end());
    signature.push_back(static_cast<int>(tree.preorder.size()));
    signature.insert(
        signature.end(), tree.preorder.begin(), tree.preorder.end());
    signature.push_back(static_cast<int>(host.tip_node_ids.size()));
    signature.insert(
        signature.end(), host.tip_node_ids.begin(), host.tip_node_ids.end());
    return signature;
}

template <typename T>
void cuda_free_if_allocated(T*& pointer)
{
    if (pointer) cudaFree(pointer);
    pointer = nullptr;
}

void update_upward_with_cached_schedule(
    DeviceTree& D, TreeBuildResult& tree, HostPacking& host,
    PlacementOpBuffer& prepared, bool& is_prepared, cudaStream_t stream)
{
    if (!is_prepared) {
        UpdateTreeClvsUpwardOnly(D, tree, host, prepared, stream);
        is_prepared = true;
    } else {
        UpdateTreeClvsUpwardOnlyPrepared(D, prepared, stream);
    }
}

} // namespace

namespace mlipper::model_optimization {

bool ModelOptimizationBackend::usesSiteBatching() const noexcept
{
    return context_.site_batched;
}

double ModelOptimizationBackend::evaluate()
{
    if (context_.device_tree.d_tipchars == nullptr) {
        throw std::runtime_error(
            "model optimization GPU state is not initialized");
    }
    if (usesSiteBatching()) {
        return context_.site_batch_workspace.evaluate(
            context_.device_tree, context_.cpu_tree,
            context_.host_packing, context_.site_batch_size,
            context_.stream);
    }
    return mlipper::likelihood::root::compute_root_loglikelihood(
        context_.device_tree, context_.device_tree.root_id,
        nullptr, 0.0, context_.stream);
}

BranchLengthOptimizationResult
ModelOptimizationBackend::optimizeBranchLengths(
    int sweeps, int newton_iterations)
{
    if (!context_.has_initialized_cpu_state ||
        context_.device_tree.d_tipchars == nullptr) {
        throw std::runtime_error(
            "optimizeBranchLengths: CPU and GPU state must be initialized");
    }
    if (context_.model_config.pinv > 0.0) {
        throw std::runtime_error(
            "optimizeBranchLengths does not support pinv > 0");
    }
    if (sweeps <= 0 || newton_iterations <= 0) return {};

    if (!usesSiteBatching()) {
        BranchOptimizationOptions options;
        options.sweeps = sweeps;
        options.newton_iterations = newton_iterations;
        options.update_scheme = EdgeUpdateScheme::Sequential;
        options.clv_retention = ClvRetention::RebuildAll;
        options.acceptance_scope = AcceptanceScope::FullTree;
        return RunAcceptedFullTreeSequentialBranchLengthOptimization(
            context_.device_tree, context_.cpu_tree,
            context_.host_packing, context_.placement_ops,
            context_.eig, context_.rate_multipliers, nullptr,
            context_.stream, options);
    }

    BranchLengthOptimizationResult result;
    result.sweeps = sweeps;
    result.newton_iterations = newton_iterations;
    result.attempted = true;
    DeviceTree& device_tree = context_.device_tree;
    const TreeBuildResult original_tree = context_.cpu_tree;
    const HostPacking original_host = context_.host_packing;
    result.log_likelihood_before = evaluate();

    const auto start = std::chrono::steady_clock::now();
    try {
        result.log_likelihood_after =
            context_.site_batch_workspace
                .optimizeSequentialBranchNewtonSweeps(
                    device_tree, context_.cpu_tree,
                    context_.host_packing,
                    context_.site_batch_size, sweeps,
                    newton_iterations, context_.stream);
    } catch (...) {
        context_.cpu_tree = original_tree;
        context_.host_packing = original_host;
        copy_host_to_device(
            device_tree.d_blen, context_.host_packing.blen.data(),
            mlipper::util::checked_allocation_bytes<fp_t>(
                static_cast<size_t>(device_tree.N),
                "branch optimization rollback lengths"));
        refreshBranchLengthTransitionMatrices(false, false);
        throw;
    }
    result.derivative_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    result.accepted = std::isfinite(result.log_likelihood_after) &&
        result.log_likelihood_after >= result.log_likelihood_before -
            optimization::branch_lengths::kCandidateAcceptanceTolerance;
    if (result.accepted) return result;

    context_.cpu_tree = original_tree;
    context_.host_packing = original_host;
    copy_host_to_device(
        device_tree.d_blen, context_.host_packing.blen.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            static_cast<size_t>(device_tree.N),
            "branch optimization rejected lengths"));
    refreshBranchLengthTransitionMatrices(false, false);
    result.log_likelihood_after = evaluate();
    return result;
}

void ModelOptimizationBackend::installSubstitutionModel(
    const std::vector<double>& frequencies,
    const std::vector<double>& rates,
    bool include_midpoint_pmats,
    bool recompute_resident_clvs,
    bool gpu_pmats_only)
{
    DeviceTree& device_tree = context_.device_tree;
    if (device_tree.states <= 0 || device_tree.rate_cats <= 0 ||
        frequencies.size() != static_cast<size_t>(device_tree.states) ||
        context_.rate_multipliers.size() !=
            static_cast<size_t>(device_tree.rate_cats)) {
        throw std::invalid_argument(
            "installSubstitutionModel received incompatible model dimensions");
    }
    parse::ModelConfig candidate_config = context_.model_config;
    candidate_config.freqs = frequencies;
    candidate_config.rates = rates;
    const std::vector<double> q = build_gtr_q_matrix(
        device_tree.states, candidate_config, frequencies);
    EigResult candidate_eig = gtr_eigendecomp_cpu(
        q.data(), frequencies.data(), device_tree.states);
    const size_t matrix_count = mlipper::util::checked_product(
        "model eigenvector elements", device_tree.states, device_tree.states);
    const size_t lambda_count = mlipper::util::checked_product(
        "model rate eigenvalues", device_tree.rate_cats, device_tree.states);
    std::vector<fp_t> lambdas(lambda_count), eigenvectors(matrix_count),
        inverse_eigenvectors(matrix_count);
    std::vector<fp_t> frequencies_fp(frequencies.begin(), frequencies.end());
    for (int r = 0; r < device_tree.rate_cats; ++r) {
        for (int state = 0; state < device_tree.states; ++state) {
            lambdas[static_cast<size_t>(r) * device_tree.states + state] =
                static_cast<fp_t>(candidate_eig.lambdas[state] *
                    context_.rate_multipliers[r]);
        }
    }
    std::transform(
        candidate_eig.V.begin(), candidate_eig.V.end(),
        eigenvectors.begin(), [](double value) { return static_cast<fp_t>(value); });
    std::transform(
        candidate_eig.Vinv.begin(), candidate_eig.Vinv.end(),
        inverse_eigenvectors.begin(),
        [](double value) { return static_cast<fp_t>(value); });
    // Construct and validate the complete CPU candidate before publishing it.
    // Device uploads that follow mirror this committed model representation.
    context_.model_config = std::move(candidate_config);
    context_.model_uses_empirical_freqs = false;
    context_.frequencies = frequencies;
    context_.eig = std::move(candidate_eig);
    copy_host_to_device(
        device_tree.d_lambdas, lambdas.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            lambdas.size(), "model eigenvalue upload"));
    copy_host_to_device(
        device_tree.d_V, eigenvectors.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            eigenvectors.size(), "model eigenvector upload"));
    copy_host_to_device(
        device_tree.d_Vinv, inverse_eigenvectors.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            inverse_eigenvectors.size(), "model inverse-eigenvector upload"));
    copy_host_to_device(
        device_tree.d_frequencies, frequencies_fp.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            frequencies_fp.size(), "model frequency upload"));
    if (gpu_pmats_only) {
        build_all_branch_pmats_gpu(
            device_tree.N, device_tree.states, device_tree.rate_cats,
            device_tree.d_blen, device_tree.d_V, device_tree.d_Vinv,
            device_tree.d_lambdas, device_tree.d_pmat,
            context_.stream);
    } else {
        refreshBranchLengthTransitionMatrices(
            include_midpoint_pmats, recompute_resident_clvs);
    }
}

void ModelOptimizationBackend::refreshBranchLengthTransitionMatrices(
    bool include_midpoint_pmats,
    bool recompute_resident_clvs)
{
    DeviceTree& device_tree = context_.device_tree;
    const bool resident = !usesSiteBatching();
    include_midpoint_pmats |= resident && recompute_resident_clvs;
    fill_pmats_in_host_packing(
        context_.cpu_tree, context_.host_packing, context_.eig,
        context_.rate_multipliers, device_tree.states,
        device_tree.rate_cats, nullptr, 0, include_midpoint_pmats);
    upload_transition_matrices(
        device_tree, context_.host_packing, include_midpoint_pmats);
    if (resident && recompute_resident_clvs) {
        UpdateTreeClvs(
            device_tree, context_.cpu_tree, context_.host_packing,
            context_.placement_ops, context_.stream);
    }
}

void ModelOptimizationBackend::installGammaAlpha(
    double alpha, bool include_midpoint_pmats, bool gpu_pmats_only)
{
    DeviceTree& device_tree = context_.device_tree;
    std::vector<double> candidate_rate_multipliers =
        model::build_gamma_rate_categories(
        alpha, device_tree.rate_cats);
    const size_t lambda_count = mlipper::util::checked_product(
        "gamma rate eigenvalues", device_tree.rate_cats, device_tree.states);
    std::vector<fp_t> lambdas(lambda_count);
    for (int r = 0; r < device_tree.rate_cats; ++r) {
        for (int state = 0; state < device_tree.states; ++state) {
            lambdas[static_cast<size_t>(r) * device_tree.states + state] =
                static_cast<fp_t>(context_.eig.lambdas[state] *
                    candidate_rate_multipliers[r]);
        }
    }
    context_.model_config.alpha = alpha;
    context_.rate_multipliers = std::move(candidate_rate_multipliers);
    copy_host_to_device(
        device_tree.d_lambdas, lambdas.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            lambdas.size(), "gamma eigenvalue upload"));
    if (gpu_pmats_only) {
        build_all_branch_pmats_gpu(
            device_tree.N, device_tree.states, device_tree.rate_cats,
            device_tree.d_blen, device_tree.d_V, device_tree.d_Vinv,
            device_tree.d_lambdas, device_tree.d_pmat,
            context_.stream);
    } else {
        refreshBranchLengthTransitionMatrices(include_midpoint_pmats);
    }
}

GlobalModelOptimizerCallbacks ModelOptimizationBackend::callbacks()
{
    GlobalModelOptimizerCallbacks result;
    result.evaluate = [this]() { return evaluate(); };
    result.install_substitution_model =
        [this](const std::vector<double>& frequencies,
               const std::vector<double>& rates) {
            const bool site_batched = usesSiteBatching();
            installSubstitutionModel(
                frequencies, rates, false, !site_batched, site_batched);
        };
    result.install_alpha = [this](double alpha) {
        const bool site_batched = usesSiteBatching();
        installGammaAlpha(alpha, false, site_batched);
        if (!site_batched) {
            UpdateTreeClvs(
                context_.device_tree, context_.cpu_tree,
                context_.host_packing, context_.placement_ops,
                context_.stream);
        }
    };
    result.optimize_branch_lengths =
        [this](int sweeps, int newton_iterations) {
            return optimizeBranchLengths(sweeps, newton_iterations);
        };
    result.branch_lengths_accepted = [this]() {
        refreshBranchLengthTransitionMatrices(false, false);
        context_.tree_newick =
            treeio::write_tree_to_newick_string(context_.cpu_tree);
    };
    return result;
}

SiteBatchSelection chooseOptimizationSiteBatchSize(
    size_t nodes,
    size_t tips,
    size_t total_sites,
    size_t rate_categories,
    size_t states,
    bool per_rate_scaling,
    OptimizationWorkspaceKind workspace_kind,
    size_t free_gpu_bytes,
    const SiteBatchSelectionConfig& config)
{
    if (nodes == 0 || tips == 0 || total_sites == 0 ||
        rate_categories == 0 || states == 0 ||
        config.min_sites == 0 || config.alignment == 0 ||
        config.target_workspace_bytes == 0 ||
        !(config.free_memory_fraction > 0.0 &&
          config.free_memory_fraction <= 1.0) ||
        !(config.safety_factor >= 1.0)) {
        throw std::invalid_argument(
            "chooseOptimizationSiteBatchSize: invalid inputs");
    }

    const size_t clv_pools =
        workspace_kind == OptimizationWorkspaceKind::AllBranchGradient
        ? 4 : 1;
    const size_t scaler_categories =
        per_rate_scaling ? rate_categories : 1;
    const long double bytes_per_site_unscaled =
        static_cast<long double>(clv_pools) * nodes * rate_categories *
            states * sizeof(fp_t) +
        static_cast<long double>(clv_pools) * nodes * scaler_categories *
            sizeof(unsigned) +
        static_cast<long double>(tips) * sizeof(unsigned char) +
        sizeof(unsigned);
    const long double guarded_bytes_per_site = std::ceil(
        bytes_per_site_unscaled * config.safety_factor);
    if (!std::isfinite(guarded_bytes_per_site) ||
        guarded_bytes_per_site <= 0.0L ||
        guarded_bytes_per_site >
            static_cast<long double>(std::numeric_limits<size_t>::max())) {
        throw std::length_error(
            "optimization bytes-per-site estimate exceeds size_t");
    }
    const size_t bytes_per_site =
        static_cast<size_t>(guarded_bytes_per_site);

    size_t memory_budget = config.target_workspace_bytes;
    if (free_gpu_bytes > 0) {
        const size_t free_memory_budget = static_cast<size_t>(
            static_cast<long double>(free_gpu_bytes) *
            config.free_memory_fraction);
        memory_budget = std::min(memory_budget, free_memory_budget);
    }
    if (bytes_per_site == 0 || memory_budget < bytes_per_site) {
        throw std::runtime_error(
            "optimization site batching has insufficient GPU memory");
    }

    const size_t raw_batch = std::min(
        total_sites, memory_budget / bytes_per_site);
    size_t batch_sites = raw_batch;
    if (raw_batch < total_sites) {
        batch_sites =
            (raw_batch / config.alignment) * config.alignment;
    }
    const size_t minimum =
        std::min(config.min_sites, total_sites);
    if (batch_sites < minimum) {
        throw std::runtime_error(
            "optimization site batching cannot fit its minimum batch");
    }

    SiteBatchSelection selection;
    selection.batch_sites = batch_sites;
    selection.batch_count = 1 + (total_sites - 1) / batch_sites;
    selection.estimated_bytes_per_site = bytes_per_site;
    selection.estimated_workspace_bytes = mlipper::util::checked_product(
        "optimization estimated workspace", bytes_per_site, batch_sites);
    selection.memory_budget_bytes = memory_budget;
    return selection;
}

RootSiteBatchWorkspace::~RootSiteBatchWorkspace() { release(); }

void RootSiteBatchWorkspace::release() {
    int original_device = -1;
    (void)cudaGetDevice(&original_device);
    if (batch_tree_.device_id >= 0) cudaSetDevice(batch_tree_.device_id);
    cuda_free_if_allocated(batch_tree_.d_edge_midpoint_clv);
    cuda_free_if_allocated(batch_tree_.d_tipchars);
    cuda_free_if_allocated(batch_tree_.d_clv_up);
    cuda_free_if_allocated(batch_tree_.d_clv_down);
    cuda_free_if_allocated(batch_tree_.d_edge_outside_clv);
    cuda_free_if_allocated(batch_tree_.d_site_scaler_up);
    cuda_free_if_allocated(batch_tree_.d_site_scaler_down);
    cuda_free_if_allocated(batch_tree_.d_edge_midpoint_scaler);
    cuda_free_if_allocated(batch_tree_.d_edge_outside_scaler);
    cuda_free_if_allocated(batch_tree_.d_pattern_weights_u);
    edge_gradient_.reset();
    edge_hessian_.reset();
    branch_sumtable_.reset();
    branch_newton_state_.reset();
    branch_newton_failure_.reset();
    free_placement_op_buffer(prepared_upward_ops_);
    free_placement_op_buffer(prepared_downward_ops_);
    upward_ops_prepared_ = false;
    full_tree_ops_prepared_ = false;
    topology_schedule_signature_.clear();
    batch_tree_ = DeviceTree{};
    site_capacity_ = 0;
    if (original_device >= 0) {
        (void)cudaSetDevice(original_device);
    }
}

void RootSiteBatchWorkspace::ensureUpwardScheduleCurrent(
    const TreeBuildResult& tree, const HostPacking& host, cudaStream_t stream)
{
    const std::vector<int> current =
        topology_schedule_signature(tree, host);
    if (current == topology_schedule_signature_) {
        return;
    }
    if (upward_ops_prepared_) {
        free_placement_op_buffer(prepared_upward_ops_, stream);
        upward_ops_prepared_ = false;
    }
    if (full_tree_ops_prepared_) {
        free_placement_op_buffer(prepared_downward_ops_, stream);
        full_tree_ops_prepared_ = false;
    }
    topology_schedule_signature_ = current;
}

void RootSiteBatchWorkspace::ensureFullTreeOperationsCurrent(
    TreeBuildResult& tree, HostPacking& host, cudaStream_t stream)
{
    ensureUpwardScheduleCurrent(tree, host, stream);
    if (full_tree_ops_prepared_) return;
    PrepareTreeClvOperations(
        tree,
        host,
        prepared_upward_ops_,
        prepared_downward_ops_,
        stream);
    upward_ops_prepared_ = true;
    full_tree_ops_prepared_ = true;
}

void RootSiteBatchWorkspace::ensureCapacity(
    const DeviceTree& source,
    size_t sites)
{
    if (sites == 0 || source.sites == 0 || source.N <= 0 || source.tips <= 0 ||
        source.states <= 0 || source.rate_cats <= 0 || source.device_id < 0) {
        throw std::invalid_argument(
            "RootSiteBatchWorkspace received an invalid device shape");
    }
    if (site_capacity_ >= sites && batch_tree_.capacity_N >= source.N &&
        batch_tree_.capacity_tips >= source.tips &&
        batch_tree_.states == source.states &&
        batch_tree_.rate_cats == source.rate_cats &&
        batch_tree_.per_rate_scaling == source.per_rate_scaling &&
        batch_tree_.device_id == source.device_id) return;
    release();
    CUDA_CHECK(cudaSetDevice(source.device_id));
    batch_tree_.device_id = source.device_id;
    batch_tree_.capacity_N = source.N;
    batch_tree_.capacity_tips = source.tips;
    batch_tree_.states = source.states;
    batch_tree_.rate_cats = source.rate_cats;
    batch_tree_.per_rate_scaling = source.per_rate_scaling;
    site_capacity_ = sites;
    const size_t clv = mlipper::util::checked_product(
        "model batch CLV", static_cast<size_t>(source.N), sites,
        static_cast<size_t>(source.rate_cats), static_cast<size_t>(source.states));
    const size_t scaler = mlipper::util::checked_product(
        "model batch scaler", static_cast<size_t>(source.N), sites,
        static_cast<size_t>(source.per_rate_scaling ? source.rate_cats : 1));
    CUDA_CHECK(cudaMalloc(&batch_tree_.d_tipchars,
        mlipper::util::checked_product("model batch tip characters",
            static_cast<size_t>(source.tips), sites)));
    CUDA_CHECK(cudaMalloc(&batch_tree_.d_clv_up,
        mlipper::util::checked_allocation_bytes<fp_t>(clv, "model batch CLV")));
    CUDA_CHECK(cudaMalloc(&batch_tree_.d_site_scaler_up,
        mlipper::util::checked_allocation_bytes<unsigned>(scaler, "model batch scaler")));
    CUDA_CHECK(cudaMalloc(&batch_tree_.d_pattern_weights_u,
        mlipper::util::checked_allocation_bytes<unsigned>(sites, "model batch pattern weights")));
}

DeviceTree& RootSiteBatchWorkspace::prepareUpwardBatch(
    const DeviceTree& source, size_t begin, size_t count, cudaStream_t stream)
{
    if (count == 0 || begin > source.sites || count > source.sites - begin ||
        count > site_capacity_) {
        throw std::out_of_range("Invalid root likelihood site batch range");
    }
    DeviceTree& D = batch_tree_;
    D.N = source.N; D.tips = source.tips; D.inners = source.inners;
    D.root_id = source.root_id; D.sites = count; D.log2_stride = source.log2_stride;
    D.d_tip_node_ids = source.d_tip_node_ids; D.d_tipmap = source.d_tipmap;
    D.d_pmat = source.d_pmat; D.d_frequencies = source.d_frequencies;
    D.d_rate_weights = source.d_rate_weights;
    CUDA_CHECK(cudaMemcpy2DAsync(D.d_tipchars, count,
        source.d_tipchars + begin, source.sites, count, source.tips,
        cudaMemcpyDeviceToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(D.d_pattern_weights_u,
        source.d_pattern_weights_u + begin,
        mlipper::util::checked_allocation_bytes<unsigned>(
            count, "model batch pattern-weight copy"),
        cudaMemcpyDeviceToDevice, stream));
    const size_t clv_count = mlipper::util::checked_product(
        "model batch CLV clear", source.N, count,
        source.rate_cats, source.states);
    const size_t scaler_count = mlipper::util::checked_product(
        "model batch scaler clear", source.N, count,
        source.per_rate_scaling ? source.rate_cats : 1);
    CUDA_CHECK(cudaMemsetAsync(D.d_clv_up, 0,
        mlipper::util::checked_allocation_bytes<fp_t>(
            clv_count, "model batch CLV clear"), stream));
    CUDA_CHECK(cudaMemsetAsync(D.d_site_scaler_up, 0,
        mlipper::util::checked_allocation_bytes<unsigned>(
            scaler_count, "model batch scaler clear"), stream));
    return D;
}

double RootSiteBatchWorkspace::evaluate(const DeviceTree& source, TreeBuildResult& tree,
    HostPacking& host, size_t batch_sites, cudaStream_t stream) {
    if (batch_sites == 0) {
        throw std::invalid_argument("Root likelihood batch size must be positive");
    }
    ensureCapacity(source, batch_sites);
    ensureUpwardScheduleCurrent(tree, host, stream);
    double total = 0.0;
    for (size_t begin = 0; begin < source.sites; begin += batch_sites) {
        const size_t count = std::min(batch_sites, source.sites - begin);
        DeviceTree& D = prepareUpwardBatch(source, begin, count, stream);
        update_upward_with_cached_schedule(
            D, tree, host, prepared_upward_ops_, upward_ops_prepared_, stream);
        total += mlipper::likelihood::root::compute_root_loglikelihood(
            D, D.root_id, nullptr, 0.0, stream);
    }
    return total;
}

double RootSiteBatchWorkspace::optimizeSequentialBranchNewtonSweeps(
    const DeviceTree& source, TreeBuildResult& tree, HostPacking& host,
    size_t batch_sites, int sweeps, int newton_iterations,
    cudaStream_t stream)
{
    if (source.sites > batch_sites) {
        throw std::runtime_error(
            "sequential sumtable BLO currently requires one site batch");
    }
    if (sweeps <= 0 || newton_iterations <= 0) {
        return evaluate(source, tree, host, batch_sites, stream);
    }

    ensureCapacity(source, batch_sites);
    ensureFullTreeOperationsCurrent(tree, host, stream);
    if (tree.root_id < 0 || tree.root_id >= source.N ||
        tree.nodes.size() < static_cast<size_t>(source.N) ||
        host.blen.size() < static_cast<size_t>(source.N)) {
        throw std::invalid_argument(
            "Sequential BLO tree shape does not match the device tree");
    }
    DeviceTree& D = prepareUpwardBatch(source, 0, source.sites, stream);
    D.downward_pmat_indexing = DownwardPmatIndexing::Rows;
    D.d_pmat_mid = source.d_pmat_mid;
    D.d_pmat_mid_prox = source.d_pmat_mid_prox;
    D.d_pmat_mid_dist = source.d_pmat_mid_dist;
    D.d_lambdas = source.d_lambdas;
    D.d_V = source.d_V;
    D.d_Vinv = source.d_Vinv;
    D.d_blen = source.d_blen;

    // A degree-2 root splits one unrooted edge into two representational
    // branches. Applying the ordinary lower bound to both halves would make
    // the effective minimum edge length twice as large. Keep one artificial
    // half at zero and optimize the combined length on the other half.
    const int artificial_root_child =
        tree.nodes[static_cast<size_t>(tree.root_id)].left;
    const int optimized_root_child =
        tree.nodes[static_cast<size_t>(tree.root_id)].right;
    if (artificial_root_child < 0 || optimized_root_child < 0 ||
        artificial_root_child >= source.N || optimized_root_child >= source.N) {
        throw std::runtime_error(
            "sequential BLO requires a bifurcating root");
    }
    const fp_t combined_root_length = static_cast<fp_t>(
        static_cast<double>(tree.nodes[static_cast<size_t>(artificial_root_child)]
                                .branch_length_to_parent) +
        static_cast<double>(tree.nodes[static_cast<size_t>(optimized_root_child)]
                                .branch_length_to_parent));
    tree.nodes[static_cast<size_t>(artificial_root_child)]
        .branch_length_to_parent = fp_t(0);
    tree.nodes[static_cast<size_t>(optimized_root_child)]
        .branch_length_to_parent = combined_root_length;
    host.blen[static_cast<size_t>(artificial_root_child)] = fp_t(0);
    host.blen[static_cast<size_t>(optimized_root_child)] = combined_root_length;
    CUDA_CHECK(cudaMemcpyAsync(
        D.d_blen + artificial_root_child,
        &host.blen[static_cast<size_t>(artificial_root_child)], sizeof(fp_t),
        cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(
        D.d_blen + optimized_root_child,
        &host.blen[static_cast<size_t>(optimized_root_child)], sizeof(fp_t),
        cudaMemcpyHostToDevice, stream));
    build_single_branch_pmat_gpu(
        artificial_root_child, D.states, D.rate_cats, D.d_blen,
        D.d_V, D.d_Vinv, D.d_lambdas, D.d_pmat, stream);
    build_single_branch_pmat_gpu(
        optimized_root_child, D.states, D.rate_cats, D.d_blen,
        D.d_V, D.d_Vinv, D.d_lambdas, D.d_pmat, stream);

    const size_t clv = mlipper::util::checked_product(
        "branch optimization CLV", static_cast<size_t>(source.N), source.sites,
        static_cast<size_t>(source.rate_cats), static_cast<size_t>(source.states));
    const size_t scaler = mlipper::util::checked_product(
        "branch optimization scaler", static_cast<size_t>(source.N), source.sites,
        static_cast<size_t>(source.per_rate_scaling ? source.rate_cats : 1));
    if (!D.d_clv_down) {
        const size_t clv_bytes = mlipper::util::checked_allocation_bytes<fp_t>(
            clv, "branch optimization CLV");
        const size_t scaler_bytes =
            mlipper::util::checked_allocation_bytes<unsigned>(
                scaler, "branch optimization scaler");
        CUDA_CHECK(cudaMalloc(&D.d_clv_down, clv_bytes));
        CUDA_CHECK(cudaMalloc(&D.d_edge_midpoint_clv, clv_bytes));
        CUDA_CHECK(cudaMalloc(&D.d_edge_outside_clv, clv_bytes));
        CUDA_CHECK(cudaMalloc(
            &D.d_site_scaler_down, scaler_bytes));
        CUDA_CHECK(cudaMalloc(
            &D.d_edge_midpoint_scaler, scaler_bytes));
        CUDA_CHECK(cudaMalloc(
            &D.d_edge_outside_scaler, scaler_bytes));
    }
    branch_sumtable_.ensureCapacity(
        mlipper::util::checked_product(
            "branch optimization sumtable", source.sites,
            static_cast<size_t>(source.rate_cats),
            static_cast<size_t>(source.states)));
    edge_gradient_.ensureCapacity(source.N);
    edge_hessian_.ensureCapacity(source.N);
    branch_newton_state_.ensureCapacity(6);
    branch_newton_failure_.ensureCapacity(1);
    CUDA_CHECK(cudaMemsetAsync(
        branch_newton_failure_.get(), 0, sizeof(int), stream));

    const size_t clv_bytes = mlipper::util::checked_allocation_bytes<fp_t>(
        clv, "branch optimization CLV clear");
    const size_t scaler_bytes =
        mlipper::util::checked_allocation_bytes<unsigned>(
            scaler, "branch optimization scaler clear");
    CUDA_CHECK(cudaMemsetAsync(D.d_clv_down, 0, clv_bytes, stream));
    CUDA_CHECK(cudaMemsetAsync(
        D.d_edge_midpoint_clv, 0, clv_bytes, stream));
    CUDA_CHECK(cudaMemsetAsync(
        D.d_edge_outside_clv, 0, clv_bytes, stream));
    CUDA_CHECK(cudaMemsetAsync(
        D.d_site_scaler_down, 0, scaler_bytes, stream));
    CUDA_CHECK(cudaMemsetAsync(
        D.d_edge_midpoint_scaler, 0, scaler_bytes, stream));
    CUDA_CHECK(cudaMemsetAsync(
        D.d_edge_outside_scaler, 0, scaler_bytes, stream));
    const size_t root_count = mlipper::util::checked_product(
        "branch optimization root seed", source.sites,
        source.rate_cats, source.states);
    root_down_seed_.assign(root_count, fp_t(1));
    CUDA_CHECK(cudaMemcpyAsync(
        D.d_clv_down +
            static_cast<size_t>(D.root_id) * D.per_node_elems(),
        root_down_seed_.data(),
        mlipper::util::checked_allocation_bytes<fp_t>(
            root_count, "branch optimization root seed upload"),
        cudaMemcpyHostToDevice, stream));
    UpdateTreeClvsPrepared(
        D, prepared_upward_ops_, prepared_downward_ops_, stream);

    std::vector<int> upward_index(static_cast<size_t>(source.N), -1);
    if (prepared_upward_ops_.upward_ops_host.size() >
            static_cast<size_t>(std::numeric_limits<int>::max()) ||
        prepared_downward_ops_.downward_ops_host.size() >
            static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error(
            "Branch optimization schedule exceeds integer indexing");
    }
    for (size_t index = 0;
         index < prepared_upward_ops_.upward_ops_host.size(); ++index) {
        const int node_id =
            prepared_upward_ops_.upward_ops_host[index].parent_id;
        if (node_id >= 0 && node_id < source.N) {
            upward_index[static_cast<size_t>(node_id)] =
                static_cast<int>(index);
        }
    }
    std::vector<int> downward_index(static_cast<size_t>(source.N), -1);
    for (size_t index = 0;
         index < prepared_downward_ops_.downward_ops_host.size(); ++index) {
        const NodeOpInfo& op =
            prepared_downward_ops_.downward_ops_host[index];
        const bool target_is_left =
            op.dir_tag == static_cast<uint8_t>(CLV_DIR_DOWN_LEFT);
        const int target_id = target_is_left ? op.left_id : op.right_id;
        if (target_id >= 0 && target_id < source.N) {
            downward_index[static_cast<size_t>(target_id)] =
                static_cast<int>(index);
        }
    }

    const std::vector<SequentialBranchTraversalStep> traversal =
        build_sequential_branch_traversal_steps(tree);
    for (int sweep = 0; sweep < sweeps; ++sweep) {
        if (sweep > 0) {
            UpdateTreeClvsDownwardOnlyPrepared(
                D, prepared_downward_ops_, stream);
        }
        for (const SequentialBranchTraversalStep& step : traversal) {
            if (step.kind ==
                SequentialBranchTraversalStepKind::RefreshParentUpward) {
                const int up_op =
                    upward_index[static_cast<size_t>(step.node_id)];
                if (up_op >= 0) {
                    UpdateSingleTreeClvUpwardWarpSitePrepared(
                        D, prepared_upward_ops_, up_op, stream);
                }
                continue;
            }
            const int child_id = step.node_id;
            const int down_op =
                downward_index[static_cast<size_t>(child_id)];
            if (down_op < 0) {
                throw std::runtime_error(
                    "sequential BLO is missing a downward operation");
            }
            // Refresh the parent-side message from all previously accepted
            // coordinates before optimizing this edge.
            BuildSingleTreeEdgeOutsideWarpSitePrepared(
                D, prepared_downward_ops_, down_op, stream);
            if (child_id != artificial_root_child) {
                OptimizeSingleTreeEdgeFromCurrentClvsWarpSite(
                    D, child_id, branch_sumtable_.get(), edge_gradient_.get(),
                    edge_hessian_.get(), branch_newton_state_.get(),
                    branch_newton_failure_.get(), newton_iterations, stream);
                build_single_branch_pmat_gpu(
                    child_id, D.states, D.rate_cats, D.d_blen,
                    D.d_V, D.d_Vinv, D.d_lambdas, D.d_pmat, stream);
            }
            if (!tree.nodes[static_cast<size_t>(child_id)].is_tip) {
                // The parent-side edge-outside message is independent of this
                // edge length and was just computed above. Reapply only the
                // accepted target PMAT; tips have no descendants that consume
                // a downward message, so they need no post-update refresh.
                RefreshSingleTreeChildDownPrepared(
                    D,
                    prepared_downward_ops_,
                    down_op,
                    child_id,
                    stream);
            }
        }
    }

    // Rebuild once for the final score. Besides matching RAxML's end-of-sweep
    // likelihood check, this prevents an incremental-message defect from
    // being mistaken for an optimizer improvement during bring-up.
    UpdateTreeClvsPrepared(
        D, prepared_upward_ops_, prepared_downward_ops_, stream);
    const double final_log_likelihood =
        mlipper::likelihood::root::compute_root_loglikelihood(
            D, D.root_id, nullptr, 0.0, stream);
    host.blen.resize(static_cast<size_t>(source.N));
    CUDA_CHECK(cudaMemcpyAsync(
        host.blen.data(), D.d_blen,
        mlipper::util::checked_allocation_bytes<fp_t>(
            static_cast<size_t>(source.N),
            "branch optimization length download"),
        cudaMemcpyDeviceToHost, stream));
    int first_newton_failure = 0;
    CUDA_CHECK(cudaMemcpyAsync(
        &first_newton_failure, branch_newton_failure_.get(), sizeof(int),
        cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    if (first_newton_failure != 0) {
        throw std::runtime_error(
            "sequential BLO encountered non-finite Newton state at node " +
            std::to_string(first_newton_failure - 1));
    }
    for (int node_id = 0; node_id < source.N; ++node_id) {
        tree.nodes[static_cast<size_t>(node_id)].branch_length_to_parent =
            host.blen[static_cast<size_t>(node_id)];
    }
    return final_log_likelihood;
}

} // namespace mlipper::model_optimization
