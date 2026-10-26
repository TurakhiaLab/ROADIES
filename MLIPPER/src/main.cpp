#include <filesystem>
#include <iostream>
#include <utility>

#include "io/CLI.hpp"
#include "gpu/gpu_admission.hpp"
#include "mlipper_session.hpp"
#include "util/msa_preprocess.hpp"
#include "workflow/dipper_starting_tree.hpp"

namespace {

int runDivideAndConquerWorkflow(
    int argc,
    char** argv,
    const std::filesystem::path& config_base)
{
    const mlipper::cli::DivideAndConquerParseResult parsed =
        mlipper::cli::loadDivideAndConquerConfigFromCommandLine(
            argc,
            argv,
            config_base);
    if (!parsed.should_run) {
        return parsed.exit_code;
    }

    const mlipper::cli::DivideAndConquerConfig& config = parsed.config;
    const int estimated_process_memory_mb =
        mlipper::gpu::estimate_dipper_and_divide_and_conquer_gpu_process_memory_mb(
            config.params.tip_budget,
            config.tree_alignment.sites,
            config.model.states,
            config.model.ncat,
            config.model.per_rate_scaling);
    mlipper::gpu::DeviceReservation gpu_reservation =
        mlipper::gpu::select_device_or_wait_or_throw(
            config.gpu_config,
            estimated_process_memory_mb);

    // DIPPER proposes only the starting topology. MLIPPER rebuilds aligned,
    // likelihood-ready CPU/GPU state before accepting any refinement move.
    mlipper::MlipperSession session;
    session.loadBackboneTree(
        mlipper::workflow::buildDipperStartingTree(
            config.tree_alignment,
            config.dipper_starting_tree_mode));
    const PreprocessedAlignments preprocessed =
        preprocess_alignments(config.tree_alignment);
    session.loadAlignment(preprocessed.tree_alignment);
    session.setPatternWeights(preprocessed.pattern_weights);
    session.loadModel(
        config.model,
        config.model_uses_empirical_freqs);
    session.initializeCPU();
    session.initializeDivideAndConquerGPU(
        config.params,
        std::move(gpu_reservation));

    const mlipper::DivideAndConquerNNIResult result =
        session.runDivideAndConquerNNI(
            config.nni_options);
    mlipper::cli::reportDivideAndConquerResult(result);
    const mlipper::DivideAndConquerFinalOptimizationResult final_result =
        session.runDivideAndConquerFinalOptimization(
            config.final_optimization);
    mlipper::cli::reportDivideAndConquerFinalOptimization(final_result);
    if (!config.write_tree_path.empty()) {
        session.writeTree(config.write_tree_path.string(), 0.0);
        std::cout << "Wrote session tree to "
                  << config.write_tree_path << "\n";
    }
    return 0;
}

int runSmallTipWorkflow(
    int argc,
    char** argv,
    const std::filesystem::path& config_base)
{
    const mlipper::cli::SmallTipParseResult parsed =
        mlipper::cli::loadSmallTipConfigFromCommandLine(
            argc,
            argv,
            config_base);
    if (!parsed.should_run) {
        return parsed.exit_code;
    }

    const mlipper::cli::SmallTipConfig& config = parsed.config;
    mlipper::MlipperSession session;
    session.loadBackboneTree(config.backbone_tree);
    session.loadAlignment(
        config.tree_alignment,
        config.query_alignment);
    session.setPatternWeights(config.pattern_weights);
    session.loadModel(config.model);
    session.initializeCPU();

    // Commit mode mutates and repeatedly rebuilds the resident tree. Placement-
    // only mode keeps the backbone immutable and can directly export jplace.
    if (config.params.commit_to_tree) {
        session.runSmallTipBatches(
            config.params,
            config.local_spr_params,
            config.gpu_config);
        if (config.final_optimization.optimize_model_parameters ||
            config.final_optimization.optimize_branch_lengths) {
            mlipper::cli::reportFinalOptimization(
                session.runFinalModelOptimization(
                    config.final_optimization));
        }
        // Collapse only numerically negligible internal branches in the
        // serialized small-tip tree; this does not affect optimization.
        session.writeTree(config.tree_output_path.string(), 1.0e-6);
        return 0;
    }

    session.initializeGPU(config.params, config.gpu_config);
    const mlipper::PlacementBatchResult placements =
        session.findBestLoadedPlacements();
    if (!config.jplace_output_path.empty()) {
        session.writeJplace(
            config.jplace_output_path.string(),
            config.invocation,
            placements);
        std::cout << "Wrote jplace to "
                  << config.jplace_output_path << "\n";
    }
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    const std::filesystem::path config_base =
        std::filesystem::current_path();
    try {
        if (mlipper::cli::selectWorkflow(argc, argv) ==
            mlipper::cli::WorkflowKind::DivideAndConquer) {
            return runDivideAndConquerWorkflow(argc, argv, config_base);
        }
        return runSmallTipWorkflow(argc, argv, config_base);
    } catch (const std::exception& e) {
        std::cout.flush();
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
