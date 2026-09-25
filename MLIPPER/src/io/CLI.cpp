#include "io/CLI.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <boost/algorithm/string/trim.hpp>
#include <boost/program_options.hpp>

#include "input_validation.hpp"
#include "parse_file.hpp"
#include "util/model_utils.hpp"
#include "util/msa_preprocess.hpp"

namespace mlipper::cli {

namespace po = boost::program_options;
namespace mlinput = mlipper::input;
namespace mlmodel = mlipper::model;

WorkflowKind selectWorkflow(int argc, char** argv)
{
    for (int index = 1; index < argc; ++index) {
        if (argv[index] != nullptr &&
            std::string(argv[index]) == "--divide-and-conquer") {
            return WorkflowKind::DivideAndConquer;
        }
    }
    return WorkflowKind::SmallTip;
}

namespace {

std::string cli_program_name(const char* argv0)
{
    if (argv0 == nullptr || *argv0 == '\0') {
        return "MLIPPER";
    }
    const std::filesystem::path argv_path(argv0);
    const std::filesystem::path filename = argv_path.filename();
    return filename.empty() ? "MLIPPER" : filename.string();
}

po::typed_value<bool>* cli_flag(bool* target)
{
    return po::value<bool>(target)->zero_tokens()->implicit_value(true);
}

int exit_with_cli_error(
    const std::string& program_name,
    const std::string& message)
{
    std::cerr << program_name << ": error: " << message << "\n";
    std::cerr << "Run '" << program_name << " --help' for usage.\n";
    return 1;
}

void require_int_at_least(const char* option, int value, int minimum)
{
    if (value < minimum) {
        throw mlinput::ValidationError(
            option,
            "must be >= " + std::to_string(minimum));
    }
}

parse::Alignment read_validated_alignment(
    const std::filesystem::path& config_base,
    const char* option,
    const std::string& path,
    int states)
{
    try {
        parse::Alignment alignment = parse::read_alignment_file(
            mlinput::normalize_cli_path(config_base, path).string());
        mlinput::validate_alignment_names(alignment, option);
        mlinput::validate_alignment_symbols(alignment, states, option);
        return alignment;
    } catch (const mlinput::ValidationError&) {
        throw;
    } catch (const std::exception& error) {
        throw mlinput::ValidationError(option, error.what());
    }
}

mlmodel::BestModelConfig read_best_model_or_validation_error(
    const std::filesystem::path& config_base,
    const std::string& path)
{
    try {
        return mlmodel::parse_best_model_file(
            mlinput::normalize_cli_path(config_base, path));
    } catch (const std::exception& error) {
        throw mlinput::ValidationError("--best-model", error.what());
    }
}

std::string read_file_or_validation_error(
    const std::filesystem::path& config_base,
    const char* option,
    const std::string& path)
{
    try {
        const std::filesystem::path resolved =
            mlinput::normalize_cli_path(config_base, path);
        std::ifstream input(resolved);
        if (!input) {
            throw std::runtime_error("Cannot open file: " + resolved.string());
        }
        std::ostringstream contents;
        contents << input.rdbuf();
        if (input.bad()) {
            throw std::runtime_error("Failed while reading file: " + resolved.string());
        }
        return contents.str();
    } catch (const std::exception& e) {
        throw mlinput::ValidationError(option, e.what());
    }
}

std::vector<double> parse_cli_double_list(
    const std::vector<std::string>& raw_tokens,
    const std::string& option_name)
{
    std::vector<double> values;
    for (const std::string& token : raw_tokens) {
        std::stringstream token_stream(token);
        std::string piece;
        while (std::getline(token_stream, piece, ',')) {
            const std::string trimmed = boost::algorithm::trim_copy(piece);
            if (trimmed.empty()) {
                throw mlinput::ValidationError(
                    option_name,
                    "contains an empty list element");
            }
            size_t parsed_chars = 0;
            double value = 0.0;
            try {
                value = std::stod(trimmed, &parsed_chars);
            } catch (const std::exception&) {
                throw mlinput::ValidationError(
                    option_name,
                    "invalid numeric value '" + trimmed + "'");
            }
            if (parsed_chars != trimmed.size()) {
                throw mlinput::ValidationError(
                    option_name,
                    "invalid numeric value '" + trimmed + "'");
            }
            values.push_back(value);
        }
    }
    if (values.empty()) {
        throw mlinput::ValidationError(
            option_name,
            "requires at least one value");
    }
    return values;
}

struct SmallTipCliArgs {
    parse::ModelConfig model;
    std::string tree_alignment_path;
    std::string query_alignment_path;
    std::string tree_path;
    std::string tree_newick;
    std::string jplace_out;
    std::string commit_tree_out;
    std::string best_model_file;
    std::vector<std::string> freqs_tokens;
    std::vector<std::string> rate_tokens;
    double filter_acc_lwr = 0.99;
    int gpu_id = 0;
    int batch_insert_size = 0;
    int local_spr_radius = 4;
    int local_spr_cluster_threshold = 3;
    int local_spr_rounds = 1;
    bool no_per_rate_scaling = false;
    bool empirical_freqs = false;
    bool no_local_spr = false;
    bool no_model_optimization = false;
    bool no_global_branch_optimization = false;
    bool gpu_auto = false;
};

} // namespace

// The loader is intentionally a front-loaded boundary: token parsing,
// cross-option checks, filesystem loading, model resolution, biological input
// validation, and site-pattern compression all finish before should_run=true.
// The session layer can therefore assume SmallTipConfig is execution-ready.
SmallTipParseResult loadSmallTipConfigFromCommandLine(
    int argc,
    char** argv,
    const std::filesystem::path& config_base)
{
    SmallTipParseResult result;
    SmallTipCliArgs args;
    args.model.states = 4;
    args.model.subst_model = "GTR";
    args.model.ncat = 4;
    args.model.alpha = 0.3;
    args.model.pinv = 0.0;
    args.model.rates = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    args.model.per_rate_scaling = true;

    po::options_description general_options("General");
    general_options.add_options()("help,h", "Show help message");
    po::options_description input_options("Input");
    input_options.add_options()
        ("tree-alignment", po::value<std::string>(&args.tree_alignment_path),
         "Reference alignment (tree MSA)")
        ("query-alignment", po::value<std::string>(&args.query_alignment_path),
         "Query alignment for placement (optional; defaults to --tree-alignment)")
        ("tree", po::value<std::string>(&args.tree_path),
         "Reference tree topology (Newick file)")
        ("tree-newick", po::value<std::string>(&args.tree_newick),
         "Reference tree topology (Newick string)");
    po::options_description output_options("Output");
    output_options.add_options()
        ("jplace-out", po::value<std::string>(&args.jplace_out),
         "Optional output path for an accumulated-LWR-filtered jplace file")
        ("filter-acc-lwr", po::value<double>(&args.filter_acc_lwr),
         "Accumulated LWR retained in jplace output (default: 0.99)")
        ("commit-to-tree", po::value<std::string>(&args.commit_tree_out),
         "Output path for the final tree after committing query placements");
    po::options_description model_options("Model");
    model_options.add_options()
        ("states", po::value<int>(&args.model.states),
         "Number of states (currently only 4-state DNA is supported)")
        ("subst-model", po::value<std::string>(&args.model.subst_model),
         "Substitution model (currently only GTR is supported)")
        ("ncat", po::value<int>(&args.model.ncat),
         "Number of discrete-Gamma rate categories (1-8)")
        ("alpha", po::value<double>(&args.model.alpha),
         "Discrete-Gamma shape alpha")
        ("pinv", po::value<double>(&args.model.pinv),
         "Proportion of invariant sites (unsupported; must be 0)")
        ("best-model", po::value<std::string>(&args.best_model_file),
         "Read a bestModel file and overwrite the corresponding model flags")
        ("freqs", po::value<std::vector<std::string>>(&args.freqs_tokens)->multitoken(),
         "Equilibrium freqs (comma-separated list)")
        ("empirical-freqs", cli_flag(&args.empirical_freqs),
         "Estimate equilibrium freqs from --tree-alignment")
        ("rates", po::value<std::vector<std::string>>(&args.rate_tokens)->multitoken(),
         "GTR rates rAC,rAG,rAT,rCG,rCT,rGT (comma-separated list)")
        ("no-per-rate-scaling", cli_flag(&args.no_per_rate_scaling),
         "Disable per-rate scaling");
    // Parse the removed option only so existing invocations receive a precise
    // model-support error instead of Boost's generic unknown-option message.
    po::options_description unsupported_options("Unsupported");
    unsupported_options.add_options()
        ("rate-weights",
         po::value<std::vector<std::string>>()->multitoken(),
         "Unsupported custom rate-category mixture weights");
    po::options_description placement_options("Placement");
    placement_options.add_options()
        ("no-local-spr", cli_flag(&args.no_local_spr),
         "Disable the default local subtree SPR refinement after each committed batch")
        ("batch-insert-size", po::value<int>(&args.batch_insert_size),
         "Insert+commit query batches of size N (0 = all at once)")
        ("local-spr-radius", po::value<int>(&args.local_spr_radius),
         "Local SPR radius")
        ("local-spr-cluster-threshold",
         po::value<int>(&args.local_spr_cluster_threshold),
         "Anchor-distance threshold for grouping inserted queries")
        ("local-spr-rounds", po::value<int>(&args.local_spr_rounds),
         "Run up to N rounds of local SPR");
    po::options_description optimization_options("Final optimization");
    optimization_options.add_options()
        ("no-model-optimization", cli_flag(&args.no_model_optimization),
         "Skip final model-parameter optimization")
        ("no-global-branch-optimization",
         cli_flag(&args.no_global_branch_optimization),
         "Skip final full-tree branch-length optimization");
    po::options_description runtime_options("Runtime");
    runtime_options.add_options()
        ("gpu-id", po::value<int>(&args.gpu_id),
         "CUDA device ordinal within the visible GPU set")
        ("gpu-auto", cli_flag(&args.gpu_auto),
         "Auto-select the visible CUDA device with the lowest projected reserved-memory ratio");

    po::options_description all_options("MLIPPER");
    all_options.add(general_options).add(input_options).add(output_options)
        .add(model_options).add(placement_options).add(optimization_options)
        .add(runtime_options);
    po::options_description parser_options;
    parser_options.add(all_options).add(unsupported_options);
    po::variables_map vm;
    const std::string program_name =
        cli_program_name(argc > 0 ? argv[0] : nullptr);
    try {
        // Phase one parses tokens only. All cross-option constraints and file
        // semantics are handled below before any GPU/session state is created.
        po::store(po::command_line_parser(argc, argv).options(parser_options).run(), vm);
        po::notify(vm);
    } catch (const po::error& e) {
        result.exit_code = exit_with_cli_error(program_name, e.what());
        return result;
    }
    if (vm.count("help") > 0) {
        std::cout << all_options << "\n";
        return result;
    }

    try {
        const bool tree_file_specified = vm.count("tree") > 0;
        const bool tree_newick_specified = vm.count("tree-newick") > 0;
        const bool freqs_specified = vm.count("freqs") > 0;
        const bool empirical_freqs_specified = vm.count("empirical-freqs") > 0;
        const bool batch_size_specified = vm.count("batch-insert-size") > 0;
        const bool gpu_id_specified = vm.count("gpu-id") > 0;
        const bool jplace_specified = vm.count("jplace-out") > 0;
        const bool filter_acc_lwr_specified =
            vm.count("filter-acc-lwr") > 0;
        const bool local_spr_tuning_requested =
            vm.count("local-spr-radius") > 0 ||
            vm.count("local-spr-cluster-threshold") > 0 ||
            vm.count("local-spr-rounds") > 0;

        if (tree_file_specified && tree_newick_specified) {
            throw mlinput::ValidationError(
                "--tree-newick", "cannot be used together with --tree");
        }
        if (freqs_specified && empirical_freqs_specified) {
            throw mlinput::ValidationError(
                "--empirical-freqs", "cannot be used together with --freqs");
        }
        if (freqs_specified) {
            args.model.freqs =
                parse_cli_double_list(args.freqs_tokens, "--freqs");
        }
        if (vm.count("rates") > 0) {
            args.model.rates =
                parse_cli_double_list(args.rate_tokens, "--rates");
        }
        if (vm.count("rate-weights") > 0) {
            throw mlinput::ValidationError(
                "--rate-weights",
                "mixture models are not supported; MLIPPER uses equal weights for discrete-Gamma categories");
        }
        require_int_at_least("--batch-insert-size", args.batch_insert_size, 0);
        require_int_at_least("--local-spr-radius", args.local_spr_radius, 0);
        require_int_at_least(
            "--local-spr-cluster-threshold",
            args.local_spr_cluster_threshold,
            0);
        require_int_at_least("--local-spr-rounds", args.local_spr_rounds, 1);
        require_int_at_least("--gpu-id", args.gpu_id, 0);
        if (!std::isfinite(args.filter_acc_lwr) ||
            args.filter_acc_lwr <= 0.0 ||
            args.filter_acc_lwr > 1.0) {
            throw mlinput::ValidationError(
                "--filter-acc-lwr", "must be in the interval (0, 1]");
        }
        if (filter_acc_lwr_specified && !jplace_specified) {
            throw mlinput::ValidationError(
                "--filter-acc-lwr", "requires --jplace-out");
        }
        if (args.gpu_auto && gpu_id_specified) {
            throw mlinput::ValidationError(
                "--gpu-auto", "cannot be used together with --gpu-id");
        }
        const bool commit_to_tree = !args.commit_tree_out.empty();
        if (commit_to_tree && jplace_specified) {
            throw mlinput::ValidationError(
                "--jplace-out",
                "cannot be used together with --commit-to-tree");
        }
        if (commit_to_tree && !args.no_local_spr &&
            args.batch_insert_size <= 0) {
            args.batch_insert_size = 5;
        }
        if (local_spr_tuning_requested &&
            (!commit_to_tree || args.no_local_spr)) {
            throw mlinput::ValidationError(
                "--no-local-spr",
                "local SPR tuning flags require an enabled --commit-to-tree workflow");
        }
        if (args.no_local_spr && !commit_to_tree) {
            throw mlinput::ValidationError(
                "--no-local-spr", "--no-local-spr requires --commit-to-tree");
        }
        if (!commit_to_tree &&
            (args.no_model_optimization || args.no_global_branch_optimization)) {
            throw mlinput::ValidationError(
                "--commit-to-tree",
                "final optimization flags require --commit-to-tree");
        }
        if (batch_size_specified && args.batch_insert_size > 0 &&
            !commit_to_tree) {
            throw mlinput::ValidationError(
                "--batch-insert-size",
                "batch insert mode requires --commit-to-tree");
        }
        if (!args.best_model_file.empty()) {
            const auto best_model = read_best_model_or_validation_error(
                config_base, args.best_model_file);
            args.model.states = best_model.model.states;
            args.model.subst_model = best_model.model.subst_model;
            args.model.ncat = best_model.model.ncat;
            args.model.alpha = best_model.model.alpha;
            args.model.freqs = best_model.model.freqs;
            args.model.rates = best_model.model.rates;
            args.empirical_freqs = best_model.empirical_freqs;
        }
        if (args.no_per_rate_scaling) {
            args.model.per_rate_scaling = false;
        }
        if (args.tree_alignment_path.empty()) {
            throw mlinput::CliError("--tree-alignment is required");
        }
        if (args.query_alignment_path.empty()) {
            args.query_alignment_path = args.tree_alignment_path;
        }
        if (args.tree_newick.empty() && args.tree_path.empty()) {
            throw mlinput::CliError(
                "one of [--tree, --tree-newick] is required");
        }
        if (!args.commit_tree_out.empty()) {
            mlinput::validate_output_path(
                config_base, "--commit-to-tree", args.commit_tree_out);
        }
        if (!args.jplace_out.empty()) {
            mlinput::validate_output_path(
                config_base, "--jplace-out", args.jplace_out);
        }
        if (!args.commit_tree_out.empty() && !args.jplace_out.empty()) {
            if (mlinput::normalize_cli_path(config_base, args.commit_tree_out) ==
                mlinput::normalize_cli_path(config_base, args.jplace_out)) {
                throw mlinput::ValidationError(
                    "--jplace-out",
                    "must not be the same path as --commit-to-tree");
            }
        }
        mlinput::validate_model_inputs(args.model);

        parse::Alignment tree_alignment = read_validated_alignment(
            config_base,
            "--tree-alignment",
            args.tree_alignment_path,
            args.model.states);
        parse::Alignment query_alignment =
            args.query_alignment_path == args.tree_alignment_path
            ? tree_alignment
            : read_validated_alignment(
                config_base,
                "--query-alignment",
                args.query_alignment_path,
                args.model.states);
        std::string tree_text = args.tree_newick.empty()
            ? read_file_or_validation_error(
                config_base, "--tree", args.tree_path)
            : args.tree_newick;
        if (query_alignment.sites != tree_alignment.sites) {
            throw mlinput::ValidationError(
                "--query-alignment",
                "sites mismatch with --tree-alignment (" +
                    std::to_string(query_alignment.sites) + " vs " +
                    std::to_string(tree_alignment.sites) + ")");
        }
        if (commit_to_tree) {
            mlinput::validate_query_reference_name_overlap(
                tree_alignment, query_alignment, "--query-alignment");
        }

        const std::vector<double> pi = args.empirical_freqs
            ? mlmodel::estimate_empirical_pi(
                tree_alignment, args.model.states)
            : mlmodel::ensure_normalized_pi(
                args.model.freqs, args.model.states);
        PreprocessedAlignments preprocessed = preprocess_alignments(
            tree_alignment, query_alignment);

        SmallTipConfig& output = result.config;
        output.backbone_tree = std::move(tree_text);
        output.tree_alignment = std::move(preprocessed.tree_alignment);
        output.query_alignment = std::move(preprocessed.query_alignment);
        output.model = args.model;
        output.model.freqs = pi;
        output.pattern_weights = std::move(preprocessed.pattern_weights);
        output.params.commit_to_tree = commit_to_tree;
        output.params.local_spr = commit_to_tree && !args.no_local_spr;
        output.params.accumulated_lwr_threshold =
            jplace_specified ? args.filter_acc_lwr : 0.0;
        if (output.query_alignment.names.size() >
            static_cast<size_t>(std::numeric_limits<int>::max())) {
            throw mlinput::ValidationError(
                "--query-alignment",
                "contains too many sequences for integer query indexing");
        }
        output.params.insertion_capacity = args.batch_insert_size > 0
            ? args.batch_insert_size
            : static_cast<int>(output.query_alignment.names.size());
        output.local_spr_params.local_spr_radius = args.local_spr_radius;
        output.local_spr_params.local_spr_cluster_threshold =
            args.local_spr_cluster_threshold;
        output.local_spr_params.local_spr_rounds = args.local_spr_rounds;
        if (args.gpu_auto) {
            output.gpu_config.acquire_mode =
                MlipperGpuAcquireMode::AutoAdmitAnyVisible;
        } else if (gpu_id_specified) {
            output.gpu_config.acquire_mode =
                MlipperGpuAcquireMode::AdmitSpecificDevice;
            output.gpu_config.gpu_id = args.gpu_id;
        }
        output.final_optimization.backend =
            GlobalOptimizationBackend::ResidentSequential;
        output.final_optimization.optimize_model_parameters =
            !args.no_model_optimization;
        output.final_optimization.optimize_branch_lengths =
            !args.no_global_branch_optimization;
        if (commit_to_tree) {
            output.tree_output_path = mlinput::normalize_cli_path(
                config_base, args.commit_tree_out);
        }
        if (!args.jplace_out.empty()) {
            output.jplace_output_path = mlinput::normalize_cli_path(
                config_base, args.jplace_out);
        }
        std::ostringstream invocation;
        for (int index = 0; index < argc; ++index) {
            if (index) invocation << ' ';
            if (argv[index] != nullptr) {
                invocation << argv[index];
            }
        }
        output.invocation = invocation.str();
        result.should_run = true;
        return result;
    } catch (const mlinput::CliError& e) {
        result.exit_code = exit_with_cli_error(program_name, e.what());
        return result;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        result.exit_code = 1;
        return result;
    }
}

void reportFinalOptimization(const FinalModelOptimizationResult& result)
{
    std::cout << "Resident final model/tree optimization: rounds="
              << result.rounds
              << " initial_log_likelihood=" << std::fixed
              << std::setprecision(12) << result.initial_log_likelihood
              << " final_log_likelihood=" << result.final_log_likelihood
              << " converged=" << (result.converged ? "true" : "false")
              << " model_updates=" << result.gtr_rate_updates.size()
              << " branch_updates=" << result.branch_updates.size()
              << '\n';
}

namespace {

struct DivideAndConquerCliArgs {
    std::string tree_alignment_path;
    std::string best_model_path;
    std::string write_tree_path;
    std::string dipper_starting_tree_mode = "nj-placement";
    int gpu_id = -1;
    int divide_and_conquer_core_edges = 80;
    int divide_and_conquer_tip_budget = 200;
    int divide_and_conquer_sweeps = 1;
    bool gpu_auto = false;
    bool gpu_id_specified = false;
};

bool parse_divide_and_conquer_args(
    int argc,
    char** argv,
    const std::string& program_name,
    DivideAndConquerCliArgs& out_args,
    int& exit_code)
{
    DivideAndConquerCliArgs args;
    po::options_description options("MLIPPER Divide and Conquer");
    options.add_options()
        ("help,h", "Show help message")
        ("tree-alignment", po::value<std::string>(&args.tree_alignment_path),
         "Full alignment used to build the DIPPER starting tree")
        ("best-model", po::value<std::string>(&args.best_model_path),
         "bestModel file used to configure the substitution model")
        ("write-tree", po::value<std::string>(&args.write_tree_path),
         "Optional path for the final tree after NNI and optimization")
        ("divide-and-conquer",
         "Run boundary-aware whole-tree NNI in independently optimized sectors")
        ("divide-and-conquer-core-edges",
         po::value<int>(&args.divide_and_conquer_core_edges),
         "Owned internal edges per D&C sector")
        ("divide-and-conquer-tip-budget",
         po::value<int>(&args.divide_and_conquer_tip_budget),
         "Real tips retained per D&C sector")
        ("divide-and-conquer-sweeps",
         po::value<int>(&args.divide_and_conquer_sweeps),
         "Maximum sector NNI sweeps")
        ("dipper-starting-tree-mode",
         po::value<std::string>(&args.dipper_starting_tree_mode)->default_value("nj-placement"),
         "DIPPER starting-tree algorithm: nj-placement or divide-and-conquer")
        ("gpu-id", po::value<int>(&args.gpu_id),
         "Visible CUDA device ordinal; waits for shared admission thresholds")
        ("gpu-auto", cli_flag(&args.gpu_auto),
         "Auto-select the visible CUDA device with the lowest projected reserved-memory ratio");

    po::variables_map vm;
    try {
        po::store(
            po::command_line_parser(argc, argv)
                .options(options)
                .run(),
            vm);
        po::notify(vm);
    } catch (const po::error& e) {
        exit_code = exit_with_cli_error(program_name, e.what());
        return false;
    }

    if (vm.count("help") > 0) {
        std::cout << options << "\n";
        exit_code = 0;
        return false;
    }

    args.gpu_id_specified = vm.count("gpu-id") > 0;
    if (vm.count("divide-and-conquer") == 0) {
        exit_code = exit_with_cli_error(
            program_name,
            "this interface requires --divide-and-conquer");
        return false;
    }
    out_args = std::move(args);
    return true;
}

DivideAndConquerConfig build_divide_and_conquer_config(
    const DivideAndConquerCliArgs& cli_args,
    const std::filesystem::path& config_base)
{
    if (cli_args.best_model_path.empty()) {
        throw mlipper::input::CliError("--best-model is required");
    }
    if (cli_args.tree_alignment_path.empty()) {
        throw mlipper::input::CliError("--tree-alignment is required");
    }
    if (cli_args.dipper_starting_tree_mode != "nj-placement" &&
        cli_args.dipper_starting_tree_mode != "divide-and-conquer") {
        throw mlipper::input::ValidationError(
            "--dipper-starting-tree-mode",
            "must be 'nj-placement' or 'divide-and-conquer'");
    }
    if (cli_args.gpu_id_specified && cli_args.gpu_id < 0) {
        throw mlipper::input::ValidationError("--gpu-id", "must be >= 0");
    }
    if (cli_args.gpu_auto && cli_args.gpu_id_specified) {
        throw mlipper::input::ValidationError(
            "--gpu-auto",
            "cannot be used together with --gpu-id");
    }
    if (cli_args.divide_and_conquer_core_edges <= 0 ||
        cli_args.divide_and_conquer_tip_budget <= 0 ||
        cli_args.divide_and_conquer_sweeps <= 0) {
        throw mlipper::input::ValidationError(
            "--divide-and-conquer", "has invalid scheduler limits");
    }
    DivideAndConquerConfig config;
    const auto best_model = read_best_model_or_validation_error(
        config_base, cli_args.best_model_path);
    config.model = best_model.model;
    config.model_uses_empirical_freqs = best_model.empirical_freqs;
    mlinput::validate_model_inputs(config.model);
    config.tree_alignment = read_validated_alignment(
        config_base,
        "--tree-alignment",
        cli_args.tree_alignment_path,
        config.model.states);
    config.dipper_starting_tree_mode =
        cli_args.dipper_starting_tree_mode == "nj-placement"
        ? workflow::DipperTreeMode::NJPlacement
        : workflow::DipperTreeMode::DivideAndConquer;
    config.params.tip_budget = cli_args.divide_and_conquer_tip_budget;
    config.nni_options.core_edges =
        cli_args.divide_and_conquer_core_edges;
    config.nni_options.max_sweeps =
        cli_args.divide_and_conquer_sweeps;
    config.final_optimization.core_edges =
        cli_args.divide_and_conquer_core_edges;
    if (cli_args.gpu_auto) {
        config.gpu_config.acquire_mode =
            mlipper::MlipperGpuAcquireMode::AutoAdmitAnyVisible;
    } else if (cli_args.gpu_id_specified) {
        config.gpu_config.acquire_mode =
            mlipper::MlipperGpuAcquireMode::AdmitSpecificDevice;
        config.gpu_config.gpu_id = cli_args.gpu_id;
    } else {
        config.gpu_config.acquire_mode =
            mlipper::MlipperGpuAcquireMode::UseCurrentDevice;
    }

    if (!cli_args.write_tree_path.empty()) {
        mlipper::input::validate_output_path(
            config_base,
            "--write-tree",
            cli_args.write_tree_path);
        config.write_tree_path =
            mlipper::input::normalize_cli_path(
                config_base,
                cli_args.write_tree_path);
    }

    return config;
}

} // namespace

// D&C follows the same execution-ready contract as the small-tip loader, but
// retains one combined alignment because DIPPER constructs the initial tree.
DivideAndConquerParseResult loadDivideAndConquerConfigFromCommandLine(
    int argc,
    char** argv,
    const std::filesystem::path& config_base)
{
    DivideAndConquerParseResult result;
    const std::string program_name =
        cli_program_name(argc > 0 ? argv[0] : nullptr);
    DivideAndConquerCliArgs cli_args;
    int early_exit_code = 0;
    if (!parse_divide_and_conquer_args(
            argc,
            argv,
            program_name,
            cli_args,
            early_exit_code)) {
        result.exit_code = early_exit_code;
        return result;
    }

    try {
        result.config = build_divide_and_conquer_config(cli_args, config_base);
        result.should_run = true;
        return result;
    } catch (const mlipper::input::CliError& e) {
        result.exit_code = exit_with_cli_error(program_name, e.what());
        return result;
    }
}

void reportDivideAndConquerResult(const DivideAndConquerNNIResult& result)
{
    std::cout << "MLIPPER D&C sector NNI done: sweeps=" << result.sweeps
              << " sectors=" << result.sectors
              << " covered=" << result.covered_edges
              << " internal_edges=" << result.final_internal_edges
              << " moves=" << result.accepted_moves
              << " converged=" << (result.converged ? "true" : "false")
              << "\n";
}

void reportDivideAndConquerFinalOptimization(
    const DivideAndConquerFinalOptimizationResult& result)
{
    const FinalModelOptimizationResult& optimization = result.optimization;
    std::cout << "MLIPPER D&C final optimization done: rounds="
              << optimization.rounds
              << " initial_log_likelihood=" << std::fixed
              << std::setprecision(12)
              << optimization.initial_log_likelihood
              << " final_log_likelihood="
              << optimization.final_log_likelihood
              << " converged="
              << (optimization.converged ? "true" : "false")
              << " model_updates="
              << optimization.gtr_rate_updates.size()
              << " branch_partition_sweeps="
              << optimization.branch_updates.size()
              << " branch_sectors=" << result.branch_sectors
              << " covered_branches=" << result.covered_branches
              << "\n";
}

} // namespace mlipper::cli
