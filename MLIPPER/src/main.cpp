#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <cuda_runtime.h>

#include "io/input_validation.hpp"
#include "io/jplace.hpp"
#include "util/precision.hpp"
#include "spr/local_spr.hpp"
#include "model_utils.hpp"
#include "io/tree_newick.hpp"
#include "likelihood/root_likelihood.cuh"
#include "tree/tree.hpp"
#include "placement/placement.cuh"
#include "io/parse_file.hpp"
#include "msa_preprocess.hpp"
#include "util/mlipper_util.h"

namespace po = boost::program_options;
namespace mlenv = mlipper::env;
namespace mlinput = mlipper::input;
namespace mljplace = mlipper::jplaceio;
namespace mlmodel = mlipper::model;
namespace mltreeio = mlipper::treeio;

namespace {

// ----- Batch commit helpers -----

struct BatchCommitTimingStats {
    double free_prev_ms = 0.0;
    double build_ms = 0.0;
    double initial_update_ms = 0.0;
    double evaluate_ms = 0.0;
    double append_ms = 0.0;
    double newick_ms = 0.0;
    int batches = 0;
    int queries = 0;
    CommitTimingStats commit{};
};

using BatchCommitClock = std::chrono::steady_clock;

double batch_commit_elapsed_ms(const BatchCommitClock::time_point& start) {
    return std::chrono::duration<double, std::milli>(
        BatchCommitClock::now() - start).count();
}

void accumulate_commit_timing(CommitTimingStats& dst, const CommitTimingStats& src) {
    dst.initial_upward_host_ms += src.initial_upward_host_ms;
    dst.initial_downward_host_ms += src.initial_downward_host_ms;
    dst.initial_upward_stage_ms += src.initial_upward_stage_ms;
    dst.initial_downward_stage_ms += src.initial_downward_stage_ms;
    dst.query_reset_stage_ms += src.query_reset_stage_ms;
    dst.query_build_clv_stage_ms += src.query_build_clv_stage_ms;
    dst.query_kernel_total_ms += src.query_kernel_total_ms;
    dst.insertion_pre_clv_ms += src.insertion_pre_clv_ms;
    dst.insertion_upward_host_ms += src.insertion_upward_host_ms;
    dst.insertion_downward_host_ms += src.insertion_downward_host_ms;
    dst.insertion_upward_stage_ms += src.insertion_upward_stage_ms;
    dst.insertion_downward_stage_ms += src.insertion_downward_stage_ms;
    dst.initial_upward_ops += src.initial_upward_ops;
    dst.initial_downward_ops += src.initial_downward_ops;
    dst.insertion_upward_ops += src.insertion_upward_ops;
    dst.insertion_downward_ops += src.insertion_downward_ops;
    dst.initial_updates += src.initial_updates;
    dst.query_evals += src.query_evals;
    dst.insertion_updates += src.insertion_updates;
}

std::string cli_program_name(const char* argv0) {
    if (argv0 == nullptr || *argv0 == '\0') {
        return "MLIPPER";
    }
    const std::filesystem::path argv_path(argv0);
    const std::filesystem::path filename = argv_path.filename();
    return filename.empty() ? "MLIPPER" : filename.string();
}

po::typed_value<bool>* cli_flag(bool* target) {
    return po::value<bool>(target)->zero_tokens()->implicit_value(true);
}

std::string trim_ascii_copy(std::string value) {
    const auto is_space = [](unsigned char ch) {
        return std::isspace(ch) != 0;
    };
    value.erase(
        value.begin(),
        std::find_if_not(
            value.begin(),
            value.end(),
            is_space));
    value.erase(
        std::find_if_not(
            value.rbegin(),
            value.rend(),
            is_space).base(),
        value.end());
    return value;
}

std::vector<double> parse_cli_double_list(
    const std::vector<std::string>& raw_tokens,
    const std::string& option_name) {
    std::vector<double> values;
    for (const std::string& token : raw_tokens) {
        std::stringstream token_stream(token);
        std::string piece;
        while (std::getline(token_stream, piece, ',')) {
            const std::string trimmed = trim_ascii_copy(piece);
            if (trimmed.empty()) {
                throw mlinput::ValidationError(option_name, "contains an empty list element");
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
        throw mlinput::ValidationError(option_name, "requires at least one value");
    }
    return values;
}

int exit_with_cli_error(const std::string& program_name, const std::string& message) {
    std::cerr << program_name << ": error: " << message << "\n";
    std::cerr << "Run '" << program_name << " --help' for usage.\n";
    return 1;
}

} // namespace

int main(int argc, char** argv) {
    auto start_gpu = std::chrono::steady_clock::time_point{};
    cudaEvent_t start = nullptr;
    cudaEvent_t stop = nullptr;
    cudaStream_t stream = nullptr;
    float gpu_ms_kernel = 0.0f;
    BuildToGpuResult res{};
    PlacementOpBuffer placement_ops{};

    const std::string program_name = cli_program_name(argc > 0 ? argv[0] : nullptr);

    // Single config object filled directly by CLI options.
    parse::RunConfig config;

    // ---- Input (files/tree) ----
    std::string tree_newick;
    std::string jplace_out;
    std::string commit_tree_out;
    // Internal-only output/control knobs, not exposed via CLI.
    double commit_collapse_internal_epsilon = 1e-6;
    bool commit_to_tree = false;

    // ---- Model ----
    // Defaults match the previous `config.yaml` defaults (pre-CLI refactor).
    config.model.states = 4;
    config.model.subst_model = "GTR";
    config.model.ncat = 4;
    config.model.alpha = 0.3;
    config.model.pinv = 0.0;
    config.model.rates = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    config.model.per_rate_scaling = true;
    std::string best_model_file;
    std::vector<std::string> freqs_tokens;
    std::vector<std::string> rate_tokens;
    std::vector<std::string> rate_weight_tokens;
    bool no_per_rate_scaling = false;
    bool empirical_freqs = false;
    bool placement_fast = false;
    bool local_spr = false;
    bool local_spr_fast = false;
    bool fast_mode = false;
    int gpu_id = 0;
    bool gpu_auto = false;
    // Internal-only local SPR tuning knobs not exposed via CLI.
    int batch_insert_size = 0;
    int local_spr_radius = 4;
    int local_spr_cluster_threshold = 3;
    int local_spr_topk_per_unit = 8;
    bool local_spr_dynamic_validation_conflicts = false;
    int local_spr_rounds = 1;
    po::options_description general_options("General");
    general_options.add_options()
        ("help,h", "Show help message");

    po::options_description input_options("Input");
    input_options.add_options()
        ("tree-alignment", po::value<std::string>(&config.files.tree_alignment),
         "Reference alignment (tree MSA)")
        ("query-alignment", po::value<std::string>(&config.files.query_alignment),
         "Query alignment for placement (optional; defaults to --tree-alignment)")
        ("tree", po::value<std::string>(&config.files.tree),
         "Reference tree topology (Newick file)")
        ("tree-newick", po::value<std::string>(&tree_newick),
         "Reference tree topology (Newick string)");

    po::options_description output_options("Output");
    output_options.add_options()
        ("jplace-out", po::value<std::string>(&jplace_out),
         "Optional output path for a top-k placement jplace file")
        ("commit-to-tree", po::value<std::string>(&commit_tree_out),
         "Commit query placements to reference tree and write the final tree in Newick (.nwk) format to this path");

    po::options_description model_options("Model");
    model_options.add_options()
        ("states", po::value<int>(&config.model.states), "Number of states")
        ("subst-model", po::value<std::string>(&config.model.subst_model), "Substitution model")
        ("ncat", po::value<int>(&config.model.ncat), "Number of rate categories")
        ("alpha", po::value<double>(&config.model.alpha), "Gamma shape alpha")
        ("pinv", po::value<double>(&config.model.pinv), "Proportion of invariant sites")
        ("best-model", po::value<std::string>(&best_model_file),
         "Read a bestModel file and overwrite the corresponding model flags")
        ("freqs", po::value<std::vector<std::string>>(&freqs_tokens)->multitoken(),
         "Equilibrium freqs (comma-separated list)")
        ("empirical-freqs", cli_flag(&empirical_freqs),
         "Estimate equilibrium freqs from --tree-alignment (distributes ambiguous DNA symbols across represented states)")
        ("rates", po::value<std::vector<std::string>>(&rate_tokens)->multitoken(),
         "GTR rates rAC,rAG,rAT,rCG,rCT,rGT (comma-separated list)")
        ("rate-weights", po::value<std::vector<std::string>>(&rate_weight_tokens)->multitoken(),
         "Rate category weights (comma-separated list)")
        ("no-per-rate-scaling", cli_flag(&no_per_rate_scaling),
         "Disable per-rate scaling");

    po::options_description placement_options("Placement");
    placement_options.add_options()
        ("placement-fast", cli_flag(&placement_fast),
         "Use fast placement scoring (1 full optimization pass instead of baseline 4)")
        ("local-spr", cli_flag(&local_spr),
         "Enable local subtree SPR after each batch insert")
        ("local-spr-fast", cli_flag(&local_spr_fast),
         "Use fast local SPR scoring (1 full optimization pass instead of baseline 4)")
        ("fast", cli_flag(&fast_mode),
         "Enable both --placement-fast and --local-spr-fast")
        ("batch-insert-size", po::value<int>(&batch_insert_size),
         "Insert+commit query batches of size N (0 = all at once)")
        ("local-spr-radius", po::value<int>(&local_spr_radius),
         "Local SPR radius (filters candidate edges and defines subtree neighborhood)")
        ("local-spr-cluster-threshold", po::value<int>(&local_spr_cluster_threshold),
         "Anchor-distance threshold for grouping inserted queries into local SPR repair units")
        ("local-spr-rounds", po::value<int>(&local_spr_rounds),
         "Run up to N rounds of local SPR, rebuilding between accepted rounds")
        ("local-spr-dynamic-validation-conflicts", cli_flag(&local_spr_dynamic_validation_conflicts),
         "Validate local SPR candidates with dynamic conflict resolution instead of one-per-unit filtering");

    po::options_description runtime_options("Runtime");
    runtime_options.add_options()
        ("gpu-id", po::value<int>(&gpu_id),
         "CUDA device ordinal within the currently visible GPU set")
        ("gpu-auto", cli_flag(&gpu_auto),
         "Auto-select a visible CUDA device by acquiring an MLIPPER reservation lock");

    po::options_description all_options("MLIPPER");
    all_options.add(general_options)
        .add(input_options)
        .add(output_options)
        .add(model_options)
        .add(placement_options)
        .add(runtime_options);

    po::variables_map vm;
    try {
        po::store(
            po::command_line_parser(argc, argv)
                .options(all_options)
                .run(),
            vm);
        po::notify(vm);
    } catch (const po::error& e) {
        return exit_with_cli_error(program_name, e.what());
    }
    if (vm.count("help") > 0) {
        std::cout << all_options << "\n";
        return 0;
    }

    const bool tree_file_specified = vm.count("tree") > 0;
    const bool tree_newick_specified = vm.count("tree-newick") > 0;
    const bool freqs_specified = vm.count("freqs") > 0;
    const bool empirical_freqs_specified = vm.count("empirical-freqs") > 0;
    const bool batch_insert_size_specified = vm.count("batch-insert-size") > 0;
    const bool gpu_id_specified = vm.count("gpu-id") > 0;
    const bool jplace_out_specified = vm.count("jplace-out") > 0;
    const bool local_spr_tuning_requested =
        vm.count("local-spr-radius") > 0 ||
        vm.count("local-spr-cluster-threshold") > 0 ||
        vm.count("local-spr-rounds") > 0 ||
        vm.count("local-spr-fast") > 0;

    const std::filesystem::path config_base = std::filesystem::current_path();

    parse::RunInputs inputs;
    try {
        if (tree_file_specified && tree_newick_specified) {
            throw mlinput::ValidationError(
                "--tree-newick",
                "cannot be used together with --tree");
        }
        if (freqs_specified && empirical_freqs_specified) {
            throw mlinput::ValidationError(
                "--empirical-freqs",
                "cannot be used together with --freqs");
        }

        if (freqs_specified) {
            config.model.freqs = parse_cli_double_list(freqs_tokens, "--freqs");
        }
        if (vm.count("rates") > 0) {
            config.model.rates = parse_cli_double_list(rate_tokens, "--rates");
        }
        if (vm.count("rate-weights") > 0) {
            config.model.rate_weights = parse_cli_double_list(rate_weight_tokens, "--rate-weights");
        }

        if (batch_insert_size < 0) {
            throw mlinput::ValidationError("--batch-insert-size", "must be >= 0");
        }
        if (local_spr_radius < 0) {
            throw mlinput::ValidationError("--local-spr-radius", "must be >= 0");
        }
        if (local_spr_cluster_threshold < 0) {
            throw mlinput::ValidationError("--local-spr-cluster-threshold", "must be >= 0");
        }
        if (local_spr_rounds <= 0) {
            throw mlinput::ValidationError("--local-spr-rounds", "must be >= 1");
        }
        if (gpu_id < 0) {
            throw mlinput::ValidationError("--gpu-id", "must be >= 0");
        }
        if (gpu_auto && gpu_id_specified) {
            throw mlinput::ValidationError(
                "--gpu-auto",
                "cannot be used together with --gpu-id");
        }

        if (fast_mode) {
            placement_fast = true;
            local_spr_fast = true;
        }
        if (local_spr && batch_insert_size <= 0) {
            batch_insert_size = 5;
        }
        commit_to_tree = !commit_tree_out.empty();

        if (local_spr_tuning_requested && !local_spr) {
            throw mlinput::ValidationError(
                "--local-spr",
                "local SPR tuning flags require --local-spr");
        }
        if (local_spr && !commit_to_tree) {
            throw mlinput::ValidationError(
                "--local-spr",
                "--local-spr requires --commit-to-tree");
        }
        if (batch_insert_size_specified && batch_insert_size > 0 && !commit_to_tree) {
            throw mlinput::ValidationError(
                "--batch-insert-size",
                "batch insert mode requires --commit-to-tree");
        }
        if (batch_insert_size > 0 && commit_to_tree && jplace_out_specified) {
            throw mlinput::ValidationError(
                "--jplace-out",
                "batch insert mode does not support --jplace-out");
        }

        if (!best_model_file.empty()) {
            try {
                const auto best_model = mlmodel::parse_best_model_file(
                    mlinput::resolve_path(config_base, best_model_file));
                config.model.states = best_model.model.states;
                config.model.subst_model = best_model.model.subst_model;
                config.model.ncat = best_model.model.ncat;
                config.model.alpha = best_model.model.alpha;
                // bestModel -> pinv override is intentionally disabled for now.
                config.model.freqs = best_model.model.freqs;
                config.model.rates = best_model.model.rates;
                empirical_freqs = best_model.empirical_freqs;
            } catch (const std::exception& e) {
                throw mlinput::ValidationError("--best-model", e.what());
            }
        }
        if (no_per_rate_scaling) {
            config.model.per_rate_scaling = false;
        }
        if (config.files.tree_alignment.empty()) {
            throw mlinput::RequiredError("--tree-alignment");
        }
        if (config.files.query_alignment.empty()) {
            config.files.query_alignment = config.files.tree_alignment;
        }
        if (tree_newick.empty() && config.files.tree.empty()) {
            throw mlinput::RequiredError("one of [--tree, --tree-newick]");
        }

        if (!commit_tree_out.empty()) {
            mlinput::validate_output_path(config_base, "--commit-to-tree", commit_tree_out);
        }
        if (!jplace_out.empty()) {
            mlinput::validate_output_path(config_base, "--jplace-out", jplace_out);
        }
        if (!commit_tree_out.empty() && !jplace_out.empty()) {
            const std::filesystem::path commit_path =
                mlinput::normalize_cli_path(config_base, commit_tree_out);
            const std::filesystem::path jplace_path =
                mlinput::normalize_cli_path(config_base, jplace_out);
            if (commit_path == jplace_path) {
                throw mlinput::ValidationError(
                    "--jplace-out",
                    "must not be the same path as --commit-to-tree");
            }
        }

        mlinput::validate_model_inputs(config.model);

        parse::Alignment tree_alignment;
        try {
            tree_alignment = parse::read_alignment_file(
                mlinput::resolve_path(config_base, config.files.tree_alignment));
        } catch (const std::exception& e) {
            throw mlinput::ValidationError("--tree-alignment", e.what());
        }

        parse::Alignment query_alignment;
        try {
            query_alignment = parse::read_alignment_file(
                mlinput::resolve_path(config_base, config.files.query_alignment));
        } catch (const std::exception& e) {
            throw mlinput::ValidationError("--query-alignment", e.what());
        }

        std::string tree_text;
        if (tree_newick.empty()) {
            try {
                tree_text = mlinput::read_file_to_string(
                    mlinput::resolve_path(config_base, config.files.tree));
            } catch (const std::exception& e) {
                throw mlinput::ValidationError("--tree", e.what());
            }
        } else {
            tree_text = tree_newick;
        }
        tree_text = parse::normalize_newick(tree_text);
        mlinput::validate_newick_with_pll(
            tree_text,
            tree_newick.empty() ? "--tree" : "--tree-newick");

        inputs = parse::RunInputs{
            std::move(config),
            std::move(tree_alignment),
            std::move(query_alignment),
            std::move(tree_text)};

        if (inputs.tree_alignment.names.empty()) {
            throw mlinput::ValidationError("--tree-alignment", "contains no sequences");
        }
        if (inputs.tree_alignment.sites == 0) {
            throw mlinput::ValidationError("--tree-alignment", "contains zero sites");
        }

        if (inputs.query_alignment.names.empty()) {
            throw mlinput::ValidationError("--query-alignment", "contains no sequences");
        }
        if (inputs.query_alignment.sites == 0) {
            throw mlinput::ValidationError("--query-alignment", "contains zero sites");
        }
        if (inputs.query_alignment.sites != inputs.tree_alignment.sites) {
            throw mlinput::ValidationError(
                "--query-alignment",
                "sites mismatch with --tree-alignment (" +
                    std::to_string(inputs.query_alignment.sites) + " vs " +
                    std::to_string(inputs.tree_alignment.sites) + ")");
        }

        mlinput::validate_alignment_names(inputs.tree_alignment, "--tree-alignment");
        mlinput::validate_alignment_names(inputs.query_alignment, "--query-alignment");
        mlinput::validate_alignment_symbols(
            inputs.tree_alignment,
            inputs.config.model.states,
            "--tree-alignment");
        mlinput::validate_alignment_symbols(
            inputs.query_alignment,
            inputs.config.model.states,
            "--query-alignment");
        if (commit_to_tree) {
            mlinput::validate_query_reference_name_overlap(
                inputs.tree_alignment,
                inputs.query_alignment,
                "--query-alignment");
        }
    } catch (const mlinput::CliError& e) {
        return exit_with_cli_error(program_name, e.what());
    }

    mlipper::gpu::DeviceReservation gpu_reservation{};
    try {
    const int visible_gpu_count = mlipper::gpu::visible_device_count();
    if (visible_gpu_count <= 0) {
        throw std::runtime_error(
            "No CUDA devices are visible to MLIPPER. "
            "Check your driver/runtime setup or CUDA_VISIBLE_DEVICES.");
    }
    int requested_gpu_id = gpu_id;
    if (gpu_auto) {
        gpu_reservation =
            mlipper::gpu::reserve_any_visible_device_or_wait_or_throw();
        requested_gpu_id = gpu_reservation.device;
    } else if (gpu_id_specified) {
        gpu_reservation =
            mlipper::gpu::reserve_specific_device_or_throw(gpu_id);
        requested_gpu_id = gpu_reservation.device;
    } else if (gpu_id >= visible_gpu_count) {
        std::ostringstream oss;
        oss << "--gpu-id " << gpu_id
            << " is out of range for the current process; "
            << visible_gpu_count << " CUDA device"
            << (visible_gpu_count == 1 ? " is" : "s are")
            << " visible.";
        throw std::runtime_error(oss.str());
    }

    mlipper::gpu::set_device_or_throw(requested_gpu_id);
    const int active_gpu_id = mlipper::gpu::current_device_or_throw();
    const cudaDeviceProp active_gpu_props =
        mlipper::gpu::current_device_properties_or_throw();
    if (gpu_auto) {
        std::cout << "Auto-reserved CUDA device " << active_gpu_id
                  << " (" << active_gpu_props.name
                  << ", PCI " << gpu_reservation.bus_id << ")"
                  << " from " << visible_gpu_count << " visible GPU"
                  << (visible_gpu_count == 1 ? "" : "s") << "\n";
    } else if (gpu_id_specified) {
        std::cout << "Using reserved CUDA device " << active_gpu_id
                  << " (" << active_gpu_props.name
                  << ", PCI " << gpu_reservation.bus_id << ")"
                  << " from " << visible_gpu_count << " visible GPU"
                  << (visible_gpu_count == 1 ? "" : "s") << "\n";
    } else {
        std::cout << "Using CUDA device " << active_gpu_id
                  << " (" << active_gpu_props.name << ")"
                  << " from " << visible_gpu_count << " visible GPU"
                  << (visible_gpu_count == 1 ? "" : "s") << "\n";
    }
    start_gpu = std::chrono::steady_clock::now();
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    CUDA_CHECK(cudaStreamCreate(&stream));

    const auto& alignment = inputs.tree_alignment;

    const auto& msa_names = alignment.names;
    std::vector<std::string> rows = alignment.sequences;
    size_t sites = alignment.sites;
    const std::string& newick = inputs.tree;
    const auto& model = inputs.config.model;
    int states = model.states;
    int rate_cats = model.ncat;
    bool per_rate_scaling = model.per_rate_scaling;
    std::vector<unsigned> pattern_weights(sites, 1u);

    std::vector<double> pi =
        empirical_freqs ? mlmodel::estimate_empirical_pi(alignment, states)
                        : mlmodel::ensure_normalized_pi(model.freqs, states);
    std::vector<double> rate_weights = build_mixture_weights(model, rate_cats);
    std::vector<double> rate_multipliers = build_gamma_rate_categories(model.alpha, rate_cats);
    std::vector<double> Q = build_gtr_q_matrix(states, model, pi);

    std::cout << "Equilibrium frequencies ("
              << (empirical_freqs ? "empirical" : (model.freqs.empty() ? "uniform" : "manual"))
              << ") =";
    for (double value : pi) {
        std::cout << ' ' << std::fixed << std::setprecision(8) << value;
    }
    std::cout << "\n";

    std::vector<NewPlacementQuery> placement_queries =
        build_placement_query(inputs.query_alignment.names, inputs.query_alignment.sequences);
    if (repetitive_column_compression_enabled()) {
        remove_repetitive_columns(rows, placement_queries, pattern_weights, sites);
        if (sites == 0) {
            throw std::runtime_error("All columns were removed after repetitive-column compression.");
        }
    }
    const bool disable_pattern_weights =
        mlenv::env_flag_enabled("MLIPPER_DISABLE_PATTERN_WEIGHTS");
    const std::vector<unsigned> no_pattern_weights;
    const std::vector<unsigned>& pattern_weights_arg =
        disable_pattern_weights ? no_pattern_weights : pattern_weights;

    if (placement_fast) {
        setenv("MLIPPER_FULL_OPT_PASSES", "1", 1);
    }

    printf("Precision mode: %s\n", FP_MODE_NAME);
    std::vector<PlacementResult> placement_results;
    std::vector<std::string> committed_query_names(placement_queries.size());

    if (commit_to_tree && batch_insert_size > 0 && !placement_queries.empty()) {
        if (!jplace_out.empty()) {
            throw std::runtime_error("Batch insert mode does not support --jplace-out.");
        }
        const bool profile_batch_timing = []() {
            const char* env = std::getenv("MLIPPER_PROFILE_COMMIT_TIMING");
            return env && std::atoi(env) != 0;
        }();
        BatchCommitTimingStats batch_commit_timing;
        cudaEventRecord(start);
        std::vector<std::string> current_names = msa_names;
        std::vector<std::string> current_rows = rows;
        std::string current_tree_newick = newick;
        const int total_queries = static_cast<int>(placement_queries.size());
        placement_results.resize(placement_queries.size());

        for (int batch_start = 0; batch_start < total_queries; batch_start += batch_insert_size) {
            const int batch_end = std::min(batch_start + batch_insert_size, total_queries);
            ++batch_commit_timing.batches;
            std::vector<NewPlacementQuery> batch_queries;
            std::vector<std::string> batch_query_names;
            std::vector<int> batch_indices;
            batch_queries.reserve(batch_end - batch_start);
            batch_query_names.reserve(batch_end - batch_start);
            batch_indices.reserve(batch_end - batch_start);
            for (int idx = batch_start; idx < batch_end; ++idx) {
                batch_queries.push_back(placement_queries[(size_t)idx]);
                batch_query_names.push_back(placement_queries[(size_t)idx].msa_name);
                batch_indices.push_back(idx);
            }
            batch_commit_timing.queries += static_cast<int>(batch_queries.size());

            if (res.dev.N != 0) {
                const auto free_prev_start = BatchCommitClock::now();
                free_placement_op_buffer(placement_ops, stream);
                cudaStreamSynchronize(stream);
                free_device_tree(res.dev);
                batch_commit_timing.free_prev_ms += batch_commit_elapsed_ms(free_prev_start);
            }
            const auto build_start = BatchCommitClock::now();
            res = BuildAllToGPU(
                current_names,
                current_rows,
                current_tree_newick,
                Q,
                pi,
                rate_multipliers,
                rate_weights,
                pattern_weights_arg,
                sites,
                states,
                rate_cats,
                per_rate_scaling,
                batch_queries,
                true);
            batch_commit_timing.build_ms += batch_commit_elapsed_ms(build_start);
            if (res.tree.nodes.empty() || res.dev.N == 0) {
                throw std::runtime_error("BuildAllToGPU returned empty tree/device structures.");
            }
            if (res.tree.root_id < 0) {
                throw std::runtime_error("BuildAllToGPU produced tree with invalid root_id.");
            }

            placement_ops = PlacementOpBuffer{};
            placement_ops.profile_commit_timing = profile_batch_timing;
            const auto initial_update_start = BatchCommitClock::now();
            UpdateTreeClvs(
                res.dev,
                res.tree,
                res.hostPack,
                placement_ops,
                stream);
            batch_commit_timing.initial_update_ms +=
                batch_commit_elapsed_ms(initial_update_start);
            std::vector<PlacementResult> batch_results;
            std::vector<std::string> inserted_names(batch_queries.size());
            PlacementCommitContext batch_ctx;
            batch_ctx.tree = &res.tree;
            batch_ctx.host = &res.hostPack;
            batch_ctx.queries = &res.queries;
            batch_ctx.placement_ops = &placement_ops;
            batch_ctx.query_names = &batch_query_names;
            batch_ctx.inserted_query_names = &inserted_names;
            const auto evaluate_start = BatchCommitClock::now();
            EvaluatePlacementQueries(
                res.dev,
                res.eig,
                rate_multipliers,
                batch_ctx,
                &batch_results,
                1,
                true,
                stream);
            batch_commit_timing.evaluate_ms += batch_commit_elapsed_ms(evaluate_start);
            if (profile_batch_timing) {
                accumulate_commit_timing(batch_commit_timing.commit, placement_ops.timing);
            }

            const auto append_start = BatchCommitClock::now();
            for (size_t i = 0; i < batch_indices.size(); ++i) {
                const int qidx = batch_indices[i];
                committed_query_names[(size_t)qidx] = inserted_names[i];
                if (i < batch_results.size()) {
                    placement_results[(size_t)qidx] = batch_results[i];
                }
                current_names.push_back(inserted_names[i]);
                current_rows.push_back(placement_queries[(size_t)qidx].msa);
            }
            batch_commit_timing.append_ms += batch_commit_elapsed_ms(append_start);

            const auto newick_start = BatchCommitClock::now();
            current_tree_newick = mltreeio::write_tree_to_newick_string(res.tree);
            batch_commit_timing.newick_ms += batch_commit_elapsed_ms(newick_start);

            if (local_spr) {
                LocalSprBatchRunContext local_spr_ctx{
                    res,
                    placement_ops,
                    stream,
                    inserted_names,
                    current_names,
                    current_rows,
                    current_tree_newick,
                    pattern_weights_arg,
                    rate_weights,
                    rate_multipliers,
                    pi,
                    sites,
                    states,
                    rate_cats,
                    per_rate_scaling,
                    profile_batch_timing,
                    local_spr_fast,
                    local_spr_radius,
                    local_spr_cluster_threshold,
                    local_spr_topk_per_unit,
                    local_spr_dynamic_validation_conflicts,
                    local_spr_rounds,
                };
                run_local_spr_batch_refinement(local_spr_ctx);
            }
        }
        if (profile_batch_timing) {
            const CommitTimingStats& stats = batch_commit_timing.commit;
            std::printf(
                "Batch commit timing: batches=%d queries=%d free_prev=%.3f ms build=%.3f ms "
                "initial_update=%.3f ms evaluate=%.3f ms append=%.3f ms newick=%.3f ms\n",
                batch_commit_timing.batches,
                batch_commit_timing.queries,
                batch_commit_timing.free_prev_ms,
                batch_commit_timing.build_ms,
                batch_commit_timing.initial_update_ms,
                batch_commit_timing.evaluate_ms,
                batch_commit_timing.append_ms,
                batch_commit_timing.newick_ms);
            std::printf(
                "Batch query timing: evals=%d reset=%.3f ms build_query_clv=%.3f ms "
                "placement_kernel_total=%.3f ms\n",
                stats.query_evals,
                stats.query_reset_stage_ms,
                stats.query_build_clv_stage_ms,
                stats.query_kernel_total_ms);
            std::printf(
                "Batch insertion timing: initial_updates=%d insertion_updates=%d "
                "initial_up_host=%.3f ms initial_up_stage=%.3f ms initial_down_host=%.3f ms initial_down_stage=%.3f ms "
                "insert_pre_clv=%.3f ms insert_up_host=%.3f ms insert_up_stage=%.3f ms insert_down_host=%.3f ms insert_down_stage=%.3f ms\n",
                stats.initial_updates,
                stats.insertion_updates,
                stats.initial_upward_host_ms,
                stats.initial_upward_stage_ms,
                stats.initial_downward_host_ms,
                stats.initial_downward_stage_ms,
                stats.insertion_pre_clv_ms,
                stats.insertion_upward_host_ms,
                stats.insertion_upward_stage_ms,
                stats.insertion_downward_host_ms,
                stats.insertion_downward_stage_ms);
            std::printf(
                "Batch op summary: initial_up_ops=%lld initial_down_ops=%lld "
                "insert_up_ops=%lld insert_down_ops=%lld\n",
                stats.initial_upward_ops,
                stats.initial_downward_ops,
                stats.insertion_upward_ops,
                stats.insertion_downward_ops);
        }
    } else {
        res = BuildAllToGPU(
            msa_names,
            rows,
            newick,
            Q,
            pi,
            rate_multipliers,
            rate_weights,
            pattern_weights_arg,
            sites,
            states,
            rate_cats,
            per_rate_scaling,
            placement_queries,
            commit_to_tree);
        if (res.tree.nodes.empty() || res.dev.N == 0) {
            throw std::runtime_error("BuildAllToGPU returned empty tree/device structures.");
        }
        if (res.tree.root_id < 0) {
            throw std::runtime_error("BuildAllToGPU produced tree with invalid root_id.");
        }

        std::cout << "Uploaded. N=" << res.dev.N << ", tips=" << res.dev.tips
                    << ", per_node_elems=" << res.dev.per_node_elems() << "\n";

        if (mlenv::env_flag_enabled("MLIPPER_DEBUG_TREE_STRUCTURE")) {
            mltreeio::print_tree_structure(res.tree);
        }

        cudaEventRecord(start);
        placement_ops.profile_commit_timing = []() {
            const char* env = std::getenv("MLIPPER_PROFILE_COMMIT_TIMING");
            return env && std::atoi(env) != 0;
        }();
        UpdateTreeClvs(
            res.dev,
            res.tree,
            res.hostPack,
            placement_ops,
            stream);
        double logL = root_likelihood::compute_root_loglikelihood_total(
            res.dev,
            res.tree.root_id,
            res.dev.d_pattern_weights_u,
            nullptr,
            0.0,
            0);
        printf("Initial tree log-likelihood = %.12f\n", logL);
        std::vector<std::string> placement_query_names;
        if (commit_to_tree) {
            placement_query_names.reserve(placement_queries.size());
            for (const NewPlacementQuery& query : placement_queries) {
                placement_query_names.push_back(query.msa_name);
            }
        }
        PlacementCommitContext commit_ctx;
        commit_ctx.placement_ops = &placement_ops;

        bool actual_commit_to_tree = commit_to_tree;
        if (commit_to_tree) {
            commit_ctx.tree = &res.tree;
            commit_ctx.host = &res.hostPack;
            commit_ctx.queries = &res.queries;
            commit_ctx.query_names = &placement_query_names;
            commit_ctx.inserted_query_names = &committed_query_names;
        }
        EvaluatePlacementQueries(
            res.dev,
            res.eig,
            rate_multipliers,
            commit_ctx,
            &placement_results,
            1,
            actual_commit_to_tree,
            stream);
    }
    if (commit_to_tree) {
        const double committed_logL = root_likelihood::compute_root_loglikelihood_total(
            res.dev,
            res.tree.root_id,
            res.dev.d_pattern_weights_u,
            nullptr,
            0.0,
            0);
        printf("Committed tree log-likelihood = %.12f\n", committed_logL);
        if (placement_ops.profile_commit_timing) {
            const CommitTimingStats& stats = placement_ops.timing;
            printf(
                "Commit timing summary: initial_updates=%d insertion_updates=%d "
                "initial_up_host=%.3f ms initial_up_stage=%.3f ms initial_down_host=%.3f ms initial_down_stage=%.3f ms "
                "query_reset=%.3f ms query_build_clv=%.3f ms query_kernel_total=%.3f ms "
                "insert_pre_clv=%.3f ms insert_up_host=%.3f ms insert_up_stage=%.3f ms insert_down_host=%.3f ms insert_down_stage=%.3f ms\n",
                stats.initial_updates,
                stats.insertion_updates,
                stats.initial_upward_host_ms,
                stats.initial_upward_stage_ms,
                stats.initial_downward_host_ms,
                stats.initial_downward_stage_ms,
                stats.query_reset_stage_ms,
                stats.query_build_clv_stage_ms,
                stats.query_kernel_total_ms,
                stats.insertion_pre_clv_ms,
                stats.insertion_upward_host_ms,
                stats.insertion_upward_stage_ms,
                stats.insertion_downward_host_ms,
                stats.insertion_downward_stage_ms);
            printf(
                "Commit op summary: query_evals=%d initial_up_ops=%lld initial_down_ops=%lld "
                "insert_up_ops=%lld insert_down_ops=%lld\n",
                stats.query_evals,
                stats.initial_upward_ops,
                stats.initial_downward_ops,
                stats.insertion_upward_ops,
                stats.insertion_downward_ops);
        }

    }

    const double final_logL = root_likelihood::compute_root_loglikelihood_total(
        res.dev,
        res.tree.root_id,
        res.dev.d_pattern_weights_u,
        nullptr,
        0.0,
        0);
    printf("Final tree log-likelihood = %.12f\n", final_logL);
    if (mlenv::env_flag_enabled("MLIPPER_DEBUG_TREE_STRUCTURE")) {
        mltreeio::print_tree_structure(res.tree);
    }

    free_placement_op_buffer(placement_ops, stream);
    cudaStreamSynchronize(stream);

    if (!commit_tree_out.empty()) {
        const size_t collapsed_internal_branches = mltreeio::write_tree_to_newick_file(
            res.tree,
            commit_tree_out,
            commit_collapse_internal_epsilon);
        std::cout << "Wrote Newick tree to " << commit_tree_out;
        if (commit_collapse_internal_epsilon >= 0.0) {
            std::cout << " after collapsing " << collapsed_internal_branches
                      << " internal branches <= " << commit_collapse_internal_epsilon;
        }
        std::cout << "\n";
    }

    if (!jplace_out.empty()) {
        if (placement_results.size() != placement_queries.size()) {
            std::cerr << "placement result count mismatch before jplace export: got "
                      << placement_results.size() << ", expected " << placement_queries.size()
                      << ". Exporting available prefix only.\n";
        }

        const mljplace::JplaceTreeExport jplace_tree =
            mljplace::build_jplace_tree_export(res.tree);
        const std::vector<mljplace::JplacePlacementRecord> jplace_records =
            mljplace::build_jplace_records(
                res.tree,
                jplace_tree,
                placement_results,
                placement_queries);

        std::ostringstream invocation;
        for (int i = 0; i < argc; ++i) {
            if (i) invocation << ' ';
            invocation << argv[i];
        }
        mljplace::write_jplace(
            jplace_out,
            jplace_tree.tree,
            jplace_records,
            invocation.str());
        std::cout << "Wrote jplace to " << jplace_out << "\n";
    }

    
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    // Get elapsed time (milliseconds)
    cudaEventElapsedTime(&gpu_ms_kernel, start, stop);
    const auto end_gpu = std::chrono::steady_clock::now();
    const double gpu_ms = std::chrono::duration<double, std::milli>(end_gpu - start_gpu).count();
    printf("GPU kernel time = %.3f ms\n", gpu_ms_kernel);
    printf("GPU Wall Clock time = %.3f ms\n", gpu_ms);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    cudaStreamDestroy(stream);

    free_device_tree(res.dev);
    mlipper::gpu::release_device_reservation(&gpu_reservation);

    return 0;
    } catch (const std::exception& e) {
        free_device_tree(res.dev);
        if (start) cudaEventDestroy(start);
        if (stop) cudaEventDestroy(stop);
        if (stream) cudaStreamDestroy(stream);
        mlipper::gpu::release_device_reservation(&gpu_reservation);
        std::cout.flush();
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
