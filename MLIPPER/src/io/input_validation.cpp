#include "input_validation.hpp"
#include "io/parse_file.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include <boost/algorithm/string/predicate.hpp>
#include <libpll/pll.h>

#include "util/checked_size.hpp"

namespace {

constexpr int kSupportedDnaStates = 4;
constexpr int kMaximumRateCategories = 8;

std::string preview_name_list(
    const std::vector<std::string>& names,
    size_t limit = 5)
{
    std::ostringstream oss;
    const size_t count = std::min(limit, names.size());
    for (size_t idx = 0; idx < count; ++idx) {
        if (idx) oss << ", ";
        oss << names[idx];
    }
    if (names.size() > limit) {
        oss << " ... (" << names.size() << " total)";
    }
    return oss.str();
}

void validate_positive_vector(
    const std::vector<double>& values,
    const std::string& option_name,
    const char* what)
{
    double sum = 0.0;
    for (double value : values) {
        if (!std::isfinite(value) || value <= 0.0) {
            throw mlipper::input::ValidationError(
                option_name,
                std::string(what) + " must be finite and > 0");
        }
        sum += value;
    }
    if (!(sum > 0.0) || !std::isfinite(sum)) {
        throw mlipper::input::ValidationError(
            option_name,
            std::string(what) + " must sum to a positive finite value");
    }
}

} // namespace

namespace mlipper {
namespace input {

std::filesystem::path normalize_cli_path(
    const std::filesystem::path& base,
    const std::string& raw_path)
{
    if (raw_path.empty()) return {};
    std::filesystem::path path(raw_path);
    if (path.is_relative()) path = base / path;
    if (!path.is_absolute()) {
        path = std::filesystem::absolute(path);
    }
    return path.lexically_normal();
}

void validate_output_path(
    const std::filesystem::path& base,
    const std::string& option_name,
    const std::string& raw_path)
{
    const std::filesystem::path path = normalize_cli_path(base, raw_path);
    if (path.empty() || path.filename().empty()) {
        throw ValidationError(option_name, "output path must name a file");
    }

    std::error_code ec;
    const bool output_exists = std::filesystem::exists(path, ec);
    if (ec) {
        throw ValidationError(option_name, "failed to inspect output path");
    }
    if (output_exists) {
        if (std::filesystem::is_directory(path, ec)) {
            throw ValidationError(option_name, "output path points to a directory");
        }
        if (ec) {
            throw ValidationError(option_name, "failed to inspect output path");
        }
    }

    std::filesystem::path ancestor = path.parent_path();
    while (!ancestor.empty()) {
        const bool exists = std::filesystem::exists(ancestor, ec);
        if (ec) {
            throw ValidationError(option_name, "failed to inspect output parent path");
        }
        if (exists) {
            if (!std::filesystem::is_directory(ancestor, ec) || ec) {
                throw ValidationError(
                    option_name,
                    "output parent path is not a directory: " + ancestor.string());
            }
            break;
        }
        ancestor = ancestor.parent_path();
    }
}

void validate_alignment_names(
    const parse::Alignment& alignment,
    const std::string& option_name)
{
    if (alignment.names.size() != alignment.sequences.size()) {
        throw ValidationError(option_name, "name/sequence count mismatch");
    }

    std::unordered_set<std::string> seen;
    std::vector<std::string> duplicates;
    seen.reserve(mlipper::util::checked_product(
        "alignment name set", alignment.names.size(), size_t{2}));
    for (const std::string& name : alignment.names) {
        if (name.empty()) {
            throw ValidationError(option_name, "contains an empty sequence name");
        }
        if (!seen.insert(name).second) {
            duplicates.push_back(name);
        }
    }
    if (!duplicates.empty()) {
        std::sort(duplicates.begin(), duplicates.end());
        duplicates.erase(
            std::unique(duplicates.begin(), duplicates.end()),
            duplicates.end());
        throw ValidationError(
            option_name,
            "contains duplicate sequence names: " + preview_name_list(duplicates));
    }
}

void validate_alignment_symbols(
    const parse::Alignment& alignment,
    int states,
    const std::string& option_name)
{
    if (states != kSupportedDnaStates) {
        throw ValidationError(
            option_name, "symbol validation currently supports only DNA4");
    }
    if (alignment.names.size() != alignment.sequences.size()) {
        throw ValidationError(option_name, "name/sequence count mismatch");
    }

    for (size_t seq_idx = 0; seq_idx < alignment.sequences.size(); ++seq_idx) {
        const std::string& name = alignment.names[seq_idx];
        const std::string& seq = alignment.sequences[seq_idx];
        if (seq.size() != alignment.sites) {
            throw ValidationError(
                option_name,
                "sequence '" + name + "' length does not match alignment sites");
        }
        for (size_t site_idx = 0; site_idx < seq.size(); ++site_idx) {
            const char c = seq[site_idx];
            if (pll_map_nt[static_cast<unsigned char>(c)] != 0) continue;
            std::ostringstream oss;
            oss << "sequence '" << name << "' has unsupported DNA symbol '"
                << c << "' at site " << (site_idx + 1);
            throw ValidationError(option_name, oss.str());
        }
    }
}

void validate_model_inputs(const parse::ModelConfig& model)
{
    if (model.states != kSupportedDnaStates) {
        throw ValidationError("--states", "currently only 4-state DNA input is supported");
    }
    if (model.ncat <= 0 ||
        model.ncat > kMaximumRateCategories) {
        throw ValidationError(
            "--ncat",
            "must be between 1 and " + std::to_string(
                kMaximumRateCategories));
    }
    if (!boost::algorithm::iequals(model.subst_model, "GTR")) {
        throw ValidationError("--subst-model", "currently only GTR is supported");
    }
    if (!std::isfinite(model.alpha) || model.alpha < 0.02) {
        throw ValidationError(
            "--alpha",
            "must be finite and >= 0.02 for stable discrete-Gamma categories");
    }
    if (!std::isfinite(model.pinv) || model.pinv != 0.0) {
        throw ValidationError(
            "--pinv",
            "currently only 0 is supported; invariant-site likelihoods are not wired end-to-end");
    }
    if (!model.freqs.empty()) {
        if (static_cast<int>(model.freqs.size()) != model.states) {
            throw ValidationError(
                "--freqs",
                "must have exactly " + std::to_string(model.states) + " values (states)");
        }
        validate_positive_vector(model.freqs, "--freqs", "equilibrium frequencies");
    }
    if (model.rates.size() != 6) {
        throw ValidationError("--rates", "must have exactly 6 values for 4-state GTR");
    }
    validate_positive_vector(model.rates, "--rates", "GTR rates");
}

void validate_query_reference_name_overlap(
    const parse::Alignment& tree_alignment,
    const parse::Alignment& query_alignment,
    const std::string& option_name)
{
    std::unordered_set<std::string> tree_names;
    tree_names.reserve(mlipper::util::checked_product(
        "reference name set", tree_alignment.names.size(), size_t{2}));
    for (const std::string& name : tree_alignment.names) {
        tree_names.insert(name);
    }

    std::vector<std::string> overlaps;
    for (const std::string& name : query_alignment.names) {
        if (tree_names.count(name)) {
            overlaps.push_back(name);
        }
    }
    if (!overlaps.empty()) {
        std::sort(overlaps.begin(), overlaps.end());
        overlaps.erase(
            std::unique(overlaps.begin(), overlaps.end()),
            overlaps.end());
        throw ValidationError(
            option_name,
            "query names overlap reference tip names, which would cause "
            "ambiguous committed tips: " + preview_name_list(overlaps));
    }
}

} // namespace input
} // namespace mlipper
