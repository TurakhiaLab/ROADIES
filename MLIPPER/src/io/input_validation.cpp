#include "input_validation.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <unistd.h>

#include <libpll/pll.h>

namespace {

std::string uppercase_ascii(std::string value) {
    for (char& ch : value) {
        ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    }
    return value;
}

std::string preview_name_list(const std::vector<std::string>& names, size_t limit = 5) {
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

bool is_valid_dna4_alignment_char(char c) {
    switch (static_cast<unsigned char>(std::toupper(static_cast<unsigned char>(c)))) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'U':
        case 'R':
        case 'Y':
        case 'S':
        case 'W':
        case 'K':
        case 'M':
        case 'B':
        case 'D':
        case 'H':
        case 'V':
        case 'N':
        case '-':
        case '.':
        case '?':
            return true;
        default:
            return false;
    }
}

void validate_positive_vector(
    const std::vector<double>& values,
    const std::string& option_name,
    const char* what) {
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

std::string read_file_to_string(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open file: " + path);
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

void validate_newick_with_pll(
    const std::string& tree_text,
    const std::string& option_name) {
    if (tree_text.empty()) {
        throw ValidationError(option_name, "tree is empty");
    }

    std::filesystem::path tree_path_template =
        std::filesystem::temp_directory_path() / "mlipper-validate-XXXXXX.nwk";
    std::string tree_path_string = tree_path_template.string();
    std::vector<char> tree_path_buffer(
        tree_path_string.begin(),
        tree_path_string.end());
    tree_path_buffer.push_back('\0');

    const int tree_fd = mkstemps(tree_path_buffer.data(), 4);
    if (tree_fd < 0) {
        throw ValidationError(option_name, "failed to create temporary file for tree validation");
    }
    close(tree_fd);

    try {
        std::ofstream ofs(tree_path_buffer.data(), std::ios::trunc);
        if (!ofs) {
            throw ValidationError(option_name, "failed to open temporary file for tree validation");
        }
        ofs << tree_text;
        ofs.close();

        pll_rtree_t* rtree = pll_rtree_parse_newick(tree_path_buffer.data());
        std::remove(tree_path_buffer.data());
        if (!rtree) {
            throw ValidationError(option_name, "invalid Newick syntax");
        }
        pll_rtree_destroy(rtree, nullptr);
    } catch (...) {
        std::remove(tree_path_buffer.data());
        throw;
    }
}

std::string resolve_path(const std::filesystem::path& base, const std::string& p) {
    std::filesystem::path candidate(p);
    if (candidate.empty()) return {};
    if (candidate.is_absolute()) return candidate.string();
    return (base / candidate).string();
}

std::filesystem::path normalize_cli_path(
    const std::filesystem::path& base,
    const std::string& raw_path) {
    if (raw_path.empty()) return {};
    std::filesystem::path path = std::filesystem::path(resolve_path(base, raw_path));
    if (!path.is_absolute()) {
        path = std::filesystem::absolute(path);
    }
    return path.lexically_normal();
}

void validate_output_path(
    const std::filesystem::path& base,
    const std::string& option_name,
    const std::string& raw_path) {
    const std::filesystem::path path = normalize_cli_path(base, raw_path);
    if (path.empty() || path.filename().empty()) {
        throw ValidationError(option_name, "output path must name a file");
    }

    std::error_code ec;
    if (std::filesystem::exists(path, ec)) {
        if (ec) {
            throw ValidationError(option_name, "failed to inspect output path");
        }
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
    const std::string& option_name) {
    if (alignment.names.size() != alignment.sequences.size()) {
        throw ValidationError(option_name, "name/sequence count mismatch");
    }

    std::unordered_set<std::string> seen;
    std::vector<std::string> duplicates;
    seen.reserve(alignment.names.size() * 2);
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
        duplicates.erase(std::unique(duplicates.begin(), duplicates.end()), duplicates.end());
        throw ValidationError(
            option_name,
            "contains duplicate sequence names: " + preview_name_list(duplicates));
    }
}

void validate_alignment_symbols(
    const parse::Alignment& alignment,
    int states,
    const std::string& option_name) {
    if (states != 4) return;

    for (size_t seq_idx = 0; seq_idx < alignment.sequences.size(); ++seq_idx) {
        const std::string& name = alignment.names[seq_idx];
        const std::string& seq = alignment.sequences[seq_idx];
        for (size_t site_idx = 0; site_idx < seq.size(); ++site_idx) {
            const char c = seq[site_idx];
            if (is_valid_dna4_alignment_char(c)) continue;
            std::ostringstream oss;
            oss << "sequence '" << name << "' has unsupported DNA symbol '"
                << c << "' at site " << (site_idx + 1);
            throw ValidationError(option_name, oss.str());
        }
    }
}

void validate_model_inputs(const parse::ModelConfig& model) {
    constexpr int kMaxRateCats = 8;
    if (model.states != 4) {
        throw ValidationError("--states", "currently only 4-state DNA input is supported");
    }
    if (model.ncat <= 0 || model.ncat > kMaxRateCats) {
        throw ValidationError(
            "--ncat",
            "must be between 1 and " + std::to_string(kMaxRateCats));
    }
    if (uppercase_ascii(model.subst_model) != "GTR") {
        throw ValidationError("--subst-model", "currently only GTR is supported");
    }
    if (!std::isfinite(model.alpha) || model.alpha <= 0.0) {
        throw ValidationError("--alpha", "must be finite and > 0");
    }
    if (!std::isfinite(model.pinv) || model.pinv < 0.0 || model.pinv > 1.0) {
        throw ValidationError("--pinv", "must be finite and between 0 and 1");
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
    if (!model.rate_weights.empty()) {
        if (static_cast<int>(model.rate_weights.size()) != model.ncat) {
            throw ValidationError(
                "--rate-weights",
                "must have exactly " + std::to_string(model.ncat) + " values (ncat)");
        }
        validate_positive_vector(model.rate_weights, "--rate-weights", "rate weights");
    }
}

void validate_query_reference_name_overlap(
    const parse::Alignment& tree_alignment,
    const parse::Alignment& query_alignment,
    const std::string& option_name) {
    std::unordered_set<std::string> tree_names;
    tree_names.reserve(tree_alignment.names.size() * 2);
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
        overlaps.erase(std::unique(overlaps.begin(), overlaps.end()), overlaps.end());
        throw ValidationError(
            option_name,
            "query names overlap reference tip names, which would cause ambiguous committed tips: " +
                preview_name_list(overlaps));
    }
}

} // namespace input
} // namespace mlipper
