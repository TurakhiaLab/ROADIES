#include "util/model_utils.hpp"

#include <cmath>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <corax/corax_model.h>

#include "tree/tree.hpp"

namespace {

constexpr double kEmpiricalFrequencyFloor = 1e-8;

void normalize_vector(std::vector<double>& vec) {
    double sum = 0.0;
    for (double value : vec) {
        sum += value;
    }
    if (!std::isfinite(sum) || sum <= 0.0) {
        throw std::runtime_error("Cannot normalize non-positive or non-finite values.");
    }
    for (double& value : vec) {
        value /= sum;
    }
}

void floor_zero_entries(std::vector<double>& vec, double floor_value) {
    bool adjusted = false;
    for (double& v : vec) {
        if (v <= 0.0) {
            v = floor_value;
            adjusted = true;
        }
    }
    if (adjusted) normalize_vector(vec);
}

double parse_double_strict(const std::string& text, const char* label) {
    size_t consumed = 0;
    double value = 0.0;
    try {
        value = std::stod(text, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Failed to parse ") + label + ": " + text);
    }
    if (consumed != text.size()) {
        throw std::runtime_error(std::string("Failed to parse ") + label + ": " + text);
    }
    return value;
}

std::vector<double> parse_slash_floats(const std::string& text, const char* label) {
    std::vector<double> values;
    std::stringstream ss(text);
    std::string part;
    while (std::getline(ss, part, '/')) {
        if (part.empty()) continue;
        values.push_back(parse_double_strict(part, label));
    }
    if (values.empty()) {
        throw std::runtime_error(std::string(label) + " is empty");
    }
    return values;
}

bool add_empirical_dna4_counts(
    char c, double observation_weight, std::vector<double>& counts)
{
    const uint8_t mask = encode_state_DNA4_mask(c);
    int bit_count = 0;
    for (int state = 0; state < 4; ++state) {
        if (mask & (1u << state)) {
            ++bit_count;
        }
    }
    // Fully ambiguous symbols (N, gap, dot, or unknown input) carry no
    // empirical base-frequency information.  Partial IUPAC ambiguities still
    // distribute one observation across the states they represent.
    if (bit_count == 0 || bit_count == 4) {
        return false;
    }

    const double share =
        observation_weight / static_cast<double>(bit_count);
    for (int state = 0; state < 4; ++state) {
        if (mask & (1u << state)) {
            counts[static_cast<size_t>(state)] += share;
        }
    }
    return true;
}

} // namespace

namespace mlipper {
namespace model {

BestModelConfig parse_best_model_file(const std::filesystem::path& path) {
    std::ifstream handle(path);
    if (!handle) {
        throw std::runtime_error("Cannot open bestModel file: " + path.string());
    }

    std::string line;
    while (std::getline(handle, line)) {
        if (!line.empty()) break;
    }
    if (line.empty()) {
        throw std::runtime_error("bestModel file is empty: " + path.string());
    }

    const std::string model_text = line.substr(0, line.find(','));
    std::smatch match;

    const std::regex free_rate_re(R"(\+R\d+(?:\{|\+|$))", std::regex::icase);
    if (std::regex_search(model_text, free_rate_re)) {
        throw std::runtime_error(
            "FreeRate (+R) models are not supported; use discrete-Gamma (+G).");
    }
    const std::regex mixture_re(
        R"((?:^|\+)MIX(?:TURE)?(?:\{|\(|\+|$))",
        std::regex::icase);
    if (std::regex_search(model_text, mixture_re)) {
        throw std::runtime_error(
            "Mixture models are not supported; use a single GTR discrete-Gamma model.");
    }
    const std::regex invariant_sites_re(
        R"(\+I(?:\{[^}]*\})?(?:\+|$))", std::regex::icase);
    if (std::regex_search(model_text, invariant_sites_re)) {
        throw std::runtime_error(
            "Invariant-site (+I) models are not supported; use a bestModel without +I.");
    }

    // Keep the accepted grammar intentionally narrower than corax/RAxML-NG:
    // one GTR model, optional DNA frequencies, and optional mean-discretized
    // Gamma rates.  Reject modifiers that MLIPPER cannot evaluate instead of
    // silently ignoring them.
    const std::regex supported_model_re(
        R"(^[A-Za-z0-9_]+\{[^}]*\}(?:\+(?:FC|(?:FU|FO|F)\{[^}]*\}))?(?:\+G\d+m\{[^}]*\})?$)",
        std::regex::icase);
    if (!std::regex_match(model_text, supported_model_re)) {
        throw std::runtime_error(
            "Unsupported bestModel syntax or modifier: " + model_text +
            ". MLIPPER supports GTR with optional +FC/+F{...} and +G<n>m{alpha}.");
    }

    const std::regex model_re(R"(^([A-Za-z0-9_]+)\{([^}]*)\})");
    if (!std::regex_search(model_text, match, model_re)) {
        throw std::runtime_error("Could not parse substitution model from bestModel: " + model_text);
    }

    BestModelConfig parsed;
    parsed.model.states = 4;
    parsed.model.subst_model = match[1].str();
    parsed.model.rates = parse_slash_floats(match[2].str(), "GTR rates");

    const std::regex gamma_re(R"(\+G(\d+)m\{([^}]*)\})", std::regex::icase);
    if (std::regex_search(model_text, match, gamma_re)) {
        parsed.model.ncat = std::stoi(match[1].str());
        parsed.model.alpha = parse_double_strict(match[2].str(), "gamma alpha");
    } else {
        parsed.model.ncat = 1;
        parsed.model.alpha = 1.0;
    }

    parsed.model.pinv = 0.0;

    const std::regex empirical_freq_re(
        R"(\+FC(?:\+|$))", std::regex::icase);
    if (std::regex_search(model_text, empirical_freq_re)) {
        parsed.empirical_freqs = true;
        parsed.model.freqs.clear();
    } else {
        const std::regex manual_freq_re(
            R"(\+(?:FU|FO|F)\{([^}]*)\})", std::regex::icase);
        if (std::regex_search(model_text, match, manual_freq_re)) {
            parsed.model.freqs = parse_slash_floats(match[1].str(), "base frequencies");
        } else {
            parsed.model.freqs.clear();
        }
    }

    const std::regex gtr_re(R"(^GTR$)", std::regex::icase);
    if (!std::regex_match(parsed.model.subst_model, gtr_re)) {
        throw std::runtime_error(
            "Unsupported substitution model in bestModel: " + parsed.model.subst_model +
            ". Currently only GTR is supported.");
    }
    if (parsed.model.rates.size() != 6) {
        throw std::runtime_error(
            "GTR bestModel must provide exactly 6 rates, got " +
            std::to_string(parsed.model.rates.size()));
    }
    if (!parsed.model.freqs.empty() && parsed.model.freqs.size() != 4) {
        throw std::runtime_error(
            "Manual DNA frequencies in bestModel must contain exactly 4 values");
    }
    for (double rate : parsed.model.rates) {
        if (!std::isfinite(rate) || rate <= 0.0) {
            throw std::runtime_error(
                "GTR rates in bestModel must be finite and positive.");
        }
    }
    if (!std::isfinite(parsed.model.alpha) || parsed.model.alpha <= 0.0) {
        throw std::runtime_error(
            "Gamma alpha in bestModel must be finite and positive.");
    }
    for (double frequency : parsed.model.freqs) {
        if (!std::isfinite(frequency) || frequency < 0.0) {
            throw std::runtime_error(
                "Manual DNA frequencies in bestModel must be finite and non-negative.");
        }
    }
    if (!parsed.model.freqs.empty()) {
        // Normalize at the parser boundary so every downstream consumer sees
        // the same stationary distribution, independent of input scale.
        normalize_vector(parsed.model.freqs);
    }

    return parsed;
}

std::vector<double> build_discrete_gamma_weights(int rate_categories) {
    if (rate_categories <= 0) return {};
    return std::vector<double>(
        static_cast<size_t>(rate_categories),
        1.0 / static_cast<double>(rate_categories));
}

std::vector<double> build_gamma_rate_categories(
    double alpha, int rate_categories)
{
    if (rate_categories <= 0) return {};

    std::vector<double> rates(static_cast<size_t>(rate_categories), 1.0);
    if (!corax_compute_gamma_cats(
            alpha,
            static_cast<unsigned int>(rate_categories),
            rates.data(),
            CORAX_GAMMA_RATES_MEAN)) {
        throw std::runtime_error(
            "Discrete-Gamma rate-category computation failed for alpha=" +
            std::to_string(alpha));
    }

    double sum = 0.0;
    for (double rate : rates) {
        if (!std::isfinite(rate) || rate <= 0.0) {
            throw std::runtime_error(
                "Discrete-Gamma rate-category computation returned an invalid rate");
        }
        sum += rate;
    }
    if (!std::isfinite(sum) || sum <= 0.0) {
        throw std::runtime_error(
            "Discrete-Gamma rate-category computation returned an invalid mean");
    }
    const double mean_scale = static_cast<double>(rate_categories) / sum;
    for (double& rate : rates) rate *= mean_scale;
    return rates;
}

std::vector<double> estimate_empirical_pi(
    const parse::Alignment& alignment,
    int states,
    const std::vector<unsigned>& pattern_weights)
{
    if (states != 4 && states != 5) {
        throw std::runtime_error(
            "Empirical frequency estimation currently supports only 4-state or 5-state DNA data.");
    }

    std::vector<double> counts(states, 0.0);
    double informative_weight = 0.0;
    if (!pattern_weights.empty() &&
        pattern_weights.size() != alignment.sites) {
        throw std::runtime_error(
            "Empirical frequency pattern weights do not match alignment sites.");
    }
    for (const std::string& seq : alignment.sequences) {
        for (size_t site = 0; site < seq.size(); ++site) {
            const char c = seq[site];
            const double weight = pattern_weights.empty()
                ? 1.0
                : static_cast<double>(pattern_weights[site]);
            if (states == 4) {
                if (add_empirical_dna4_counts(c, weight, counts)) {
                    informative_weight += weight;
                }
                continue;
            }

            switch (c) {
                case 'A':
                case 'a':
                    counts[0] += weight;
                    break;
                case 'C':
                case 'c':
                    counts[1] += weight;
                    break;
                case 'G':
                case 'g':
                    counts[2] += weight;
                    break;
                case 'T':
                case 't':
                case 'U':
                case 'u':
                    counts[3] += weight;
                    break;
                case '-':
                case '.':
                    counts[4] += weight;
                    break;
                default:
                    continue;
            }
            informative_weight += weight;
        }
    }
    if (informative_weight <= 0.0) {
        throw std::runtime_error(
            "Cannot estimate empirical frequencies: alignment contains no informative states.");
    }
    normalize_vector(counts);
    floor_zero_entries(counts, kEmpiricalFrequencyFloor);
    return counts;
}

std::vector<double> ensure_normalized_pi(std::vector<double> pi, int states) {
    if (states <= 0) {
        throw std::runtime_error("State count must be positive.");
    }
    if (pi.size() != static_cast<size_t>(states)) {
        pi.assign(
            static_cast<size_t>(states),
            1.0 / static_cast<double>(states));
    }
    for (double frequency : pi) {
        if (!std::isfinite(frequency) || frequency < 0.0) {
            throw std::runtime_error(
                "State frequencies must be finite and non-negative.");
        }
    }
    normalize_vector(pi);
    return pi;
}

} // namespace model
} // namespace mlipper
