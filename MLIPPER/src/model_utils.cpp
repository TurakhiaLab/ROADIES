#include "model_utils.hpp"

#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "tree/tree.hpp"

namespace {

constexpr double kEmpiricalFrequencyFloor = 1e-8;

void normalize_vector(std::vector<double>& vec) {
    double sum = 0.0;
    for (double v : vec) sum += v;
    if (sum <= 0.0) return;
    for (double& v : vec) v /= sum;
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

bool add_empirical_dna4_counts(char c, std::vector<double>& counts) {
    const uint8_t mask = encode_state_DNA4_mask(c);
    int bit_count = 0;
    for (int state = 0; state < 4; ++state) {
        if (mask & (1u << state)) {
            ++bit_count;
        }
    }
    if (bit_count == 0) {
        return false;
    }

    const double share = 1.0 / static_cast<double>(bit_count);
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
    const std::regex model_re(R"(^([A-Za-z0-9_]+)\{([^}]*)\})");
    if (!std::regex_search(model_text, match, model_re)) {
        throw std::runtime_error("Could not parse substitution model from bestModel: " + model_text);
    }

    BestModelConfig parsed;
    parsed.raw_text = model_text;
    parsed.model.states = 4;
    parsed.model.subst_model = match[1].str();
    parsed.model.rates = parse_slash_floats(match[2].str(), "GTR rates");

    const std::regex gamma_re(R"(\+G(\d+)[^{}]*\{([^}]*)\})");
    if (std::regex_search(model_text, match, gamma_re)) {
        parsed.model.ncat = std::stoi(match[1].str());
        parsed.model.alpha = parse_double_strict(match[2].str(), "gamma alpha");
    } else {
        parsed.model.ncat = 1;
        parsed.model.alpha = 1.0;
    }

    // bestModel -> pinv parsing is intentionally disabled for now.
    // Keep the runtime on the existing explicit/default CLI path until
    // nonzero pinv is wired and validated end-to-end.
    parsed.model.pinv = 0.0;

    const std::regex empirical_freq_re(R"(\+FC(?:\+|$))");
    if (std::regex_search(model_text, empirical_freq_re)) {
        parsed.empirical_freqs = true;
        parsed.model.freqs.clear();
    } else {
        const std::regex manual_freq_re(R"(\+(?:FU|FO|F)\{([^}]*)\})");
        if (std::regex_search(model_text, match, manual_freq_re)) {
            parsed.model.freqs = parse_slash_floats(match[1].str(), "base frequencies");
        } else {
            parsed.model.freqs.clear();
        }
    }

    if (parsed.model.subst_model != "GTR") {
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
        throw std::runtime_error("Manual DNA frequencies in bestModel must contain exactly 4 values");
    }

    return parsed;
}

std::vector<double> estimate_empirical_pi(const parse::Alignment& alignment, int states) {
    if (states != 4 && states != 5) {
        throw std::runtime_error(
            "Empirical frequency estimation currently supports only 4-state or 5-state DNA data.");
    }

    std::vector<double> counts(states, 0.0);
    double informative_weight = 0.0;
    for (const std::string& seq : alignment.sequences) {
        for (char c : seq) {
            if (states == 4) {
                if (add_empirical_dna4_counts(c, counts)) {
                    informative_weight += 1.0;
                }
                continue;
            }

            switch (c) {
                case 'A':
                case 'a':
                    counts[0] += 1.0;
                    break;
                case 'C':
                case 'c':
                    counts[1] += 1.0;
                    break;
                case 'G':
                case 'g':
                    counts[2] += 1.0;
                    break;
                case 'T':
                case 't':
                case 'U':
                case 'u':
                    counts[3] += 1.0;
                    break;
                case '-':
                case '.':
                    counts[4] += 1.0;
                    break;
                default:
                    continue;
            }
            informative_weight += 1.0;
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
    if (static_cast<int>(pi.size()) != states) pi.assign((size_t)states, 1.0 / states);
    normalize_vector(pi);
    return pi;
}

} // namespace model
} // namespace mlipper
