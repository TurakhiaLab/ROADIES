#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "io/parse_file.hpp"

namespace mlipper {
namespace model {

struct BestModelConfig {
    parse::ModelConfig model;
    bool empirical_freqs = false;
};

// Parse the supported RAxML-NG bestModel subset. empirical_freqs records +FC;
// unsupported rate heterogeneity, mixtures, or invariant sites are rejected.
BestModelConfig parse_best_model_file(const std::filesystem::path& path);

// Gamma categories are normalized to mean one; their mixture weights are
// uniform because only mean-discretized +G models are accepted.
std::vector<double> build_discrete_gamma_weights(int rate_categories);
std::vector<double> build_gamma_rate_categories(
    double alpha, int rate_categories);

// Pattern weights restore the multiplicity lost during column compression.
std::vector<double> estimate_empirical_pi(
    const parse::Alignment& alignment,
    int states,
    const std::vector<unsigned>& pattern_weights = {});
// Validate positivity/finite values, apply the supported state dimension, and
// return frequencies normalized to sum exactly one within floating precision.
std::vector<double> ensure_normalized_pi(std::vector<double> pi, int states);

} // namespace model
} // namespace mlipper
