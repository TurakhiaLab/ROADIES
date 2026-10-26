#pragma once

#include <functional>
#include <vector>
#include "optimize/optimization_types.hpp"

namespace mlipper {

// Result objects retain both sides of an attempted coordinate update so logs
// and callers can distinguish rejected proposals from accepted no-op changes.
struct EquilibriumFrequencyOptimizationResult {
    std::vector<double> initial_frequencies;
    std::vector<double> frequencies;
    double true_log_likelihood_before = 0.0;
    double true_log_likelihood_after = 0.0;
    int evaluations = 0;
    bool accepted = false;
};

struct GtrRateOptimizationResult {
    std::vector<double> initial_rates;
    std::vector<double> rates;
    double log_likelihood_before = 0.0;
    double log_likelihood_after = 0.0;
    int evaluations = 0;
    bool accepted = false;
};

struct GammaAlphaOptimizationResult {
    double initial_alpha = 1.0;
    double alpha = 1.0;
    double log_likelihood_before = 0.0;
    double log_likelihood_after = 0.0;
    int evaluations = 0;
    bool accepted = false;
};

// Automatic selects resident or site-batched execution from available VRAM.
// Explicit modes are primarily reproducibility and diagnostics controls.
enum class GlobalOptimizationBackend {
    Automatic,
    ResidentSequential,
    SiteBatchedSequential,
};

// Alternates model coordinates and branch-length sweeps until the full-tree
// likelihood improvement is at most likelihood_tolerance or max_rounds ends.
struct FinalModelOptimizationOptions {
    int max_rounds = 100;
    int branch_length_sweeps = 8;
    int branch_length_newton_iterations = 30;
    double likelihood_tolerance = 1.0e-6;
    bool optimize_model_parameters = true;
    bool optimize_branch_lengths = true;
    GlobalOptimizationBackend backend = GlobalOptimizationBackend::Automatic;
};

struct FinalModelOptimizationResult {
    int rounds = 0;
    double initial_log_likelihood = 0.0;
    double final_log_likelihood = 0.0;
    bool converged = false;
    std::vector<EquilibriumFrequencyOptimizationResult> frequency_updates;
    std::vector<GtrRateOptimizationResult> gtr_rate_updates;
    std::vector<GammaAlphaOptimizationResult> alpha_updates;
    std::vector<mlipper::BranchLengthOptimizationResult> branch_updates;
};

} // namespace mlipper

namespace mlipper::model_optimization {

struct GlobalModelOptimizerCallbacks {
    // Candidate installation is stateful. The optimizer always reinstalls the
    // accepted candidate or the original values before a normal return.
    std::function<double()> evaluate;
    std::function<void(
        const std::vector<double>&, const std::vector<double>&)>
        install_substitution_model;
    std::function<void(double)> install_alpha;
    std::function<mlipper::BranchLengthOptimizationResult(int, int)>
        optimize_branch_lengths;
    std::function<void()> branch_lengths_accepted;
};

// CPU search coordinator. Numerical likelihood state lives behind callbacks;
// accepted values are written into the referenced frequency/rate/alpha inputs.
class GlobalModelOptimizer {
public:
    GtrRateOptimizationResult optimizeSubstitutionRates(
        const std::vector<double>& frequencies,
        const std::vector<double>& rates,
        const GlobalModelOptimizerCallbacks&) const;
    EquilibriumFrequencyOptimizationResult optimizeFrequencies(
        const std::vector<double>& frequencies,
        const std::vector<double>& rates,
        const GlobalModelOptimizerCallbacks&) const;
    GammaAlphaOptimizationResult optimizeAlpha(
        double alpha, int rate_categories,
        const GlobalModelOptimizerCallbacks&) const;
    FinalModelOptimizationResult optimize(
        std::vector<double>& frequencies,
        std::vector<double>& rates,
        double& alpha,
        int states,
        int rate_categories,
        const FinalModelOptimizationOptions&,
        const GlobalModelOptimizerCallbacks&) const;
};

} // namespace mlipper::model_optimization
