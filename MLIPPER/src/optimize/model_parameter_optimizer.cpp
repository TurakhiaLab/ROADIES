#include "model_parameter_optimizer.hpp"
#include "optimization_types.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>

#include <corax/corax_optimize.h>

namespace mlipper::model_optimization {

using namespace mlipper::optimization::model_parameters;

namespace {
struct BoundedOptimizationResult {
    std::vector<double> parameters;
    double objective = 0.0;
    int evaluations = 0;
};

struct OneDimensionalOptimizationResult {
    double parameter = 0.0;
    double objective = 0.0;
    double second_derivative = 0.0;
    int evaluations = 0;
};

struct BoundedContext {
    const std::function<double(const std::vector<double>&)>* objective;
    std::size_t dimensions;
    int evaluations = 0;
    std::exception_ptr exception;
};

double evaluateBounded(void* opaque, double* parameters)
{
    auto& context = *static_cast<BoundedContext*>(opaque);
    ++context.evaluations;
    if (context.exception) {
        return std::numeric_limits<double>::infinity();
    }
    try {
        return (*context.objective)(std::vector<double>(
            parameters, parameters + context.dimensions));
    } catch (...) {
        // Never unwind a C++ exception through the C optimizer. Preserve the
        // first failure and make this candidate maximally unattractive until
        // control returns to the C++ boundary below.
        if (!context.exception) {
            context.exception = std::current_exception();
        }
        return std::numeric_limits<double>::infinity();
    }
}

struct OneDimensionalContext {
    const std::function<double(double)>* objective;
    int evaluations = 0;
    std::exception_ptr exception;
};

double evaluateOneDimensional(void* opaque, double parameter)
{
    auto& context = *static_cast<OneDimensionalContext*>(opaque);
    ++context.evaluations;
    if (context.exception) {
        return std::numeric_limits<double>::infinity();
    }
    try {
        return (*context.objective)(parameter);
    } catch (...) {
        if (!context.exception) {
            context.exception = std::current_exception();
        }
        return std::numeric_limits<double>::infinity();
    }
}

BoundedOptimizationResult minimizeBoundedWithCorax(
    const std::vector<double>& initial,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::function<double(const std::vector<double>&)>& objective,
    double factr = 1.0e7,
    double projected_gradient_tolerance = 0.001)
{
    if (initial.empty() || initial.size() != lower.size() ||
        initial.size() != upper.size()) {
        throw std::invalid_argument(
            "minimizeBoundedWithCorax: inconsistent dimensions");
    }
    if (!objective ||
        initial.size() >
            static_cast<size_t>(std::numeric_limits<unsigned int>::max())) {
        throw std::invalid_argument(
            "minimizeBoundedWithCorax: invalid objective or dimension");
    }
    BoundedOptimizationResult result;
    result.parameters = initial;
    std::vector<double> mutable_lower = lower;
    std::vector<double> mutable_upper = upper;
    std::vector<int> bounds(initial.size(), 2);
    BoundedContext context{&objective, initial.size(), 0, nullptr};
    result.objective = corax_opt_minimize_lbfgsb(
        result.parameters.data(), mutable_lower.data(), mutable_upper.data(),
        bounds.data(), static_cast<unsigned int>(initial.size()),
        factr, projected_gradient_tolerance, &context, evaluateBounded);
    if (context.exception) {
        std::rethrow_exception(context.exception);
    }
    result.evaluations = context.evaluations;
    return result;
}

OneDimensionalOptimizationResult minimizeOneDimensionalWithCorax(
    double lower,
    double initial,
    double upper,
    double tolerance,
    const std::function<double(double)>& objective)
{
    if (!(lower <= initial && initial <= upper) || !(tolerance > 0.0)) {
        throw std::invalid_argument(
            "minimizeOneDimensionalWithCorax: invalid interval");
    }
    OneDimensionalOptimizationResult result;
    result.parameter = initial;
    if (!objective) {
        throw std::invalid_argument(
            "minimizeOneDimensionalWithCorax: empty objective");
    }
    OneDimensionalContext context{&objective, 0, nullptr};
    result.parameter = corax_opt_minimize_brent(
        lower, initial, upper, tolerance,
        &result.objective, &result.second_derivative,
        &context, evaluateOneDimensional);
    if (context.exception) {
        std::rethrow_exception(context.exception);
    }
    result.evaluations = context.evaluations;
    return result;
}

} // namespace


GtrRateOptimizationResult GlobalModelOptimizer::optimizeSubstitutionRates(
    const std::vector<double>& frequencies,
    const std::vector<double>& rates,
    const GlobalModelOptimizerCallbacks& callbacks) const
{
    if (rates.size() != 6) {
        throw std::invalid_argument(
            "RAxML-NG 0.9.0 GTR optimization requires six rates");
    }
    GtrRateOptimizationResult result;
    result.initial_rates = rates;
    result.log_likelihood_before = callbacks.evaluate();
    const double reference_rate = rates.back();
    if (!(reference_rate > 0.0)) {
        throw std::invalid_argument("invalid GTR reference rate");
    }
    std::vector<double> initial(5);
    for (int i = 0; i < 5; ++i) {
        initial[i] = std::clamp(
            rates[i] / reference_rate,
            kMinimumRateRatio, kMaximumRateRatio);
    }
    auto expand_rates = [](const std::vector<double>& parameters) {
        std::vector<double> candidate(6, 1.0);
        std::copy(parameters.begin(), parameters.end(), candidate.begin());
        return candidate;
    };
    BoundedOptimizationResult optimization;
    try {
        optimization = minimizeBoundedWithCorax(
            initial, std::vector<double>(5, kMinimumRateRatio),
            std::vector<double>(5, kMaximumRateRatio),
            [&](const std::vector<double>& parameters) {
                callbacks.install_substitution_model(
                    frequencies, expand_rates(parameters));
                return -callbacks.evaluate();
            });
    } catch (...) {
        const std::exception_ptr failure = std::current_exception();
        try {
            callbacks.install_substitution_model(frequencies, rates);
        } catch (...) {
            // Preserve the objective failure; device-level partial update is
            // documented separately because rollback can fail for the same reason.
        }
        std::rethrow_exception(failure);
    }
    result.evaluations = optimization.evaluations;
    result.rates = expand_rates(optimization.parameters);
    result.log_likelihood_after = -optimization.objective;
    if (std::isfinite(result.log_likelihood_after) &&
        result.log_likelihood_after >=
            result.log_likelihood_before + kModelAcceptanceTolerance) {
        result.accepted = true;
        callbacks.install_substitution_model(frequencies, result.rates);
    } else {
        result.rates = result.initial_rates;
        result.log_likelihood_after = result.log_likelihood_before;
        callbacks.install_substitution_model(frequencies, result.initial_rates);
    }
    return result;
}

EquilibriumFrequencyOptimizationResult GlobalModelOptimizer::optimizeFrequencies(
    const std::vector<double>& frequencies,
    const std::vector<double>& rates,
    const GlobalModelOptimizerCallbacks& callbacks) const
{
    if (frequencies.size() != 4) {
        throw std::invalid_argument(
            "RAxML-NG-compatible frequency optimization currently requires DNA");
    }
    EquilibriumFrequencyOptimizationResult result;
    result.initial_frequencies = frequencies;
    result.true_log_likelihood_before = callbacks.evaluate();
    const int reference_state = static_cast<int>(
        std::distance(frequencies.begin(),
            std::max_element(frequencies.begin(), frequencies.end())));
    std::vector<double> initial;
    initial.reserve(3);
    for (int state = 0; state < 4; ++state) {
        if (state != reference_state) {
            initial.push_back(std::clamp(
                frequencies[state] / frequencies[reference_state],
                kMinimumFrequencyRatio, kMaximumFrequencyRatio));
        }
    }
    auto expand_frequencies =
        [reference_state](const std::vector<double>& ratios) {
            std::vector<double> candidate(4, 0.0);
            candidate[reference_state] = 1.0;
            double sum = 1.0;
            for (int state = 0, parameter = 0; state < 4; ++state) {
                if (state == reference_state) continue;
                candidate[state] = ratios[parameter++];
                sum += candidate[state];
            }
            for (double& frequency : candidate) frequency /= sum;
            return candidate;
        };
    BoundedOptimizationResult optimization;
    try {
        optimization = minimizeBoundedWithCorax(
            initial, std::vector<double>(3, kMinimumFrequencyRatio),
            std::vector<double>(3, kMaximumFrequencyRatio),
            [&](const std::vector<double>& ratios) {
                callbacks.install_substitution_model(
                    expand_frequencies(ratios), rates);
                return -callbacks.evaluate();
            });
    } catch (...) {
        const std::exception_ptr failure = std::current_exception();
        try {
            callbacks.install_substitution_model(frequencies, rates);
        } catch (...) {
            // Preserve the objective failure if restoring the same backend also fails.
        }
        std::rethrow_exception(failure);
    }
    result.evaluations = optimization.evaluations;
    const std::vector<double> candidate =
        expand_frequencies(optimization.parameters);
    const double candidate_logl = -optimization.objective;
    if (std::isfinite(candidate_logl) &&
        candidate_logl >= result.true_log_likelihood_before +
            kModelAcceptanceTolerance) {
        result.frequencies = candidate;
        result.true_log_likelihood_after = candidate_logl;
        result.accepted = true;
        callbacks.install_substitution_model(candidate, rates);
    } else {
        result.frequencies = result.initial_frequencies;
        result.true_log_likelihood_after = result.true_log_likelihood_before;
        callbacks.install_substitution_model(result.initial_frequencies, rates);
    }
    return result;
}

GammaAlphaOptimizationResult GlobalModelOptimizer::optimizeAlpha(
    double alpha, int rate_categories,
    const GlobalModelOptimizerCallbacks& callbacks) const
{
    GammaAlphaOptimizationResult result;
    result.initial_alpha = alpha;
    result.alpha = alpha;
    result.log_likelihood_before = callbacks.evaluate();
    result.log_likelihood_after = result.log_likelihood_before;
    if (rate_categories <= 1) {
        return result;
    }

    OneDimensionalOptimizationResult optimization;
    try {
        optimization = minimizeOneDimensionalWithCorax(
            kMinimumAlpha, std::clamp(alpha, kMinimumAlpha, kMaximumAlpha),
            kMaximumAlpha, kAlphaTolerance,
            [&](double candidate) {
                callbacks.install_alpha(candidate);
                return -callbacks.evaluate();
            });
    } catch (...) {
        const std::exception_ptr failure = std::current_exception();
        try {
            callbacks.install_alpha(alpha);
        } catch (...) {
            // Preserve the objective failure if restoring the same backend also fails.
        }
        std::rethrow_exception(failure);
    }
    result.evaluations = optimization.evaluations;
    const double candidate_logl = -optimization.objective;
    if (std::isfinite(candidate_logl) &&
        candidate_logl >= result.log_likelihood_before +
            kAlphaAcceptanceTolerance) {
        result.alpha = optimization.parameter;
        result.log_likelihood_after = candidate_logl;
        result.accepted = true;
        callbacks.install_alpha(result.alpha);
    } else {
        callbacks.install_alpha(result.initial_alpha);
    }
    return result;
}

// Coordinate-descent driver. Each sub-optimizer is responsible for restoring
// rejected values through its callbacks. A round is committed only after the
// complete model/branch state has a finite nondecreasing full-tree likelihood.
FinalModelOptimizationResult GlobalModelOptimizer::optimize(
    std::vector<double>& frequencies,
    std::vector<double>& rates,
    double& alpha,
    int states,
    int rate_categories,
    const FinalModelOptimizationOptions& options,
    const GlobalModelOptimizerCallbacks& callbacks) const
{
    if (states != 4 || options.max_rounds < 1 ||
        options.branch_length_sweeps < 1 ||
        options.branch_length_newton_iterations < 1 ||
        !(options.likelihood_tolerance >= 0.0)) {
        throw std::invalid_argument("GlobalModelOptimizer: invalid inputs");
    }
    FinalModelOptimizationResult result;
    result.initial_log_likelihood = callbacks.evaluate();
    result.final_log_likelihood = result.initial_log_likelihood;
    for (int round = 0; round < options.max_rounds; ++round) {
        const double round_start = result.final_log_likelihood;
        if (options.optimize_model_parameters) {
            result.gtr_rate_updates.push_back(
                optimizeSubstitutionRates(frequencies, rates, callbacks));
            rates = result.gtr_rate_updates.back().rates;
            result.frequency_updates.push_back(
                optimizeFrequencies(frequencies, rates, callbacks));
            frequencies = result.frequency_updates.back().frequencies;
            result.alpha_updates.push_back(
                optimizeAlpha(alpha, rate_categories, callbacks));
            alpha = result.alpha_updates.back().alpha;
        }
        if (options.optimize_branch_lengths) {
            result.branch_updates.push_back(
                callbacks.optimize_branch_lengths(
                    options.branch_length_sweeps,
                    options.branch_length_newton_iterations));
            if (result.branch_updates.back().accepted &&
                callbacks.branch_lengths_accepted) {
                callbacks.branch_lengths_accepted();
            }
        }
        if (options.optimize_branch_lengths &&
            !result.branch_updates.empty()) {
            result.final_log_likelihood =
                result.branch_updates.back().log_likelihood_after;
        } else {
            result.final_log_likelihood = callbacks.evaluate();
        }
        result.rounds = round + 1;
        const double improvement =
            result.final_log_likelihood - round_start;
        if (improvement < -kLikelihoodDecreaseTolerance) {
            throw std::runtime_error(
                "global model optimization decreased likelihood");
        }
        // The branch optimizer's convergence flag is local to one inner
        // invocation.  Rebuilt CLVs can expose additional coordinate updates
        // on the next outer round, so only the audited full-tree likelihood
        // improvement may terminate this loop early.
        if (improvement <= options.likelihood_tolerance) {
            result.converged = true;
            break;
        }
    }
    return result;
}

} // namespace mlipper::model_optimization
