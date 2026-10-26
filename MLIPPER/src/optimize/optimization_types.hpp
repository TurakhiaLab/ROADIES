#pragma once

namespace mlipper {

namespace optimization::model_parameters {

// Bounds keep the GTR parameterization away from singular zero/infinite ratios.
// Acceptance tolerances are log-likelihood deltas, not relative errors.
inline constexpr double kMinimumRateRatio = 1.0e-3;
inline constexpr double kMaximumRateRatio = 1000.0;
inline constexpr double kMinimumFrequencyRatio = 1.0e-3;
inline constexpr double kMaximumFrequencyRatio = 100.0;
inline constexpr double kMinimumAlpha = 0.0201;
inline constexpr double kMaximumAlpha = 100.0;
inline constexpr double kAlphaTolerance = 0.001;
inline constexpr double kModelAcceptanceTolerance = 1.0e-10;
inline constexpr double kAlphaAcceptanceTolerance = 1.0e-8;
inline constexpr double kLikelihoodDecreaseTolerance = 1.0e-8;

} // namespace optimization::model_parameters

namespace optimization::branch_lengths {

// Branch lengths are expected substitutions per site. Placement uses a larger
// floor for newly inserted edges; topology refinement may collapse internal
// edges down to kTreeMinimumLength.
inline constexpr double kPlacementMinimumLength = 1.0e-4;
inline constexpr double kTreeMinimumLength = 1.0e-6;
inline constexpr double kMaximumLength = 100.0;
inline constexpr double kDefaultLength = 0.10536051565782628; // -log(0.9)
inline constexpr double kPlacementNewtonTolerance = 1.0e-5;
inline constexpr double kNewtonTolerance = 1.0e-7;
inline constexpr double kCandidateAcceptanceTolerance = 1.0e-8;
inline constexpr double kFullTreeAuditDecreaseTolerance = 1.0e-8;
inline constexpr double kLengthChangeTolerance = 1.0e-12;

} // namespace optimization::branch_lengths

// Jacobi scores a group from one frozen likelihood state; Sequential commits
// each accepted coordinate before scoring the next one.
enum class EdgeUpdateScheme {
    Jacobi,
    Sequential,
};

// Selects whether an optimization rebuild may retain immutable tip CLVs.
enum class ClvRetention {
    RebuildAll,
    PreserveTips,
};

// Controls whether a proposal is audited on its local subtree or against the
// complete resident tree likelihood.
enum class AcceptanceScope {
    LocalSubtree,
    FullTree,
};

// Shared control surface for resident branch-length optimizers.
struct BranchOptimizationOptions {
    int sweeps = 8;
    int newton_iterations = 30;
    double likelihood_tolerance = 1.0e-6;
    EdgeUpdateScheme update_scheme = EdgeUpdateScheme::Sequential;
    ClvRetention clv_retention = ClvRetention::RebuildAll;
    AcceptanceScope acceptance_scope = AcceptanceScope::FullTree;
};

struct BranchLengthOptimizationResult {
    int sweeps = 0;
    int newton_iterations = 0;
    double log_likelihood_before = 0.0;
    double log_likelihood_after = 0.0;
    double derivative_seconds = 0.0;
    bool attempted = false;
    bool accepted = false;
};

} // namespace mlipper
