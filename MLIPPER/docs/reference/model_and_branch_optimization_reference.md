# Model and Branch Optimization Design Reference

This document explains why MLIPPER's global model and branch-length optimizers
use their current parameterization and initial constants. The original design
was based on RAxML-NG 0.9.0, so that version is recorded below as historical
evidence. It is not the source-of-truth for current MLIPPER runtime behavior.

Use these sources for current behavior:

- [`src/optimize/README.md`](../../src/optimize/README.md): optimizer ownership,
  entrypoints, and GPU workspaces
- [`optimization_types.hpp`](../../src/optimize/optimization_types.hpp): shared
  bounds and acceptance tolerances
- [`supported_model_scope.md`](supported_model_scope.md): models accepted by
  the public workflow

## Pinned Upstream Reference

The initial comparison used the following revisions. A newer RAxML-NG or
Coraxlib checkout may make different optimization choices.

| Component | Git revision |
| --- | --- |
| RAxML-NG tag | `0.9.0` |
| RAxML-NG | `0a064e9a40f2e00828662795141659d946440c81` |
| pll-modules | `840a19525ed133ad243afb23124911b849d4685d` |
| libpll-2 | `6874327c9ca9544e27b426c0e7205c8b6f4ed297` |
| terraphast-one | `6dcb061ceb6ed41a7ebe1a8af1f3354e0126125c` |

## Upstream Choices Retained by MLIPPER

MLIPPER retains the following optimizer structure:

1. optimize GTR substitution rates
2. optimize equilibrium frequencies
3. optimize discrete-Gamma alpha
4. optimize all branch lengths

The complete cycle repeats until the MLIPPER convergence policy stops it.

### GTR substitution rates

- Use bounded L-BFGS-B.
- Express the first five rates relative to the sixth rate.
- Fix the sixth rate at `1.0`.
- Optimize raw ratios in `[1e-3, 1000]`, not log-rates.
- Use L-BFGS-B factor `1e7` and projected-gradient tolerance `0.001`.
- Rebuild the eigendecomposition and its dependent GPU state for each
  candidate.

### Equilibrium frequencies

- Use bounded L-BFGS-B.
- Select the largest current frequency as the reference state.
- Optimize the other three frequency/reference ratios in `[1e-3, 100]`.
- Normalize all four frequencies inside the objective callback.
- Use L-BFGS-B factor `1e7` and projected-gradient tolerance `0.001`.

### Discrete-Gamma alpha

- Use bounded one-dimensional Brent minimization.
- Optimize alpha directly in `[0.0201, 100]`, not in log space.
- Use tolerance `0.001`.
- Recompute the discrete-Gamma rate categories and dependent GPU state for
  every candidate.

### Branch lengths

- Use per-branch Newton updates over all tree edges.
- Bound tree branch lengths to `[1e-6, 100]`.
- Use Newton tolerance `1e-7` and at most 30 iterations per branch.
- Use eight full-tree smoothing sweeps in each model-optimization cycle. This
  came from the RAxML-NG 0.9.0 smoothing count of 32 multiplied by `0.25`.
- Preserve the 0.9.0 raw clipped Newton displacement semantics. Do not infer
  behavior from a newer RAxML-NG release.

## Intentional MLIPPER Differences

MLIPPER follows the upstream parameterization, but it does not promise
bit-for-bit or stopping-rule equivalence with RAxML-NG 0.9.0.

| Area | RAxML-NG 0.9.0 reference | Current MLIPPER policy |
| --- | --- | --- |
| Outer likelihood epsilon | `0.1` | Configurable; default `1e-6` |
| Outer round limit | Controlled by the upstream workflow | Configurable; default 100 |
| Candidate acceptance | Upstream treeinfo behavior | Accept only a finite improving candidate using MLIPPER tolerances |
| Rejected candidate | Upstream treeinfo behavior | Restore the previously accepted model state |
| Branch acceptance | FAST Newton behavior; no newer whole-tree FAST-to-SAFE retry | Audit rebuilt likelihood and reject a decreasing result |
| GPU evaluation | Not applicable | Resident or site-batched CLVs must produce the same weighted objective |

The current tolerances belong in `optimization_types.hpp`, not in this
historical reference. If a value in this document and the implementation
differ, first determine whether the difference is intentional; do not change
the implementation solely to match this file.

## Candidate Evaluation Invariant

Regardless of resident or site-batched execution, every parameter candidate
must follow this state transition:

```text
install candidate parameter
  -> rebuild dependent model state
  -> rebuild required PMATs
  -> recompute required CLVs and scalers
  -> compute the full pattern-weighted likelihood
  -> accept the candidate or restore the accepted state
```

This invariant is more important than matching an upstream callback layout.
A likelihood computed from stale PMATs or CLVs is not a valid optimizer
objective.

## Compatibility Boundary

The original comparison covered GTR rates, DNA equilibrium frequencies,
discrete-Gamma alpha, and iterative branch lengths. It did not validate PINV,
joint alpha/PINV optimization, FreeRate rates or weights, mixture models, or
multi-partition proportional branch-length scalers. See
[`supported_model_scope.md`](supported_model_scope.md) before expanding the
public model surface.

RAxML-NG optimizes only parameters marked as free/ML. Two programs reading
similar model text can still use different fixed/free masks, so a future
compatibility experiment must compare that mask explicitly as well as the
starting values.
