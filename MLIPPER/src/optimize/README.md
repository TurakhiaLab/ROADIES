# Model and Branch-Length Optimization

This directory separates numerical search from the GPU/tree state needed to evaluate a candidate model.

## Responsibility Boundary

```text
GlobalModelOptimizer
  numerical parameterization, search, acceptance, outer convergence
        |
        v callbacks
ModelOptimizationBackend
  install candidate model, rebuild PMAT/CLV state, evaluate likelihood,
  invoke the selected branch optimizer
```

`GlobalModelOptimizer` must not know about `MlipperSession`, CUDA buffers, or resident versus site-batched CLVs. `ModelOptimizationBackend` adapts a supplied `ModelOptimizationContext`; it does not own the session.

## Files

- `model_parameter_optimizer.hpp/.cpp`: GTR-rate, equilibrium-frequency, and Gamma-alpha numerical optimization plus the outer iterative loop.
- `model_optimization_backend.hpp/.cpp`: GPU-aware model installation, likelihood evaluation, site batching, and branch-optimizer connection.
- `optimization_types.hpp`: shared constants, typed branch-update choices, options, and results used by placement, Local SPR, small-tip, and D&C.

## Key Numerical Functions

- `GlobalModelOptimizer::optimizeSubstitutionRates()`: searches for five relative GTR rates while keeping one reference rate fixed. It decides which candidate rates to test, while `ModelOptimizationBackend` installs each candidate and evaluates its likelihood.
- `optimizeFrequencies()`: optimize relative frequency parameters and normalize them to a valid probability vector.
- `optimizeAlpha()`: bounded one-dimensional optimization of Gamma alpha and regeneration of discrete rate categories.
- `optimize()`: runs the outer loop shown by the [responsibility boundary](#responsibility-boundary), alternating parameter groups and branch-length sweeps until the likelihood tolerance or maximum rounds is reached.

The internal `evaluateBounded()` and `evaluateOneDimensional()` functions are C-compatible adapters used by the numerical libraries. They convert raw parameter arrays/scalars back into the C++ objective stored in their context. They are not workflow entrypoints.

## Key Backend Functions

- `ModelOptimizationBackend::evaluate()`: implements the backend side of the [responsibility boundary](#responsibility-boundary) by evaluating the installed model using resident or site-batched CLVs.
- `installSubstitutionModel()`: applies one candidate model and follows the [PMAT rebuild dependency](../pmatrix/README.md#when-pmats-must-be-rebuilt): update frequencies/rates, rebuild the GTR Q matrix and eigendecomposition, upload model data, rebuild PMATs, and refresh CLVs when requested.
- `installGammaAlpha()`: rebuild Gamma rates and dependent GPU model state.
- `optimizeBranchLengths()`: dispatch to resident sequential or site-batched sequential full-tree branch optimization.
- `callbacks()`: bind these operations into the narrow callback interface used by `GlobalModelOptimizer`.

## Site-Batched Workspace

### Why Site Batching Exists

A resident likelihood evaluation keeps the CLVs for every alignment site on the
GPU at the same time. Its dominant CLV storage therefore grows approximately as

```text
number of CLV slots x sites x rate categories x states x sizeof(fp_t)
```

For a long alignment or a large tree, that allocation may not fit even though
the likelihood calculation itself is valid. Model optimization makes this
pressure more important because it evaluates many candidate rate, frequency,
and Gamma-alpha values and must repeatedly rebuild the likelihood state.

Site batching bounds that memory use by keeping CLVs for only `batch_sites`
contiguous sites in a reusable GPU workspace. It processes all batches and then
combines their contributions, so batching is a memory-execution policy rather
than a different statistical objective.

```text
full alignment:  [ sites 0 ................................ sites N-1 ]
                          |
                          v
GPU workspace:   [ one contiguous site batch ]
                          |
             rebuild batch CLVs and evaluate
                          |
                          v
CPU total:       batch 0 + batch 1 + ... + final batch
```

The transition matrices are indexed by branch and rate category, not by site,
so PMATs normally remain resident. The large site-dependent inputs and CLVs are
the data moved or rebuilt one batch at a time.

The batching fields in `ModelOptimizationContext` connect this policy to the
backend without exposing it to `GlobalModelOptimizer`:

- `site_batched` selects the resident or site-batched evaluation path.
- `site_batch_size` records how many sites the reusable workspace can hold.
- `site_batch_workspace` owns the temporary GPU buffers used for each range.

`GlobalModelOptimizer` still asks only to install a candidate or evaluate its
likelihood. `ModelOptimizationBackend` decides whether that request is executed
once over resident CLVs or repeatedly over site batches.

### Workspace Selection and Evaluation

- `chooseOptimizationSiteBatchSize()`: applies the [site-batched workspace policy](#site-batched-workspace): estimate bytes per site and select a batch that fits the target/free-memory budget.
- `RootSiteBatchWorkspace::evaluate()`: copy successive site ranges into a reusable upward-only workspace and sum their root likelihoods.
- `optimizeSequentialBranchNewtonSweeps()`: accumulate edge derivatives across site batches, update one edge, and continue in coordinate order.

`OptimizationWorkspaceKind` distinguishes global-parameter evaluation from the larger all-branch-gradient workspace.

## GPU Work and Cross-Site Accumulation

Site batching changes memory residency, not the objective. Each batch installs a contiguous site range, runs the same site-parallel CLV and likelihood kernels, and returns one batch log likelihood. The CPU adds batch results to obtain the full-alignment objective.

For sequential branch optimization, each site batch contributes a partial first derivative and partial second derivative for the current edge. Those batch contributions are accumulated before one Newton update is accepted, so a branch is never updated from only the last batch.

```text
for one edge
  batch 0 sites -> partial df, partial ddf
  batch 1 sites -> partial df, partial ddf
  ...
  CPU accumulation -> full-site df, full-site ddf
  one Newton update -> new branch length
```

Resident single-edge optimization can instead reduce sites entirely on the GPU. Threads first accumulate site subsets, each block emits a partial `(df, ddf)` pair, and the cooperative grid leader combines block partials before updating the branch. See [Existing-tree branch optimization](../placement/README.md#existing-tree-branch-optimization) for the block and reduction mapping and [How GPU Work Is Assigned](../pmatrix/README.md#how-gpu-work-is-assigned) for PMAT rebuild mapping.

## Typed Branch Options

- `EdgeUpdateScheme`: Jacobi or sequential/Gauss-Seidel updates.
- `ClvRetention`: rebuild all CLVs or preserve installed tip CLVs for virtual D&C boundary tips.
- `AcceptanceScope`: validate likelihood on a local subtree or the full tree.
- `BranchOptimizationOptions`: shared run-time policy.
- `BranchLengthOptimizationResult`: likelihood, timing, attempt, and acceptance report returned to workflows.

When changing acceptance tolerances or branch bounds, update `optimization_types.hpp` instead of introducing a local constant in a CUDA file. When changing the optimizer parameterization, consult the design origin recorded in `docs/reference/model_and_branch_optimization_reference.md`.
