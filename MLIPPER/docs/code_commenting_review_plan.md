# Code Commenting and Style Review Plan

## Purpose

Make the end-to-end design flow understandable before polishing individual
statements. Review MLIPPER-owned production code line by line, normalize locally
inconsistent style, and add comments where they preserve information that cannot
be recovered by simply reading the statement itself. This review does not include
vendored code under `third_party/` and must not mix in unrelated algorithm changes.

## Architecture-first review

Follow control and data in execution order rather than reviewing directories in
alphabetical order:

```text
CLI parsing and validation
  -> workflow selection in main
  -> MlipperSession input/state transitions
  -> CPU tree/model construction
  -> host packing and GPU ownership
  -> placement and likelihood evaluation
  -> topology commit/refinement
  -> model and branch optimization
  -> Newick or jplace serialization
```

At every transition, establish which object owns the data, which representations
are derived caches, what invalidates them, and whether execution crosses a
CPU/GPU or external-library boundary. Add comments at those boundaries; local
implementation comments come afterward.

## Commenting standard

Add a comment when it explains at least one of the following:

- an invariant, ownership rule, or resource lifetime;
- why an apparently simpler implementation is incorrect;
- a CPU/GPU synchronization or memory-layout contract;
- a numerical-stability constraint or domain-specific convention;
- a non-obvious error-handling or concurrency decision;
- a public API's preconditions, postconditions, or units.

Do not add comments that:

- restate the next line of code;
- narrate loops, assignments, or obvious control flow;
- preserve dead code or commented-out implementations;
- make promises that are not enforced by code or tests;
- describe temporary implementation details likely to drift.

Comments should sit immediately above the contract they describe. Use complete
sentences for design comments and short phrases only for field/layout labels.
Public interfaces use `///` only when documenting an externally useful contract;
implementation notes use `//`.

## Style review standard

Because the repository has no checked-in formatter or linter configuration,
preserve the dominant local style rather than performing a repository-wide
reformat. During each file review, check:

- include grouping and unused includes;
- naming and namespace consistency;
- braces, indentation, line length, and one-statement-per-line consistency;
- signed/unsigned conversions and narrowing at API/kernel boundaries;
- constness, `noexcept`, ownership, and exception safety;
- duplicated literals and unexplained units;
- stale, misleading, redundant, and commented-out code;
- error messages that omit the failing operation or relevant dimensions;
- CUDA launch, synchronization, allocation, and host/device layout contracts.

Behavioral changes discovered during review must be small, directly justified,
and covered by an existing or new test. Larger defects should be recorded for a
separate change instead of being hidden in comment-only work.

## Deferred correctness findings

- `commit_placement_result()` mutates `TreeBuildResult` before rebuilding
  `HostPacking`, PMATs, device slices, and CLVs. An exception after insertion can
  leave those representations out of sync. A separate change should define a
  rollback/commit transaction and fault-injection tests; comments must not claim
  strong exception safety until that exists.
- `ensure_local_spr_placement_scratch_capacity()` releases the old scratch and
  then performs several raw `cudaMalloc` calls. A failure partway through can
  leave a partially allocated workspace. Converting the scratch fields to an
  owning aggregate constructed off to the side would provide a transactional
  replacement; comments must not claim that guarantee with the current layout.
- `RootSiteBatchWorkspace::ensureCapacity()` publishes capacity metadata and
  performs several raw `cudaMalloc` calls directly into the live workspace.
  A failed allocation can leave partial storage. An owning candidate workspace
  should be fully constructed before replacing the live batch state.
- Model installation validates and constructs the CPU candidate before commit,
  but uploads several independent GPU arrays afterward. A CUDA failure between
  uploads can still leave CPU and GPU model representations out of sync; full
  transactionality requires double-buffered device model storage or rollback.

## Review order and progress

Each file moves through `pending -> reviewed -> verified`. A directory is only
complete after all of its production headers and sources are verified.

| Area | Scope | Status |
| --- | --- | --- |
| Utilities | `src/util/*` | Verified |
| GPU resource/admission | `src/gpu/*` | Verified |
| I/O and CLI | `src/io/*` | Verified |
| Model and preprocessing | `src/util/model_utils.*`, `src/util/msa_preprocess.*` | Verified |
| Tree and divide-and-conquer | `src/tree/*` | Verified |
| Probability matrices | `src/pmatrix/*` | Verified |
| Likelihood kernels | `src/likelihood/*` | Verified |
| Placement kernels | `src/placement/*` | Verified |
| Topology refinement | `src/spr/*` | Verified |
| Optimization | `src/optimize/*` | Verified |
| Session and executable | `src/mlipper_session.*`, `src/main.cpp` | Verified |
| Workflow | `src/workflow/*` | Verified |

## Per-file procedure

1. Read the header and implementation together, including direct callers.
2. Identify public contracts, hidden invariants, resource ownership, and device
   memory layouts before editing comments.
3. Remove or rewrite misleading comments before adding new ones.
4. Apply only locally consistent style and low-risk correctness improvements.
5. Re-read the diff without surrounding implementation context; every added
   comment must still be accurate and useful.
6. Compile the affected translation units and run the narrowest relevant tests.
7. Mark the area verified only after its complete file set has passed review.

## Completion gates

- Every MLIPPER-owned production file in the scope table has been reviewed.
- No stale TODO, FIXME, dead-code comment, or unexplained public contract remains
  without an explicit decision.
- `git diff --check` passes.
- `make test` passes.
- A full `make MLIPPER` build passes.
- The final diff is inspected separately for comment accuracy and accidental
  behavior changes.

## Completion evidence

- All production areas in the review matrix reached `Verified`.
- The final owned-source scan found no unresolved `TODO`, `FIXME`, `XXX`, or
  `HACK` markers.
- `git diff --check`, `make test -j2`, and `make MLIPPER -j2` passed.
- The final comment audit removed a misleading strong-exception-safety claim
  from placement commit documentation; remaining transactional limitations are
  recorded under deferred correctness findings.
