# MLIPPER Documentation Map

This directory holds design references and historical evidence. It is not the starting point for building or running MLIPPER.

## New Maintainer: Read in This Order

1. Read the root [`README.md`](../README.md) for the project overview, dependencies, build and test commands, and complete workflow examples.
2. Use [`STARTER_GUIDE.md`](../STARTER_GUIDE.md) for the design motivation, architecture, repository map, and contributor reading order.
3. Read [`src/README.md`](../src/README.md) for the small-tip and D&C call paths and the source-directory ownership map.
4. Read the README in the component you will change, such as [`src/optimize/README.md`](../src/optimize/README.md) or [`src/tree/README.md`](../src/tree/README.md).
5. Use the documents below only when the change needs their specific design or historical context.

The source code and public CLI are authoritative. A recorded experiment does not define current behavior.

## Current Design References

These documents describe constraints that still affect implementation work.

| Document | Read when | Status |
| --- | --- | --- |
| [`reference/supported_model_scope.md`](reference/supported_model_scope.md) | Changing model parsing, rate categories, or likelihood model support | Current supported boundary |
| [`reference/model_and_branch_optimization_reference.md`](reference/model_and_branch_optimization_reference.md) | Understanding why the GTR, Gamma, and branch optimizers use their current parameterization and initial constants | Historical design reference; current differences are explicit |
| [`D&C_SPR_BUILD_PLAN.md`](D&C_SPR_BUILD_PLAN.md) | Designing the proposed whole-tree radius-based SPR stage after sector NNI | Proposed work; not implemented and not a current runtime promise |

The optimizer reference explains which RAxML-NG 0.9.0 design choices MLIPPER retained. Current code and `src/optimize/README.md` remain authoritative; the reference does not promise exact stopping-rule or bit-for-bit compatibility.

## Performance Evidence

| Document | Question it answers | Status |
| --- | --- | --- |
| [`performance/full_tree_blo_optimization_profile.md`](performance/full_tree_blo_optimization_profile.md) | Where time was spent during sequential full-tree branch optimization and why the present kernel path was selected | Historical measurement on gene638 |

Performance reports describe one recorded machine, dataset, and revision. Run a new benchmark before making a current performance claim.

## Historical Experiments

| Document | Question it investigated | Status |
| --- | --- | --- |
| [`experiments/alisim10k.md`](experiments/alisim10k.md) | Which placement, Local SPR, sector NNI, and branch/model experiments changed AliSim topology error | Historical experiment catalog |
| [`experiments/rnasim10k.md`](experiments/rnasim10k.md) | Whether conclusions from AliSim transferred to a gappy RNASim 10K dataset | Historical cross-dataset record |

These files intentionally preserve old paths, commands, negative results, and intermediate designs for provenance. Do not copy their commands into a new workflow without checking the current root README and `MLIPPER --help`.

## Completed Maintenance Records

| Document | What it records | Status |
| --- | --- | --- |
| [`code_commenting_review_plan.md`](code_commenting_review_plan.md) | Architecture-first source commenting/style review, deferred findings, and verification evidence | Completed maintenance record; not a current implementation specification |

## Where New Documentation Belongs

- Project identity and public build, CLI, input, or workflow behavior: root `README.md`.
- Design overview, repository orientation, and contributor workflow: `STARTER_GUIDE.md`.
- Source ownership and call paths: `src/README.md` or the owning component README.
- Stable external algorithm compatibility: `docs/reference/`.
- Reproducible profiling that explains a design decision: `docs/performance/`.
- Completed experiment evidence: `docs/experiments/`.

Every new performance or experiment record must identify the MLIPPER commit, DIPPER source revision, dataset, precision, GPU, exact command, result, and conclusion. Large generated trees, logs, profiles, and datasets stay outside the repository.
