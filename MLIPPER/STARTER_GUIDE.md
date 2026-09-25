# Understanding the MLIPPER Repository

This guide explains why MLIPPER is structured as it is, how its workflows move
through the codebase, and where to begin when making a change. For installation
and commands, use the [repository README](README.md).

## Motivation and Design

Large gene-tree workflows face a recurring tradeoff. Distance-based methods
can construct or extend trees quickly, but their topology and branch lengths
are not optimized under the evolutionary likelihood model used by later
phylogenetic analysis. A conventional whole-tree maximum-likelihood search can
perform that refinement, but repeatedly running it for many large gene trees
is expensive in both time and memory.

MLIPPER bridges these stages. It starts from a tree that is cheap to obtain—a
backbone with missing taxa or a tree constructed by DIPPER—and spends
likelihood computation where it is useful:

```text
fast distance/backbone construction
  -> likelihood-based query placement
  -> local or sector topology repair
  -> model and whole-tree branch-length optimization
```

Query taxa are placed and committed one at a time. After each small batch of
committed queries, MLIPPER applies Local SPR to refine the affected
neighborhoods. Large de novo trees are refined as bounded NNI sectors. GPU
conditional likelihood vectors (CLVs), transition matrices, and site batching
make this iterative workflow practical within one process.

MLIPPER uses fast methods to narrow the search and propose candidate changes,
but accepts topology changes only when they improve likelihood under the
supported GTR+Gamma model.

## Design Lineage and Related Papers

| Earlier work                                                                                                              | Idea used by MLIPPER                                                            | Difference in MLIPPER                                                                                          |
| ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| [DIPPER](https://doi.org/10.1101/2025.08.12.669583)                                                                        | GPU distance calculation, placement construction, and D&C starting trees        | DIPPER proposes the starting topology; MLIPPER evaluates and refines it with sequence likelihood               |
| [EPA-ng](https://doi.org/10.1093/sysbio/syy054)                                                                            | Evolutionary placement and likelihood weight ratios                             | MLIPPER can emit compatible output or commit queries and refine the topology                                   |
| [MAPLE](https://doi.org/10.1038/s41588-023-01368-0) and [UShER/matOptimize](https://doi.org/10.1093/bioinformatics/btac401) | Fast placement followed by local SPR-style improvement                          | MLIPPER uses dense GTR+Gamma likelihood CLVs for aligned gene data                                             |
| [Sectorial search](https://doi.org/10.1111/j.1096-0031.1999.tb00278.x)                                                     | Reoptimize bounded parts of a large tree                                        | MLIPPER preserves omitted context with directional boundary CLVs and accepts NNI changes by likelihood         |
| [RAxML-NG](https://doi.org/10.1093/bioinformatics/btz305) and libpll-style computation                                     | Felsenstein pruning, GTR parameterization, Gamma rates, and branch optimization | MLIPPER supplies its own CUDA execution and restricted model surface; it is not a general RAxML-NG replacement |

## Vocabulary

- **reference/backbone taxa:** taxa already present in the input tree;
- **query taxa:** taxa added by the small-tip workflow;
- **placement:** scoring candidate tree edges for a query;
- **commit:** attaching a query to its selected edge;
- **SPR:** subtree-prune-and-regraft, used after small-tip commits;
- **NNI:** nearest-neighbor interchange, used in D&C sectors;
- **CLV:** conditional likelihood vector stored for tree likelihood evaluation;
- **model optimization:** updating supported GTR+Gamma parameters;
- **branch optimization:** updating branch lengths without changing topology.

## Tree Representation

MLIPPER accepts a binary tree in rooted or unrooted Newick form. For example:

```text
((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);
```

libpll parses the Newick and validates that the topology is binary. If the
input is unrooted, libpll places a temporary root on one edge so MLIPPER can
orient every relationship as parent, left child, and right child. This root is
a computational device for traversals and likelihood evaluation; it does not
change the biological interpretation of the tree as unrooted.

MLIPPER then maintains three related representations:

1. `TreeBuildResult` is the authoritative mutable CPU topology. It stores one
   `TreeNode` per tip or internal node, the parent/child links, branch length to
   the parent, root ID, preorder and postorder traversals, and tip-name lookup.
2. `HostPacking` converts that topology and the alignment into contiguous host
   arrays, including links, branch lengths, tip states, traversal operations,
   transition matrices, and pattern weights.
3. `OwnedDeviceTree` owns the corresponding CUDA allocations. Its
   `DeviceTree` base is a non-owning descriptor passed to GPU kernels.

The representations flow in one direction:

```text
Newick + aligned tip names
  -> libpll parse and rooted traversal orientation
  -> TreeBuildResult (mutable CPU topology)
  -> HostPacking (contiguous evaluation arrays)
  -> OwnedDeviceTree / DeviceTree (GPU storage and kernel view)
```

Topology changes are made to `TreeBuildResult` first. MLIPPER then rebuilds the traversals and refreshes the derived host and GPU state before evaluating the changed tree. Existing full-tree node IDs remain array indices and are not renumbered by `rebuild_traversals()`; new query nodes are appended. D&C extraction creates a separate sector-local ID namespace, while a stable node label preserves logical identity when code must cross such rebuild or remapping boundaries. See [Node IDs, traversal order, and stable labels](src/tree/README.md#node-ids-traversal-order-and-stable-labels).

## Workflow Architecture

With the tree representations established, `src/main.cpp` is the executable
entrypoint and selects the workflow. `MlipperSession` owns the active CPU tree,
derived state, and GPU allocations for that workflow; it is a stateful engine,
not a second executable entrypoint.

### Small-tip commit

```text
parse the backbone, reference/query alignments, and GTR+Gamma model
  -> validate matching taxa and alignment lengths
  -> compress repeated columns across both alignments
  -> build MLIPPER's CPU tree and likelihood state
  -> place and commit each query individually
  -> after each small batch, refine affected neighborhoods with Local SPR
  -> optimize model parameters and full-tree branch lengths
  -> write the final Newick tree
```

Placement-only shares parsing, validation, and placement scoring with this
path, but writes jplace records without committing queries or refining the
tree.

### Divide and conquer

```text
parse the full alignment and GTR+Gamma model
  -> select and reserve a GPU
  -> use DIPPER to construct the starting topology
  -> compress repeated alignment columns
  -> build MLIPPER's CPU tree and GPU likelihood state
  -> partition internal edges into bounded sectors
  -> refine each sector with boundary-aware NNI
  -> optimize model parameters with site batching
  -> optimize branch lengths in sectors
  -> write the final Newick tree
```

DIPPER only proposes the starting topology. MLIPPER rebuilds aligned state,
partitions internal edges into sectors, preserves omitted context with
boundary CLVs, and evaluates NNI and optimization work itself.

Both workflows compress repeated alignment columns. Small-tip preprocessing
is performed while `loadSmallTipConfigFromCommandLine()` constructs its
validated config. D&C keeps the original full alignment long enough to build
the DIPPER starting tree, then preprocesses it in `runDivideAndConquerWorkflow()`
before initializing MLIPPER's likelihood state.

The corresponding implementation entry points are
`workflow::buildDipperStartingTree()`,
`MlipperSession::initializeDivideAndConquerGPU()`,
`MlipperSession::runDivideAndConquerNNI()`, and
`MlipperSession::runDivideAndConquerFinalOptimization()`. The orchestration is
kept in `runDivideAndConquerWorkflow()` in `src/main.cpp`.

## Repository Map

| Path                      | Responsibility                                                      |
| ------------------------- | ------------------------------------------------------------------- |
| `src/main.cpp`          | Workflow selection and orchestration                                |
| `src/mlipper_session.*` | Shared state, lifecycle, and high-level operations                  |
| `src/io/`               | CLI, validation, parsing, Newick, and jplace                        |
| `src/tree/`             | CPU tree representation, topology, packing, upload, and D&C sectors |
| `src/placement/`        | Candidate-edge scoring and placement branch lengths                 |
| `src/spr/`              | Local SPR and sector NNI                                            |
| `src/optimize/`         | Model and branch-length optimization                                |
| `src/likelihood/`       | CUDA likelihood kernels                                             |
| `src/pmatrix/`          | Transition matrices and eigendecomposition                          |
| `src/gpu/`              | GPU selection and multi-process admission                           |
| `src/workflow/`         | Narrow adapters to external workflow components                     |
| `tests/`                | CPU regressions and GPU integration fixture                         |
| `scripts/`              | Maintained ROADIES single-gene wrapper                              |
| `docker/`               | Development and runtime containers                                  |
| `docs/`                 | Design references and historical evidence                           |
| `third_party/dipper/`   | Vendored MLIPPER-compatible DIPPER source                           |

## Recommended Reading Order

Do not begin by reading `mlipper_session.cpp` from top to bottom. Trace one
workflow from its public boundary and descend only into the relevant layer:

1. Read [`src/README.md`](src/README.md) for the source map and call flows.
2. Read `src/main.cpp` for the complete workflow order.
3. Read [`src/io/README.md`](src/io/README.md) and `src/io/CLI.hpp` for typed configuration.
4. Read the public section of `src/mlipper_session.hpp` for engine operations and lifecycle.
5. Read [`src/tree/README.md`](src/tree/README.md) and `src/tree/tree.hpp` for tree and GPU-state vocabulary.
6. Search `src/mlipper_session.cpp` for the public operation being changed.
7. Enter `placement/`, `spr/`, `optimize/`, `likelihood/`, `gpu/`, or `pmatrix/` through its README.
8. Use [`docs/README.md`](docs/README.md) only when a change needs a design reference or historical experiment.

As a first exercise, trace `runSmallTipWorkflow()` to `runSmallTipBatches()`,
then find where a placement becomes a committed attachment and where Local SPR
is invoked. For D&C, trace `runDivideAndConquerWorkflow()` through DIPPER tree
construction, sector NNI, and final optimization.

## Where to Make a Change

| Change                                     | Start here                                                    |
| ------------------------------------------ | ------------------------------------------------------------- |
| CLI flags, defaults, or validation         | `src/io/CLI.cpp`, `src/io/input_validation.*`             |
| Workflow order or output dispatch          | `src/main.cpp`                                              |
| Shared state, ownership, and lifecycle     | `src/mlipper_session.*`                                     |
| Alignment parsing or column compression    | `src/io/parse_file.*`, `src/util/msa_preprocess.*`        |
| Tree parsing, topology, packing, or upload | `src/tree/`                                                 |
| Candidate scoring and placement lengths    | `src/placement/`                                            |
| Local SPR                                  | `src/spr/local_spr*`                                        |
| Sector NNI                                 | `src/spr/nni.*`, `src/tree/divide_and_conquer.*`          |
| Model and branch optimization              | `src/optimize/`, `src/mlipper_session.cpp`                |
| Likelihood kernels                         | `src/likelihood/`                                           |
| Transition matrices/eigendecomposition     | `src/pmatrix/`                                              |
| Newick or jplace output                    | `src/io/tree_newick.*`, `src/io/jplace.*`                 |
| DIPPER adapter                             | `src/workflow/dipper_starting_tree.*`                       |
| GPU admission                              | `src/gpu/gpu_admission.*`                                   |
| ROADIES invocation                         | `README_ROADIES.md`                                      |

## Development Principles

For a first contribution, prefer a parser/validation fix or CPU regression
test. Changes to placement scoring, CLV/scaler layout, topology commits,
boundary messages, or optimization acceptance require GPU regression of both
workflows.

For a normal change:

1. Identify the affected workflow and owning layer.
2. Make the smallest coherent change.
3. Build and run the CPU tests described in the [README](README.md).
4. Run the GPU integration fixture when CUDA is available.
5. Check the relevant `--help` output after CLI changes.
6. Compare the affected real-data workflow with an approved golden.
7. Regress both workflows after shared session, tree, likelihood, model, PMAT,
   or GPU-resource changes.

For algorithm changes, record the command and commit, dataset identity,
precision and visible GPU, output artifact, initial/final likelihood where
applicable, topology comparison against a named golden, runtime, and peak GPU
memory when performance is relevant. Do not call a result “golden” without
recording its reference tree and normalization procedure.

## Current Boundaries

- Model support is restricted to four-state DNA GTR, 1–8 equal-weight
  discrete-Gamma categories, and `pinv = 0`.
- The small-tip helper assumes split reference/query alignments.
- One CUDA device is selected per invocation.
- Real-data scientific regressions and approved goldens are external.
- DIPPER updates must retain the integration changes documented in
  [`third_party/dipper/README.md`](third_party/dipper/README.md).
- Documents under `docs/experiments/` record historical experiments and are
  not the source of truth for the current CLI.
