# Source Guide

This directory contains the MLIPPER executable and its maintained C++/CUDA implementation. Start here after reading the repository-level `README.md`. For a higher-level explanation of the tree representations and workflow design, read [`STARTER_GUIDE.md`](../STARTER_GUIDE.md) first.

## Recommended Reading Order

1. `main.cpp` selects the small-tip or divide-and-conquer workflow.
2. `io/README.md` explains how command-line input becomes validated config.
3. `mlipper_session.hpp` shows the workflow-level API and long-lived state.
4. `tree/README.md` explains the CPU tree, packed arrays, and GPU tree views.
5. Read the README for the algorithm you are changing: `placement/`, `spr/`, `likelihood/`, or `optimize/`.

Do not begin with a large `.cu` file. First identify the public host function that launches its kernels and the workflow that calls that host function.

## End-to-End Call Flow

Small-tip placement and commit:

```text
main.cpp
  -> cli::loadSmallTipConfigFromCommandLine()
       -> parse and validate the inputs
       -> preprocess_alignments() for reference and query alignments
  -> MlipperSession::load*()
  -> MlipperSession::initializeCPU()
  -> MlipperSession::runSmallTipBatches()
       -> place and commit queries individually
       -> Local SPR after each completed batch
  -> MlipperSession::runFinalModelOptimization()
  -> MlipperSession::writeTree()
```

Divide and conquer:

```text
main.cpp
  -> cli::loadDivideAndConquerConfigFromCommandLine()
  -> gpu::select_device_or_wait_or_throw()
  -> workflow::buildDipperStartingTree()
  -> preprocess_alignments() for the full alignment
  -> MlipperSession::load*() / initializeCPU()
  -> MlipperSession::initializeDivideAndConquerGPU()
  -> MlipperSession::runDivideAndConquerNNI()
  -> MlipperSession::runDivideAndConquerFinalOptimization()
  -> MlipperSession::writeTree()
```

## Root Files

- `main.cpp`: the only executable entrypoint. It should coordinate workflows, not implement likelihood or topology algorithms.
- `mlipper_session.hpp/.cpp`: owns one run's CPU tree, model, GPU buffers, placement state, topology-refinement workspaces, and GPU reservation.
- `util/model_utils.hpp/.cpp`: parses RAxML-NG `bestModel` input, normalizes base frequencies, estimates empirical frequencies, and builds Gamma categories.
- `util/msa_preprocess.hpp/.cpp`: validates compatible alignments and compresses repeated site patterns into one column plus a pattern weight.

## Key `MlipperSession` Functions

- `loadBackboneTree()`, `loadAlignment()`, `loadModel()`: install validated inputs from the [I/O validation boundary](io/README.md#data-structures).
- `initializeCPU()`: parse the Newick tree, build host packing, construct the GTR eigendecomposition, and create the [CPU/host tree representations](tree/README.md#the-short-version).
- `initializeGPU()`: create the resident [`OwnedDeviceTree`](tree/README.md#devicetree-and-owneddevicetree) for placement.
- `runSmallTipBatches()`: place and commit queries individually, rebuild the required GPU state, and run the [Local SPR round](spr/README.md#round-structure) after each completed query batch.
- `runDivideAndConquerNNI()`: create [bounded sectors with boundary messages](tree/README.md#dc-boundary-tips) and run NNI refinement.
- `runFinalModelOptimization()`: alternate model-parameter and global branch length optimization using the [optimizer/backend separation](optimize/README.md#responsibility-boundary).
- `runDivideAndConquerFinalOptimization()`: run the corresponding site-batched and sector-partitioned final optimization for D&C.
- `writeTree()` and `writeJplace()`: serialize final workflow results.

## Directory Responsibilities

- `gpu/`: GPU selection, cross-process admission, and small RAII buffers.
- `io/`: CLI parsing, validation, alignment input, Newick, and jplace output.
- `likelihood/`: partial/CLV kernels and root or placement likelihoods.
- `optimize/`: numerical model optimization and branch-optimization policy.
- `placement/`: candidate-edge scoring and branch-length derivatives.
- `pmatrix/`: GTR eigendecomposition and transition-matrix construction.
- `spr/`: shared topology-refinement engine for Local SPR and sector NNI.
- `tree/`: tree representation, packing, device allocation, CLV schedules, subtree construction, and D&C boundary messages.
- `util/`: precision primitives and legacy shared CUDA utility structures.
- `workflow/`: adapter for DIPPER starting-tree construction.

CUDA work is documented beside the data it produces: [PMAT construction](pmatrix/README.md#how-gpu-work-is-assigned), [CLV and likelihood kernels](likelihood/README.md#gpu-parallelism-and-outputs), [placement and derivative reductions](placement/README.md#gpu-parallelism-during-placement), and [tree traversal assignment](tree/README.md#how-tree-work-is-assigned-on-the-gpu).

## State and Ownership Rules

- `TreeBuildResult` owns the mutable CPU topology.
- `HostPacking` owns host arrays derived from that topology.
- `OwnedDeviceTree` owns CUDA allocations; `DeviceTree` is only a non-owning, trivially-copyable view passed to kernels.
- `PlacementOpBuffer` owns/prepares traversal operations but does not own the tree itself.
- `gpu::DeviceReservation` owns one shared admission-table entry, not CUDA memory. Its destructor releases that reservation.
- After changing topology, rebuild traversals and refresh every affected host and device representation before evaluating likelihood.

## Before Changing Shared Code

Run `make test`. Changes to tree packing, likelihood, PMAT, optimization, session state, or GPU resource lifetime also require the opt-in GPU workflow regression described in `tests/README.md`.
