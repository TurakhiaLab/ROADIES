# Regression Tests

A regression test fixes an already approved behavior and reports when a later code change alters it. Passing a regression does not prove that the algorithm is scientifically correct; it proves that the tested behavior did not change unexpectedly.

## Test Map

- `msa_preprocess_test.cpp`: repeated-column compression, pattern weights, and reference/query alignment consistency.
- `model_utils_test.cpp`: model parsing, normalization, empirical frequencies, and Gamma-category helpers.
- `input_validation_test.cpp`: invalid names, symbols, model shapes, paths, and CLI option combinations.
- `tree_io_test.cpp`: Newick parsing/serialization and jplace-facing tree I/O.
- `topology_refinement_test.cpp`: topology-only SPR/NNI legality, traversal, and candidate-selection behavior without GPU scoring.
- `gpu_admission_test.cpp`: estimator behavior and move-only reservation API; it links CUDA/NVML but does not initialize a GPU.
- `compare_workflow_outputs_test.py`: positive and negative tests for tolerant jplace values, exact edge assignments, canonical tree topology, and branch lengths.
- `run_workflow_regression.sh`: opt-in placement and D&C end-to-end GPU check.
- `test_gpu_auto_contention.sh`: opt-in multi-process host test for automatic GPU selection and shared reservation behavior.

When fixing a bug, place the smallest reproduction in the test matching the owning component. Add an end-to-end case only when the behavior cannot be verified below the workflow boundary.

## CPU Tests

Run the host-side unit tests with:

```sh
make test
```

This also checks the topology-only SPR/NNI primitives and GPU-admission memory accounting. It links the CUDA/NVML libraries but does not initialize a GPU or run a dataset benchmark.

Use these tests while developing. Add a focused test here whenever the bug or behavior can be reproduced without GPU execution.

The `Makefile` owns the list of binaries run by `make test`. Creating a test source file alone does not add it to the suite.

## Small GPU Integration Regression

The repository includes a six-taxon fixture under `tests/fixtures/gpu_small`. It runs placement-only and D&C through the production executable. The D&C case covers every internal edge and accepts an NNI move, then runs final model and branch optimization.

Run it on an available CUDA GPU with:

```sh
make test-gpu-integration
```

The comparison ignores only jplace `metadata.invocation`. Query ordering, candidate edge assignments, and D&C topology must match. Floating-point likelihoods and branch lengths use relative tolerance `1e-7` and absolute tolerance `1e-6`; comparator unit tests verify that meaningful edge, topology, and branch-length changes fail.

This fixture detects integration regressions but is too small to establish scientific accuracy, large-tree memory behavior, or performance.

## External GPU Workflow Regression

For real-data sign-off, override every fixture path with a provenance-recorded dataset and reviewed expected outputs, then run `make test-workflows`:

```sh
export MLIPPER_TEST_REF_MSA=/path/to/reference.fa
export MLIPPER_TEST_QUERY_MSA=/path/to/query.fa
export MLIPPER_TEST_BACKBONE_TREE=/path/to/backbone.nwk
export MLIPPER_TEST_MODEL=/path/to/model.bestModel
export MLIPPER_TEST_EXPECTED_JPLACE=/path/to/golden.jplace
export MLIPPER_TEST_DNC_MSA=/path/to/full_alignment.fa
export MLIPPER_TEST_EXPECTED_DNC_TREE=/path/to/golden.dnc.nwk
export MLIPPER_TEST_GPU_ID=0
make test-workflows
```

Optional D&C overrides are `MLIPPER_TEST_DIPPER_TREE_MODE`, `MLIPPER_TEST_DNC_CORE_EDGES`, `MLIPPER_TEST_DNC_TIP_BUDGET`, and `MLIPPER_TEST_DNC_SWEEPS`.

`MLIPPER_TEST_REF_MSA` is the placement backbone alignment; `MLIPPER_TEST_DNC_MSA` is the separate full-taxon alignment used to construct the D&C tree. Do not point both variables at a split reference-only alignment.

## Golden-Output Policy

A golden is an approved expected artifact, not simply the newest output. Store or record the following together:

- MLIPPER commit and compatible DIPPER commit
- exact reference/query alignment, backbone tree, and model identities
- workflow command and non-default flags
- precision mode and GPU model
- expected jplace or D&C tree
- for scientific sign-off, the named reference species/gene tree and the tool and options used to calculate RF/nRF

Do not overwrite a golden merely to make a failing test pass. First classify the difference as one of:

- formatting-only and intentionally accepted
- expected numerical drift with unchanged selected edges/topology
- intended algorithm change with reviewed evidence
- unintended regression

Only the first three cases justify approving a replacement. Keep the old and new comparison in the experiment record that explains the decision.

The script always runs both placement-only and D&C, so use a complete set of overrides for shared-code sign-off. A workflow-specific local command may be used during development, but it does not replace the combined regression before merging a change to shared tree, likelihood, PMAT, model, session, or GPU-resource code.

To exercise cross-process GPU admission independently, build MLIPPER and run:

```sh
tests/test_gpu_auto_contention.sh --gene-dir /path/to/gene-inputs
```

Use more workers than visible GPUs when the purpose is to observe contention.

## Reading Failures

- Assertion failure in one unit binary: inspect that component before running a full GPU workflow.
- jplace mismatch: compare selected edge, optimized lengths, likelihood/LWR, and metadata separately; the script ignores only invocation text.
- Newick comparison failure: the comparator reports taxon-set, canonical topology, or branch-length mismatch separately.
- CUDA failure: retain the first CUDA error; later cleanup errors are often a consequence rather than the cause.
