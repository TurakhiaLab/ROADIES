# Small GPU Integration Fixture

This six-taxon fixture checks that both maintained GPU workflows cross their full executable boundary. It is deliberately small enough to run on every developer change; it is not a scientific accuracy benchmark.

The placement case uses four reference taxa (`A`--`D`) and two query taxa (`E`, `F`). The D&C case uses the corresponding six-taxon `full.fa`. Its baseline run covers all three internal edges and accepts one NNI move, so the test exercises candidate scoring, exact likelihood validation, topology replay, and final optimization rather than only CUDA startup.

The checked-in expected outputs were generated from MLIPPER commit `13e57c5` plus the source-only `model_utils`/`msa_preprocess` directory move, using the default double-precision build on an NVIDIA RTX A6000. The model and all input files are stored beside the outputs.

`tests/compare_workflow_outputs.py` requires jplace structure, query ordering, candidate edge assignments, and tree topology to match. It compares floating values with relative tolerance `1e-7` and absolute tolerance `1e-6`, and ignores only `metadata.invocation`. This accommodates insignificant GPU arithmetic drift without accepting a different placement or topology.

Do not replace these expected files merely because the regression fails. Classify and review the difference first. Real-data scientific sign-off remains external because this fixture is too small to establish RF/nRF quality or large-tree performance.
