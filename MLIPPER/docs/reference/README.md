# Design References

Documents here record stable implementation constraints or compatibility decisions. Read the owning `src/*/README.md` first, then use these references to understand why the current boundary or algorithm was selected.

Update a reference only when the accepted design changes. Ordinary benchmark results belong under `docs/performance/` or `docs/experiments/`.

Current files:

- [`supported_model_scope.md`](supported_model_scope.md) defines the accepted public model boundary.
- [`model_and_branch_optimization_reference.md`](model_and_branch_optimization_reference.md) records the RAxML-NG/Coraxlib design lineage and current deviations.
