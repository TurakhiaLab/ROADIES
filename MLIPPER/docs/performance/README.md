# Performance Evidence

Documents here preserve profiling results that explain an implementation decision. They are tied to the recorded dataset, hardware, precision, and Git revision and are not current performance guarantees.

Before replacing a result, retain the old and new measurements and explain the code change that caused the difference.

Current record:

- [`full_tree_blo_optimization_profile.md`](full_tree_blo_optimization_profile.md) compares the single-block and cooperative multi-block sequential branch-optimization paths on gene638.
