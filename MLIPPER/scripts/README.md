# Reproducibility Scripts

This directory contains experiment drivers. General multi-gene scheduling,
species-tree analysis, and downstream metrics belong in ROADIES or an
experiment repository.

## ROADIES yeast iteration reproducibility

`run_roadies_yeast_iterations.py` runs complete extracted yeast iterations with one sequential worker per GPU. By default it runs the remaining Iterations 2–5 on GPUs 0 and 1 with the established small-tip configuration:

```sh
scripts/run_roadies_yeast_iterations.py --dry-run
scripts/run_roadies_yeast_iterations.py
```

Each iteration receives a manifest, per-gene tree and log, failure table, aggregate tree file, checksums, and a JSON summary. Existing non-empty trees are reused by default, so the same command resumes an interrupted run. Use `--no-resume` to recompute everything.

Byte-for-byte equality requires an existing MLIPPER reference tree for every gene. Point `--reference-root` at this layout:

```text
REFERENCE_ROOT/
  iter_2/genes/gene_N/mlipper_gene_tree.nwk
  iter_3/genes/gene_N/mlipper_gene_tree.nwk
  ...
```

When a reference root is supplied, the driver writes `byte_comparison.tsv` and records matches, mismatches, and missing references separately. Without one, the summary says `not_checked`; successful execution is never presented as byte identity.
