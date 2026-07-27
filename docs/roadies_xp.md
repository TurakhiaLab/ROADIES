# ROADIES_XP: Placement Mode & GPU Acceleration

!!! Note
    A manuscript describing ROADIES_XP is in preparation as a publication separate from (but building on) the original [ROADIES PNAS paper](https://doi.org/10.1073/pnas.2500553122).

## Introduction

ROADIES_XP extends the de novo ROADIES pipeline described elsewhere in this Wiki with two independent, opt-in capabilities:

- **Placement mode**: grow or update an existing ("backbone") species tree with new query genomes, without re-inferring the whole tree from scratch.
- **GPU acceleration**: swap in GPU-accelerated tools for the alignment and placement stages, usable in both de novo and placement modes.

De novo mode (`accurate` / `balanced` / `fast`, CPU-only) remains ROADIES' default behavior. Nothing about the pipeline you already know changes unless you opt into `--mode placement` and/or `--gpu` explicitly.

## Tool substitutions by mode and backend

| Stage | De novo (CPU) | De novo (GPU) | Placement (CPU) | Placement (GPU) |
| --- | --- | --- | --- | --- |
| Pairwise alignment | [LASTZ](https://lastz.github.io/lastz/) | [KegAlign](https://github.com/galaxyproject/KegAlign) | LASTZ | KegAlign |
| Multiple sequence alignment | [PASTA](https://github.com/smirarab/pasta) | [TWILIGHT](https://github.com/TurakhiaLab/TWILIGHT) | TWILIGHT (aligns query onto backbone MSA/tree) | TWILIGHT (aligns query onto backbone MSA/tree) |
| Gene tree building | RAxML-NG (unconstrained) | RAxML-NG (unconstrained) | RAxML-NG, `--tree-constraint` on the backbone gene tree | [MLIPPER](https://github.com/TurakhiaLab/MLIPPER) (GPU-accelerated placement onto the backbone gene tree) |
| Species tree estimation | ASTRAL-Pro3 | ASTRAL-Pro3 | ASTRAL-Pro3 (optionally `--constraint` on the backbone species tree, see `--grow` below) | ASTRAL-Pro3 (optionally `--constraint`) |

!!! Note
    `--gpu` is supported with `--mode accurate` (default), `--mode balanced`, and `--mode placement`. `--mode fast` (MashTree-based) has no GPU variant.

## Placement mode

Placement mode aligns and places a directory of **query** genomes onto an existing **backbone** species/gene tree, instead of building a tree from scratch. It's the mode to use when you already have a ROADIES/ROADIES_XP species tree and want to add more genomes to it cheaply.

To run it directly with `run_roadies.py`, set two things in `config.yaml` first:

- `GENOMES`: path to the directory of **query** genome assemblies (same `.fa`/`.fa.gz` naming rules as de novo mode).
- `REF_DIR`: path to the **backbone**'s output directory — i.e. the `OUT_DIR` from a previous de novo (or placement) run of ROADIES on the reference genomes. This must contain that run's `genes/` alignments, gene trees, and models.

Then run:

```bash
python run_roadies.py --cores 16 --mode placement
```

The resulting species tree (`roadies.nwk`) is saved in the current run's `OUT_DIR`, same as de novo mode.

### Grow vs. update

By default (no `--grow`), placement mode re-infers the combined species tree freely from the backbone's and query's gene trees together — the final topology can revise relationships from the original backbone tree ("**update**").

Adding `--grow` instead constrains ASTRAL-Pro3's search to the backbone tree's topology, so the backbone relationships are preserved exactly and only the new query taxa are attached ("**grow**"):

```bash
python run_roadies.py --cores 16 --mode placement --grow
```

## GPU acceleration

Add `--gpu N` (where `N` is the number of GPU devices to use) to any supported mode:

```bash
python run_roadies.py --cores 16 --gpu 1                       # de novo, GPU-accelerated
python run_roadies.py --cores 16 --mode placement --gpu 1       # placement, GPU-accelerated
```

GPU mode requires the GPU-specific tools to be built, which `roadies_env.sh` handles automatically when the necessary build dependencies are present on your system:

- **TWILIGHT**: built automatically; uses CUDA if `nvcc` is available, falls back to a CPU build otherwise.
- **MLIPPER**: only built if both `nvcc` (CUDA) and `libpll` (`pll.h`) are found. If they're missing at setup time, `roadies_env.sh` skips the build with a warning — GPU *placement* specifically needs MLIPPER, so build it once CUDA/libpll are available:
  ```bash
  bash MLIPPER/install/setup_host.sh
  ```
- **KegAlign**: installed via the `kegalign` conda environment (`workflow/envs/kegalign.yaml`), used automatically by Snakemake's `--use-conda` for the `kegalign` rule.

## Iterative backbone + placement convergence

`placement_converge.py` automates repeated backbone-then-placement runs until the tree stabilizes — analogous to ROADIES' own de novo convergence loop, but across separate backbone and query genome sets. Each iteration: (1) builds/reuses a de novo backbone tree from the backbone genomes, (2) places the query genomes onto it, (3) combines gene trees across all iterations so far into a running ASTRAL-Pro3 estimate, and (4) checks the same high-support convergence criterion used by `converge.py`.

```bash
python placement_converge.py \
    --backbone /path/to/backbone_genomes \
    --query /path/to/query_genomes \
    --cores 16 \
    --gpu 1
```

Key arguments (run `python placement_converge.py --help` for the full list):

| Argument | Description | Default |
| --- | --- | --- |
| `--backbone` | Directory of backbone/reference genome assemblies. | *(required)* |
| `--query` | Directory of query genome assemblies to place onto the backbone. | *(required)* |
| `--config` | Path to `config.yaml`; will be updated in place with `OUT_DIR`/`GENOMES`/`REF_DIR` for each iteration. | `config/config.yaml` next to the script |
| `--out-base-dir` | Base directory where `iter_<n>_backbone/`, `iter_<n>_placement/`, and the final combined `roadies_final/` outputs are written. | `roadies_iterations/` next to the script |
| `--cores` | Number of CPU cores. | 48 |
| `--gpu` | Number of GPU devices to use. | 0 (CPU) |
| `--support-threshold` | Local posterior probability threshold for convergence. | 0.95 |
| `--max-iterations` | Maximum number of backbone+placement iterations. | 9 |
| `--resume` | Resume from the last completed iteration (reads `time_stamps.csv` in `--out-base-dir`). | off |
| `--deep` | Enable deep-phylogeny mode for both backbone and placement stages. | off |

The final converged tree is written to `<out-base-dir>/roadies_final/roadies.nwk` (and `roadies_stats.nwk` with support values), alongside per-iteration outputs in `<out-base-dir>/iter_<n>_backbone/` and `<out-base-dir>/iter_<n>_placement/`.

## New `config.yaml` parameters

These parameters are used by ROADIES_XP in addition to the ones documented in the [User Guide](usage.md):

| Parameter | Description | Default |
| --- | --- | --- |
| **REF_DIR** | Path to the backbone iteration's output directory, used only in `--mode placement`. Not used in de novo modes. | `null` |
| **GROUP_CSV** | Optional path to a CSV file with `species,group` columns, used to bias gene sampling towards under-represented lineages/clades rather than sampling genomes uniformly. Leave unset (or point to a nonexistent path) to sample uniformly, as in de novo mode. | `""` (unset) |
