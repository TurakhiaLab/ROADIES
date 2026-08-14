# ROADIES_XP: Placement mode with GPU Acceleration

!!! Note
    A manuscript describing ROADIES_XP is in preparation as a publication separate from (but building on) the original [ROADIES PNAS paper](https://doi.org/10.1073/pnas.2500553122).

## Introduction

ROADIES_XP extends the de novo ROADIES pipeline described elsewhere in this Wiki with **placement mode**: grow or update an existing ("backbone") species tree with new query genomes, without re-inferring the whole tree from scratch.

Placement mode runs in two variants:

- **CPU placement** (default)
- **GPU placement** (add `--gpu N`) — swaps in GPU-accelerated tools for the alignment and placement stages

De novo mode (`accurate` / `balanced` / `fast`) remains ROADIES' default behavior and is CPU-only — it has no GPU variant. Nothing about the pipeline you already know changes unless you opt into `--mode placement` explicitly.

## Tool substitutions by backend

| Stage | De novo (CPU) | Placement (CPU) | Placement (GPU) |
| --- | --- | --- | --- |
| Pairwise alignment | [LASTZ](https://lastz.github.io/lastz/) | LASTZ | [KegAlign](https://github.com/galaxyproject/KegAlign) |
| Multiple sequence alignment | [PASTA](https://github.com/smirarab/pasta) | [TWILIGHT](https://github.com/TurakhiaLab/TWILIGHT) (aligns query onto backbone MSA/tree) | TWILIGHT (aligns query onto backbone MSA/tree) |
| Gene tree building | RAxML-NG (unconstrained) | RAxML-NG, `--tree-constraint` on the backbone gene tree | [MLIPPER](https://github.com/TurakhiaLab/MLIPPER) (GPU-accelerated placement onto the backbone gene tree) |
| Species tree estimation | ASTRAL-Pro3 | ASTRAL-Pro3 (optionally `--constraint` on the backbone species tree, see `--grow` below) | ASTRAL-Pro3 (optionally `--constraint`) |

!!! Note
    `--gpu` is only supported with `--mode placement`. De novo modes (`accurate`, `balanced`, `fast`) have no GPU variant.

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

!!! Note
    A plain (non-`--cluster`) placement run just needs `REF_DIR/samples/out.fa`, which any de novo or placement backbone build already produces automatically — no extra setup. Placement **with** `--cluster` is different: it requires the backbone's samples to already be pre-split into `REF_DIR/samples/out_batch_1.fa..out_batch_N.fa` (see `BATCH_SIZE` below), which only a cluster-scale backbone build produces. See [Run ROADIES on a SLURM cluster](usage.md#run-roadies-on-a-slurm-cluster) for details.

### Grow vs. update

By default (no `--grow`), placement mode re-infers the combined species tree freely from the backbone's and query's gene trees together — the final topology can revise relationships from the original backbone tree ("**update**").

Adding `--grow` instead constrains ASTRAL-Pro3's search to the backbone tree's topology, so the backbone relationships are preserved exactly and only the new query taxa are attached ("**grow**"):

```bash
python run_roadies.py --cores 16 --mode placement --grow
```

## GPU placement

Add `--gpu N` (where `N` is the number of GPU devices to use) to run placement mode on GPU:

```bash
python run_roadies.py --cores 16 --mode placement --gpu 1
```

GPU placement requires the GPU-specific tools to be built, which `roadies_env.sh` handles automatically when the necessary build dependencies are present on your system:

- **TWILIGHT**: built automatically; uses CUDA if `nvcc` is available, falls back to a CPU build otherwise.
- **MLIPPER**: not shipped as a prebuilt binary — always built from source on your machine. `roadies_env.sh` attempts a best-effort build only if both `nvcc` (CUDA) and `libpll` (`pll.h`) are already found; it skips with a warning otherwise, since these two aren't installed automatically. Either way, MLIPPER's build also needs `gfortran`, `libblas-dev`, `liblapack-dev`, and `libtbb-dev` (not installed by `roadies_env.sh`). Once CUDA/libpll/these are available, build (or rebuild) it explicitly:
  ```bash
  bash MLIPPER/install/setup_host.sh
  ```
  Pass `--skip-apt` if you've already installed the apt dependencies yourself (the script's default apt install needs sudo); see `bash MLIPPER/install/setup_host.sh --help` for all options (e.g. pointing at a non-default `libpll` install location).

  !!! Note
      If MLIPPER fails at *runtime* with `undefined symbol: ATL_dGetNB` (not a build failure), that's an unrelated, pre-existing broken BLAS/LAPACK alternative on your system, not a MLIPPER or ROADIES bug — see [Troubleshooting: Error 7](troubleshooting.md#error-7-mlipper-fails-with-undefined-symbol-atl_dgetnb-gpu-placement-mode).
- **KegAlign**: installed via the `kegalign` conda environment (`workflow/envs/kegalign.yaml`), used automatically by Snakemake's `--use-conda` for the `kegalign` rule.

## New `config.yaml` parameters

These parameters are used by ROADIES_XP in addition to the ones documented in the [User Guide](usage.md):

| Parameter | Description | Default |
| --- | --- | --- |
| **REF_DIR** | Path to the backbone iteration's output directory, used only in `--mode placement`. Not used in de novo modes. | `null` |
| **GROUP_CSV** | Optional path to a CSV file with `species,group` columns, used to bias gene sampling towards under-represented lineages/clades rather than sampling genomes uniformly. Leave unset (or point to a nonexistent path) to sample uniformly, as in de novo mode. | `""` (unset) |
| **BATCH_SIZE** | *(`--mode placement` with `--cluster` only)* Number of loci per LASTZ batch. `REF_DIR/samples` must already contain `out_batch_1.fa..out_batch_N.fa` split at this size, i.e. `N * BATCH_SIZE` must equal `GENE_COUNT`. Ignored for plain (non-cluster) placement runs — see the note below. | 250 |
