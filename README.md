<div align="center">
    
# Reference-free Orthology-free Annotation-free DIscordance aware Estimation of Species tree (ROADIES)

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/TurakhiaLab/ROADIES/blob/main/LICENSE

[![License][license-badge]][license-link]
[![Build Status](https://github.com/TurakhiaLab/ROADIES/actions/workflows/ci.yml/badge.svg)](https://github.com/TurakhiaLab/ROADIES/actions)
[<img src="https://img.shields.io/badge/Made with-Snakemake-green.svg?logo=snakemake">](https://snakemake.readthedocs.io/en/v7.19.1/index.html)
[<img src="https://img.shields.io/badge/Install with-Biooconda-brightgreen.svg?logo=conda">](http://bioconda.github.io/recipes/roadies/README.html)
[<img src="https://img.shields.io/badge/Install with-DockerHub-informational.svg?logo=Docker">](https://hub.docker.com/r/ang037/roadies)
[<img src="https://img.shields.io/badge/Published in-PNAS-informational.svg?logo=LOGO">](https://doi.org/10.1073/pnas.2500553122)
[<img src="https://img.shields.io/badge/DOI-10.5061/dryad.tht76hf73-yellowgreen.svg?logo=LOGO">](https://doi.org/10.5061/dryad.tht76hf73)
[<img src="https://img.shields.io/badge/Watch it on-Youtube-FF0000.svg?logo=YouTube">](https://youtu.be/1sR741TvZnM?si=vVNAnonvzNEzrLKq)

<div align="center">

<img src="docs/images/ROADIES_logo.png" style="height: 350px; width: auto;">

</div>

</div>

## Table of Contents
- [Introduction](#overview)
- [Quick Install](#usage)
- [Quick Start](#start)
- [Running ROADIES on your own data](#runpipeline)
- [ROADIES_XP: Placement mode with GPU Acceleration](#xp)
- [Citing ROADIES](#citation)

<br>

## <a name="overview"></a> Introduction

**ROADIES** is a fully automated, scalable pipeline for inferring phylogenetic species trees directly from raw genomic assemblies, eliminating manual steps of gene annotation and orthology inference.

**ROADIES_XP** (ROADIES with eXtented Placer) adds *placement mode* for growing or updating an existing species tree with new genomes, without recomputing it from scratch — with an optional GPU-accelerated variant for faster alignment/tree-building — see [ROADIES_XP](#xp) below. It's opt-in; de novo tree inference (the original ROADIES) remains the default and CPU-only.

**Requirements:** Linux (tested on Ubuntu 20.04/22.04). GPU placement (ROADIES_XP, optional) additionally needs an NVIDIA GPU with CUDA.

### 🟡 For a detailed overview of ROADIES' features and configuration options, please visit our [Wiki](https://turakhialab.github.io/ROADIES/).

### 🟡 If you encounter issues while running the pipeline, please refer to [this page](https://turakhialab.github.io/ROADIES/troubleshooting/) for common errors and troubleshooting tips.
<br>

<div align="center">

  <figure>
    <img src="docs/images/drawing_github.png" alt="ROADIES Pipeline Stages">
    <figcaption>Figure: ROADIES Pipeline Stages</figcaption>
  </figure>

</div>

<br>

## <a name="usage"></a> Quick Install

Four ways to install ROADIES — pick one. Full step-by-step instructions (including troubleshooting notes) for every option are on the [Install wiki page](https://turakhialab.github.io/ROADIES/install/).

**Option 1: Bioconda (recommended)**
```bash
conda create -n roadies_env -c bioconda -c conda-forge python=3.9 ete3 seaborn
conda activate roadies_env
conda install roadies=0.1.10
cd $CONDA_PREFIX/ROADIES   # the installed package contents live here
```
See [Install: Bioconda](https://turakhialab.github.io/ROADIES/install/#option-1-install-via-bioconda-recommended) for the remaining steps.

**Option 2: DockerHub** — pull and run the prebuilt image:
```bash
docker pull ang037/roadies:latest
docker run -it ang037/roadies:latest
```

**Option 3: Local Docker build** — clone this repo, then:
```bash
docker build -t roadies_image .
docker run -it roadies_image
```

**Option 4: Install from source** — clone this repo, install the [system dependencies](https://turakhialab.github.io/ROADIES/install/#option-4-install-via-source-script) (Java, Python 3.9+, GCC, cmake, Boost, zlib), then:
```bash
source roadies_env.sh
```
This builds/activates the `roadies_env` conda environment with everything ROADIES (and ROADIES_XP) needs.

Once installed, jump to [Quick Start](#start).

<br>

## <a name="start"></a> Quick Start

**1. Download the test dataset** (11 Drosophila genomes, run from the repo root):
```bash
mkdir -p test/test_data && cat test/input_genome_links.txt | xargs -I {} sh -c 'wget -O test/test_data/$(basename {}) {}'
```

**2. Run the pipeline.** By default ROADIES runs multiple iterations for the most accurate tree; add `--noconverge` for a quick single-iteration test run:
```bash
python run_roadies.py --cores 16              # full run (default)
python run_roadies.py --cores 16 --noconverge # quick single-iteration test
```

**3. Get the tree.** The final (unrooted) species tree is written to `OUT_DIR/roadies.nwk` (`OUT_DIR` in `config.yaml`), kept up to date after every iteration — no need to hunt through `ALL_OUT_DIR/iteration_<n>` for the latest one. Reroot it yourself, or with the provided `workflow/scripts/reroot.py` (given a rooted reference tree).

<br>

## <a name="runpipeline"></a> Running ROADIES on your own data

1. Edit `config/config.yaml`: point `GENOMES` at a directory of `.fa`/`.fa.gz` assemblies, one species per file, named after the species (e.g. `Aardvark.fa`; split multi-species files first with `faSplit byname <input_dir> <output_dir>`). Adjust any other parameters — see the [User Guide](https://turakhialab.github.io/ROADIES/usage/) for the full list.
2. Run it, optionally picking a mode (`accurate` is the default) to trade off accuracy vs. runtime:
```bash
python run_roadies.py --cores 16 --mode accurate   # or: balanced, fast
```

Per-iteration trees land in `ALL_OUT_DIR/iteration_<n>/`; the final tree is always `OUT_DIR/roadies.nwk` (`ALL_OUT_DIR`/`OUT_DIR` in `config.yaml`).

### For contributing to the code, or running on a SLURM cluster, see the [User Guide](https://turakhialab.github.io/ROADIES/usage/#run-roadies-on-a-slurm-cluster) and [Contribution guide](https://turakhialab.github.io/ROADIES/contribution/)

<br>

## <a name="xp"></a> ROADIES_XP: Placement mode with GPU Acceleration

ROADIES_XP adds **placement mode**: grow or update an existing ("backbone") species tree with new query genomes, instead of re-inferring the whole tree from scratch. It's opt-in — de novo mode (`accurate` / `balanced` / `fast`, CPU-only) remains ROADIES' default behavior, and nothing changes unless you opt into placement mode explicitly.

Placement mode runs in two variants:
- **CPU placement** (default)
- **GPU placement** (add `--gpu N`) — GPU-accelerated alignment and tree-building, needs an NVIDIA GPU with CUDA

De novo mode does not have a GPU variant.

**Run in placement mode** (add `--mode placement`, and point `GENOMES`/`REF_DIR` in `config.yaml` at your query genomes and existing backbone output directory respectively):
```bash
python run_roadies.py --cores 16 --mode placement
```

**Run placement mode on GPU:**
```bash
python run_roadies.py --cores 16 --mode placement --gpu 1
```

**Grow vs. update an existing tree in placement mode**: by default, placement mode re-infers the combined species tree freely from backbone + query gene trees ("update"). Add `--grow` to instead constrain the result to the existing backbone topology while attaching the new query taxa ("grow"):
```bash
python run_roadies.py --cores 16 --mode placement --grow
```

Building the GPU/placement dependencies is handled automatically by `roadies_env.sh` when the required build dependencies (CUDA, libpll) are available.

### For full details on placement mode, GPU requirements, and new `config.yaml` parameters (`REF_DIR`, `GROUP_CSV`, `BATCH_SIZE`), refer to the [ROADIES_XP Wiki page](https://turakhialab.github.io/ROADIES/roadies_xp/)

<br>

## <a name="citation"></a> Citing ROADIES

If you use ROADIES in your research or publications, please cite the following paper:

A. Gupta, S. Mirarab, & Y. Turakhia, Accurate, scalable, and fully automated inference of species trees from raw genome assemblies using ROADIES, Proc. Natl. Acad. Sci. U.S.A. 122 (19) e2500553122, [https://doi.org/10.1073/pnas.2500553122](https://doi.org/10.1073/pnas.2500553122) (2025).

A manuscript describing ROADIES_XP (placement mode with GPU acceleration) is in preparation as a separate publication. Citation details will be added here once available.

### Accessing ROADIES output files

The output files with the gene trees and species trees generated by ROADIES in the manuscript are deposited to [Dryad](https://datadryad.org/stash). To access it, please refer to the following:

Gupta, Anshu; Mirarab, Siavash; Turakhia, Yatish (2024). Accurate, scalable, and fully automated inference of species trees from raw genome assemblies using ROADIES [Dataset]. Dryad. [https://doi.org/10.5061/dryad.tht76hf73](https://doi.org/10.5061/dryad.tht76hf73).


