# Estimating Ultra-Large Phylogenies from Raw Genome Sequences (wga-phylo)

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Authors](#authors)
- [Acknowledgements](#acknowledgements)

## <a name="overview"></a> Overview

Most phylogeny estimation methods use gene markers/gene trees which require computationally expensive, error-prone, and/or semi-automated steps to infer orthology.Our technique uses raw genomes to build phylogenies. <br>

WGA Pipeline steps:
- Take random sequences of fixed sized sampled from different genomes and treat them as genes
- Cluster species through alignments with those genes
- Build gene trees using these clusters
- Reconstruct a species tree using these gene trees using ASTRAL-PRO
- ASTRAL-pro can work with homology only, and does not require orthology

## <a name="installation"></a> Installation

### Required Installations

- Snakemake
- ASTRAL-Pro

### Install Snakemake

[Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and Snakedeploy are best installed via the Mamba package manager (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). 

Given that Mamba is installed, run 

```
mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
``` 

to install both Snakemake and Snakedeploy in an isolated environment. 

Once snakemake is installed, it needs to be activated via

```
conda activate snakemake
```
### Install ASTRAL-Pro

In order to run the pipeline, a manual installation of ASTRAL-PRO (latest version) is required. If previous version of ASTER-Linux is installed, replace the ASTER-Linux directory with the ASTRAL-PRO git repository as mentioned below. 

```
rmdir ASTER-Linux
wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip
unzip Linux.zip
cd ASTER-Linux
make
cd ..
```

## <a name="usage"></a> Usage

### Input Requirements

All input genomic sequences should be in fasta format. The path for input dataset, along with the input configuration parameters should be provided in `config/config.yaml` file.

### Running Snakemake

Once snakemake environment is activated and `config/config.yaml` is configured, run
```
snakemake --core [number of cores] --use-conda
```
For starting the run from an incomplete previous session instead of starting a fresh run everytime, run
```
snakemake --core [number of cores] --use-conda --rerun-incomplete
```
### Output options


### Requirements for running the convergence:\

- Working Snakemake installation
- ete3 in Snakemake conda environment

```
conda activate snakemake
python converge.py {args}
```

#### List of arguments for convergence\
--ref {reference tree to compare to}\
--input_gt {input gene trees newick for A-Pro}\
--input_map {input mapping file for A-Pro}\
-c {cpu cores to use}\
-l {gene length}\
-k {number of genes sampled*}\
-t {TreeDistance threshold for stopping}\
--bootstrap {number of bootstrap trees to create for comparisons}\
--max_iter {maximum number of converge runs}\
--stop_ter {number of consecutive bootstrapped self_dists and ref_dists satisfying threshold before stopping}\
--out_dir {converge output directory}\
--roadies_dir {roadies output directory}

Running drosophila dataset\
Default Settings are set to avian dataset to change genomes:
edit config/config.yaml GENOMES to "/home/roadies-datasets/drosophila"

Suggested command:\
`python converge.py -c 16 --
ref trees/refTree.nwk --stop_iter 1 `

Reference Trees:\
Avian: trees/cn48.nwk\
Drosphila: trees/refTree.nwk

*after filtering number of genes might be less than k per iteration

## <a name="authors"></a> Authors

## <a name="acknowledgements"></a> Acknowledgements




