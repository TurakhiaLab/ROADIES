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

After completing the run, the output files (along with all intermediate output files for each stage of the pipeline) will be saved in a separate `results` folder, which contains the following subfolders:

- `alignments`
- `genes`
- `msa`
- `plots`
- `samples`
- `statistics`

### Running the convergence

To rerun the workflow iteratively for checking the convergence, follow the steps below:

- Activate conda environment (make sure that ete3 is present in snakemake conda environment)
- Configure input parameters like `max_iter`, `stop_iter` in `config.yaml` file
- Provide input reference tree for the species in `trees` folder (make sure that the reference tree is in newick format)
- Run
```
python converge.py {args}
```
## Example

- `conda activate snakemake` (activates snakemake environment)
- `snakemake --core 16 --use-conda` (run snakemake on 16 cores)
- Results will be saved in `results` folder
- `python workflow/scripts/converge.py --ref trees/flies_ref_11.nwk -c 16` (run the convergence on 16 cores by taking 11 drosophila species tree as reference)

### Available reference trees

`trees` folder contain the following reference trees (along with their sources mentioned below)

- 11 Drosophila species - `flies_ref_11.nwk` ([Source]())
- 100 Drosophila species - `flies_ref_100.nwk` ([Source]())
- 


## <a name="authors"></a> Authors

## <a name="acknowledgements"></a> Acknowledgements




