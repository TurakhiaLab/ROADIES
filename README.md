# Estimating Ultra-Large Phylogenies from Raw Genome Sequences (wga-phylo)

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Authors](#authors)
- [Acknowledgements](#acknowledgements)

## <a name="overview"></a> Overview

## <a name="installation"></a> Installation

Snakemake and Snakedeploy are best installed via the Mamba package manager (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via Mambaforge. 

Given that Mamba is installed, run

mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
to install both Snakemake and Snakedeploy in an isolated environment. For all following commands ensure that this environment is activated via
conda activate snakemake
SnakeMake workflow for sequence selection to lastz, change your configuration data in configuration file

`snakemake --core [number of cores] --use-conda`

Snakemake will automatically detect the main Snakefile in the workflow subfolder and execute the workflow module that has been defined by the deployment in step 2.

In order to run the pipeline, a manual installation of ASTRAL-PRO(latest version) is required. Replace the ASTER-Linux directory with the ASTRAL-PRO git repository. 

`rmdir ASTER-Linux`\
`wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip `\
`unzip Linux.zip`\
`cd ASTER-Linux`\
`make`\
`cd ..`

Running converge:\
REQUIREMENTS:

Working Snakemake installation\
ete3 in Snakemake conda environment

`conda activate snakemake`\
`python converge.py {args}`

## <a name="usage"></a> Usage

list of arguments\
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




