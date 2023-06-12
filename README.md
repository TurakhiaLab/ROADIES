![ROADIES_logo](https://github.com/TurakhiaLab/wga-phylo/assets/114828525/05cd206e-542c-4ee4-bfd6-d4c03fed5984)
# Reference free Orthology free Alignment free DIscordant Aware Estimation of Species Tree (ROADIES)

## Table of Contents
- [Overview](#overview)
- [Installation](#gettingstarted) 
- [Using ROADIES](#usage)
  - [Step1: Configuring ROADIES](#configuration)
  - [Step2: Running the pipeline](#run)
  - [Step3: Analyzing output files](#output)
- [Example run](#example)
- [References](#references)

## <a name="overview"></a> Overview

Welcome to the official repository of ROADIES, a novel pipeline designed for phylogenetic tree inference of the species directly from their raw genomic assemblies. Our pipeline offers a fully-automated, easy-to-use, scalable solution, eliminating any reference bias and provides unique flexility in adjusting the tradeoff between accuracy and runtime. 
<br>
#### Key Features
- **Automation**: ROADIES automates the process of species tree inference without requiring any intermediate gene annotations or orthologous groups, making it effortless for users to generate accurate species tree.
- **Scalability**: ROADIES handles both small-scale or large-scale datasets efficiently, including diverse life forms such as mammals, flies, and birds. ROADIES also scales efficiently with multiple cores and produce faster results.
- **Reference Free**: ROADIES ensures unbiased results by eliminating reference bias, enabling accurate species tree inference from raw genome assemblies.
- **Flexibility**: ROADIES provides users to tune the tradeoff between accuracy and runtime by configuring the parameters, tailoring the pipeline to their specific needs.
- **Debugging options**: ROADIES provides multiple plots as output for graphical analysis, making easier for user to debug. 

#### ROADIES Pipeline
ROADIES pipeline consists of multiple stages from raw genome assemblies to species tree estimation, with several user configurable parameters in each stages. ROADIES samples subsequences from input genomic assemblies as genes which is then pairwise aligned with all assemblies using LASTZ. Next, ROADIES filter the alignments and perform multiple sequence alignment using PASTA for individual genes across all species. Lastly, ROADIES estimates gene trees from MSA using IQTREE and eventually estimates species tree from gene trees using ASTRAL-Pro.

## <a name="gettingstarted"></a> Installation

This section provides detailed instruction on how to install and set up ROADIES to get started.

### Dependencies
- Snakemake
- ASTRAL-Pro

### Install Snakemake

To install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), it is recommended to install via the Mamba package manager (a drop-in replacement for conda)

If you have neither Conda nor Mamba, install it by following the steps [here](https://github.com/conda-forge/miniforge#mambaforge). 

Given that Mamba is installed, run
```
mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
``` 
to install both Snakemake and Snakedeploy in an isolated environment. 

Once snakemake is installed, activate the environment via
```
conda activate snakemake
```
### Install ASTRAL-Pro

In order to run ROADIES, a manual installation of ASTRAL-Pro (latest version) is required (if previous version of ASTER-Linux is installed, replace the ASTER-Linux directory with the ASTRAL-Pro git repository as mentioned below)

```
rmdir ASTER-Linux
wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip
unzip Linux.zip
cd ASTER-Linux
make
cd ..
```

## <a name="usage"></a> Using ROADIES

This section provides the detailed instruction on how to configure, run and analyze the output of ROADIES for species tree inference. Once the required installation process is complete, follow the steps below for using ROADIES.

### <a name="configuration"></a> Step 1: Configure ROADIES

ROADIES provides multiple option for the user to configure the pipeline specific to their requirements before running the pipeline. Following is the list of available input configurations, provided in `config/config.yaml` (Note: ROADIES has default values for some of the parameters, users are required to modify the values specific to their needs).

- **Input Files**: 
- **Gene Length**: 
- **Number of genes**: 
- **Output directory**: 
- **Maximum iterations**: 
- **LASTZ parameters**: 
- **Mininum number of species in gene fasta**:
- **Reference tree**:
- **Convergence parameters**:
- **Weighted Sampling**:

### <a name="run"></a> Step 2: Running the pipeline

Once the required installations are completed and the pipeline is configured, follow the steps below:

**Step 1**: Activate snakemake
```
conda activate snakemake
```
**Step 2**: Run the pipeline by executing the following command
```
python workflow/scripts/converge.py --cores [number of cores] --jobs [number of jobs]
```
**Note**: For starting the run from an incomplete previous session instead of starting a fresh run everytime, run
```
snakemake --core [number of cores] --use-conda --rerun-incomplete
```
### <a name="output"></a> Step 3: Analyzing output files

After completing the run, the output files (along with all intermediate output files for each stage of the pipeline) will be saved in a separate `results` folder mentioned in `OUT_DIR` parameter, which contains the following subfolders:

- `alignments` - contains the sampled output of each species in `<species_name>.maf` format generated by the sampling step
- `genes`- contains the fasta files of all highest scoring genes which is filtered after LASTZ step in `gene_<id>.fa` format
- `geneTree` - contains all gene trees merged together in a single `gene_tree_merged.nwk` file 
- `msa`- contains the MSA filtered output for all gene fasta files in `gene_aln_<id>.fa` format
- `plots` - 
- `samples` - 
- `statistics`- contains the list of input species in `num_genes.csv`, number of gene trees in `num_gt.txt`, number of homologues with corresponding genes in `homologues.csv`

`trees` folder contain the following reference trees (along with their sources mentioned below)

- 11 Drosophila species - `flies_ref_11.nwk` ([Source](http://timetree.org/))
- 100 Drosophila species - `flies_ref_100.nwk` ([Source](https://github.com/flyseq/drosophila_assembly_pipelines/blob/master/figure_data/figure5/busco_species_astral.tree))

## <a name="example"></a> Example run
```
python workflow/scripts/converge.py --ref trees/flies_ref_11.nwk -c 16
```
The above command runs the convergence script on 16 cores by taking 11 drosophila species tree as reference. 

The output after the converge run is saved in a separate `converge` folder, which has the following files:

## <a name="references"></a> References





