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

ROADIES provides multiple option for the user to configure the pipeline specific to their requirements before running the pipeline. Following is the list of available input configurations, provided in `config/config.yaml` (Note: ROADIES has default values for some of the parameters which gives the best results, users can modify the values specific to their needs).

- **Input Files**: Specify the path to your input files which includes raw genome assembiles of the species with `--INPUT_DATASET`. All genome assemblies should be in fasta format. The files should be named according to the species' names (for example, Aardvark's genome assembly to be named as `Aardvark.fa`)
- **Gene Length**: Configure the lengths of each of the sampled subsequence or genes with `--LENGTH` parameter (default is 500).
- **Number of genes**: Configure the number of genes to be sampled at each iteration with `--KREG` parameter (default is 100).
- **Output directory**: Specify the path where you want ROADIES to store the output files
- **Maximum iterations**: Provide the maximum number of iterations for ROADIES to run with `--MAX_ITER` parameter. Set high `--MAX_ITER` if you want to run the pipeline longer to generate accurate results. Provide `--MAX_ITER` a small value if you want quicker estimate of species tree (such as guide trees for other phylogenetic tools). 
- **LASTZ parameters**: LASTZ tool comes with several user-configurable. For ROADIES, we only configure three LASTZ parameters. 1. `--COVERAGE` which sets the percentage of input sequence included in the alignment (default is 85), 2. `--CONTINUITY` which defines the allowable percentage of non-gappy alignment columns (default is 85), 3. `--IDENTITY` which sets the percentage of the aligned base pairs (default is 65). Modify these parameters based on your use-cases
- **Mininum number of species in gene fasta**: Specify the minimum number of allowed species to exist in gene fasta files using `--MIN_ALIGN` parameter (default is 4). This parameter is used for filtering gene fasta files which has very less species representation. It is recommended to set the value more than the default value since ASTRAL-Pro follows quartet-based topology for species tree inference. 
- **Reference tree**: Specify path for the reference tree in Newick format using `--REFERENCE` parameter. This is used if user wants to compare ROADIES' results with a state-of-the-art approach. 
- **Convergence parameters**: ROADIES provides few parameters to configure the convergence criteria. The default convegrence criteria checks the absolute difference in the mean of two windows of bootstrapped distances. This value is compared with a threshold to converge and terminate the pipeline. Set the number of times the final species needs to be bootstrapped with `--BOOSTRAP` parameter (default is 10). Set the window size with `--STOP_ITER` parameter (default is 5). Threshold can be configured using `--DIST_THRESHOLD` parameter (default is 0.01). It is recommended to follow the default convergence criteria for best results, however, users can modify it according to their datasets.
- **Weighted Sampling**: ROADIES follows a weighted sampling approach to select some species out of all species every iteration where the selection of species is based on past iterations' performance. User can specify the number of species to be aligned with, using `--TO_ALIGN` parameter (default is set as half of the total number of input species). Set lower `TO_ALIGN` if you want speedy results, set `--TO_ALIGN` with a high value if you want slower and stable results. 

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

## <a name="example"></a> Example run
```
python workflow/scripts/converge.py --ref trees/flies_ref_11.nwk -c 16
```
The above command runs the convergence script on 16 cores by taking 11 drosophila species tree as reference. 

The output after the converge run is saved in a separate `converge` folder, which has the following files:

## <a name="references"></a> References





