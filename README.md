

<div align="center">

![ROADIES_logo](https://github.com/TurakhiaLab/wga-phylo/assets/114828525/05cd206e-542c-4ee4-bfd6-d4c03fed5984)

# Reference free Orthology free Alignment free DIscordant Aware Estimation of Species Tree (ROADIES)

</div>

## Table of Contents
- [Introduction](#overview)
- [Environment Setup](#gettingstarted) 
- [Using ROADIES](#usage)
  - [Step1: Configuring ROADIES](#configuration)
  - [Step2: Running the pipeline](#run)
  - [Step3: Analyzing output files](#output)
- [Example run](#example)
- [Contributions and Support](#support)
- [Citing ROADIES](#citation)

## <a name="overview"></a> Introduction

Welcome to the official repository of ROADIES, a novel pipeline designed for phylogenetic tree inference of the species directly from their raw genomic assemblies. Our pipeline offers a fully automated, easy-to-use, scalable solution, eliminating any error-prone manual steps and providing unique flexibility in adjusting the tradeoff between accuracy and runtime. 
<br>
#### Key Features
- **Automation**: ROADIES automates the process of species tree inference without requiring any intermediate gene annotations or orthologous groups, making it effortless for users to generate accurate species trees.
- **Scalability**: ROADIES handles both small-scale and large-scale datasets efficiently, including diverse life forms such as mammals, flies, and birds. ROADIES also scales efficiently with multiple cores and produces faster results.
- **Reference Free**: ROADIES ensures unbiased results by eliminating reference bias, enabling accurate species tree inference by randomly sampling genes from raw genome assemblies.
- **Flexibility**: ROADIES allows users to tune the tradeoff between accuracy and runtime by configuring the parameters and tailoring the pipeline to their specific needs.
- **Debugging options**: ROADIES provides multiple plots as output for graphical analysis, making it easier for the user to debug. 

#### ROADIES Pipeline Overview
ROADIES pipeline consists of multiple stages, from raw genome assemblies to species tree estimation, with several user-configurable parameters in each stage. ROADIES randomly samples subsequences from input genomic assemblies as genes which are then aligned with all individual assemblies using [LASTZ](https://lastz.github.io/lastz/). Next, ROADIES filters the alignments, gathers all homology data per gene, and performs multiple sequence alignments for every gene using [PASTA](https://github.com/smirarab/pasta). Lastly, ROADIES estimates gene trees from MSA using [IQTREE](http://www.iqtree.org/) and eventually estimates species trees from gene trees using [ASTRAL-Pro](https://github.com/chaoszhang/A-pro). 

## <a name="gettingstarted"></a> Environment Setup

This section provides detailed instructions on how to install and set up the environment to run ROADIES in your system.

ROADIES is built on Snakemake (workflow parallelization tool). It also requires various tools (PASTA, LASTZ, IQTREE, ASTRAL-Pro) to be installed before performing the analysis. To ease the process, instead of individually installing the tools, we provide a script to automatically download all dependencies into the user system. 

### Linux user

Execute bash script `roadies_env.sh` by following the commands below:

```
chmod +x roadies_env.sh
./roadies_env.sh
```
This will install and build all required tools and dependencies required by the user to get started. Once setup is complete, it will print "Setup complete" in the terminal. On its completion, a snakemake environment named "roadies_env" will be activated with all conda packages installed in it. 

## <a name="usage"></a> Using ROADIES

This section provides detailed instructions on how to configure, run, and analyze the output of ROADIES for species tree inference. Once the required environment setup process is complete, follow the steps below for using ROADIES.

### <a name="configuration"></a> Step 1: Configuring ROADIES

ROADIES provides multiple options for the user to configure the pipeline specific to their requirements before running the pipeline. Following is the list of available input configurations, provided in `config/config.yaml` (Note: ROADIES has default values for some of the parameters that give the best results, users can modify the values specific to their needs).



- **GENOMES**: Specify the path to your input files which includes raw genome assembiles of the species with `--GENOMES`. All genome assemblies should be in fasta format. The files should be named according to the species' names (for example, Aardvark's genome assembly to be named as `Aardvark.fa`)
- **REFERENCE**: Specify path for the reference tree in Newick format using `--REFERENCE` parameter. This is used if user wants to compare ROADIES' results with a state-of-the-art approach. 
- **LENGTH**: Configure the lengths of each of the sampled subsequence or genes with `--LENGTH` parameter (default is 500).
- **GENE_COUNT**: Configure the number of genes to be sampled at each iteration with `--KREG` parameter (default is 100).
- **UPPER_CASE**: ROADIES samples the genes only if the percentage of upper cases in each gene is more than a particular value, set by `--UPPER_CASE` parameter (default is 0.9). 
- **OUT_DIR**: Specify the path where you want ROADIES to store the output files.
- **MIN_ALIGN**: Specify the minimum number of allowed species to exist in gene fasta files using `--MIN_ALIGN` parameter (default is 4). This parameter is used for filtering gene fasta files which has very less species representation. It is recommended to set the value more than the default value since ASTRAL-Pro follows quartet-based topology for species tree inference.
- **LASTZ parameters**: LASTZ tool comes with several user-configurable. For ROADIES, we only configure three LASTZ parameters.
  - **COVERAGE**: Set the percentage of input sequence included in the alignment (default is 85).
  - **CONTINUITY**: Define the allowable percentage of non-gappy alignment columns (default is 85).
  - **IDENTITY**: Set the percentage of the aligned base pairs (default is 65).
- **MAX_DUP**:
- **STEPS**:
- **NUM_BOOTSTRAP**:
- **INPUT_GENE_TREES**:
- **INPUT_MAP**:
- **ITERATIONS**: Provide the maximum number of iterations for ROADIES to run with `--MAX_ITER` parameter. Set high `--MAX_ITER` if you want to run the pipeline longer to generate accurate results. Provide `--MAX_ITER` a small value if you want quicker estimate of species tree (such as guide trees for other phylogenetic tools).
- **WEIGHTED**:
- **TO_ALIGN**: ROADIES follows a weighted sampling approach to select some species out of all species every iteration where the selection of species is based on past iterations' performance. User can specify the number of species to be aligned with, using `--TO_ALIGN` parameter (default is set as half of the total number of input species). Set lower `--TO_ALIGN` if you want speedy results, set `--TO_ALIGN` with a high value if you want slower and stable results.
- **PASTA parameters**:
  - **FILTERFRAGMENTS**:
  - **MASKSITES**:
 - **MSA**:
 - **CORES**:

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

After pipeline finishes running, you can analyze the output files (along with all intermediate output files for each stage of the pipeline) saved in a separate `results` folder mentioned in `--OUT_DIR` parameter containing the following subfolders:

- `alignments` - contains the sampled output of each species in `<species_name>.maf` format generated by the sampling step
- `genes`- contains the fasta files of all highest scoring genes which is filtered after LASTZ step in `gene_<id>.fa` format
- `geneTree` - contains all gene trees merged together in a single `gene_tree_merged.nwk` file 
- `msa`- contains the MSA filtered output for all gene fasta files in `gene_aln_<id>.fa` format
- `plots` - 
- `samples` - 
- `statistics`- contains the list of input species in `num_genes.csv`, number of gene trees in `num_gt.txt`, number of homologues with corresponding genes in `homologues.csv`

## <a name="example"></a> Example run

To help you get started with ROADIES, we provide an example command below. By default, `--GENOMES` points to 48 Avian species datasets and `--MAX_ITER` is set to 50. To execute our pipeline with 16 cores and 16 jobs in parallel, run
```
python workflow/scripts/converge.py --cores 16 --jobs 16
```
This command will execute the converge script to iterate ROADIES 50 times and estimates phylogeny of 48 Avian species. It takes around XXX hours to complete the execution. The output files will be saved in separate `results` directory. The final inferred species tree in Newick format will be saved as `roadies.nwk` file in the same directory.

## <a name="support"></a> Contributions and Support

We welcome contributions from the community to enhance the capabilities of ROADIES. If you encounter any issues or have suggestions for improvement, please open an issue on GitHub. For general inquiries and support, reach out to our team.

## <a name="citation"></a> Citing ROADIES

If you use the ROADIES pipeline for species tree inference in your research or publications, we kindly request that you cite the following paper:



