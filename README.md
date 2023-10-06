

<div align="center">

![ROADIES_logo](https://github.com/TurakhiaLab/wga-phylo/assets/114828525/05cd206e-542c-4ee4-bfd6-d4c03fed5984)

# Reference free Orthology free Alignment free DIscordant Aware Estimation of Species Tree (ROADIES)

</div>

## Table of Contents
- [Introduction](#overview)
- [Environment Setup](#gettingstarted) 
- [Using ROADIES](#usage)
  - [Step1: Configuring parameters](#configuration)
  - [Step2: Running the pipeline](#run)
  - [Step3: Analyzing output files](#output)
- [Contributions and Support](#support)
- [Citing ROADIES](#citation)

## <a name="overview"></a> Introduction

Welcome to the official repository of ROADIES, a novel pipeline designed for phylogenetic tree inference of the species directly from their raw genomic assemblies. ROADIES pipeline offers a fully automated, easy-to-use, scalable solution, eliminating any error-prone manual steps and providing unique flexibility in adjusting the tradeoff between accuracy and runtime. 
<br>
#### Key Features
- **Automation**: ROADIES automates the process of species tree inference from their raw genome assemblies without requiring any intermediate gene annotations or orthologous groups, making it effortless for users to generate accurate species trees.
- **Scalability**: ROADIES handles both small-scale and large-scale datasets efficiently, including diverse life forms such as mammals, flies, and birds. ROADIES also scales efficiently with multiple cores and produces faster results.
- **Reference Free**: ROADIES ensures unbiased results by eliminating reference bias, enabling accurate species tree inference by randomly sampling genes from raw genome assemblies.
- **Flexibility**: ROADIES allows users to tune the tradeoff between accuracy and runtime by configuring the parameters and tailoring the pipeline to their specific needs.
- **Debugging options**: ROADIES provides multiple plots as output for graphical analysis, making it easier for the user to debug. 

#### ROADIES Pipeline Overview
ROADIES pipeline consists of multiple stages, from raw genome assemblies to species tree estimation, with several user-configurable parameters in each stage. ROADIES randomly samples subsequences from input genomic assemblies as genes which are then aligned with all individual assemblies using [LASTZ](https://lastz.github.io/lastz/). Next, ROADIES filters the alignments, gathers all homology data per gene, and performs multiple sequence alignments for every gene using [PASTA](https://github.com/smirarab/pasta). Lastly, ROADIES estimates gene trees from MSA using [IQTREE](http://www.iqtree.org/) and eventually estimates species trees from gene trees using [ASTRAL-Pro](https://github.com/chaoszhang/A-pro). 

<div align="center">

<img src="drawing_github.png">

</div>

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

### <a name="configuration"></a> Step 1: Configuring parameters

ROADIES provides multiple options for the user to configure the pipeline specific to their requirements before running the pipeline. Following is the list of available input configurations, provided in `config/config.yaml` (Note: ROADIES has default values for some of the parameters that give the best results, users can modify the values specific to their needs).

| Parameters | Description | Default value |
| --- | --- | --- |
| **GENOMES** | Specify the path to your input files, which includes raw genome assemblies of the species. All input genome assemblies should be in fasta format. The genome assembly files should be named according to the species' names (for example, Aardvark's genome assembly is to be named `Aardvark.fa`). Each file should contain the genome assembly of one unique species. If a file contains multiple species, split it into individual genome files (fasplit can be used for this: `faSplit byname <input_dir> <output_dir>`)| |
| **REFERENCE** (optional) | Specify the path for the reference tree (state-of-the-art) in Newick format to compare ROADIES' results with a state-of-the-art approach. If you don't want to specify any reference tree, set it to `null`. | `null` |
| **LENGTH** | Configure the lengths of each of the randomly sampled subsequences or genes. | 500 |
| **GENE_COUNT** | Configure the number of genes to be sampled across all input genome assemblies. | 750 |
| **UPPER_CASE** | Configure the lower limit threshold of upper cases for valid sampling. ROADIES samples the genes only if the percentage of upper cases in each gene is more than this value. | 0.9 (Recommended) |
| **OUT_DIR** | Specify the path for ROADIES output files. | |
| **MIN_ALIGN** | Specify the minimum number of allowed species to exist in gene fasta files after LASTZ. This parameter is used for filtering gene fasta files which has very less species representation. It is recommended to set the value more than the default value since ASTRAL-Pro follows a quartet-based topology for species tree inference. | 4 (Recommended) |
| **COVERAGE** | Set the percentage of input sequence included in the alignment for LASTZ. | 85 (Recommended) |
| **CONTINUITY** | Define the allowable percentage of non-gappy alignment columns for LASTZ. | 85 (Recommended) |
| **IDENTITY** | Set the percentage of the aligned base pairs for LASTZ. | 65 (Recommended) | 
| **MAX_DUP** | Specify max number of allowed gene copies from one input genome in an alignment. | 10|
| **STEPS** |Specify the number of steps in the LASTZ sampling (increasing number speeds up alignment but decreases LASTZ accuracy).|1 (Recommended)|
| **FILTERFRAGMENTS** | Specify the portion so that sites with less than the specified portion of non-gap characters in PASTA alignments will be masked out. If it is set to 0.5, then sites with less than 50% of non-gap characters will be masked out. | 0.5 (Recommended)|
| **MASKSITES** | Specify the portion so that sequences with less than the specified portion of non-gap sequences will be removed in PASTA alignment. If it is set to 0.05, then sequences having less than 5% of non-gap characters (i.e., more than 95% gaps) will be masked out.| 0.05 (Recommended)|
| **MSA** | Specify the modes of operation for the ROADIES pipeline. Set it to `iqtree` if you want to run in accurate mode, `fasttree` for balanced mode, and `mashtree` for fast mode of operation. | `iqtree` |
| **CORES** | Specify the number of cores in the system for ROADIES to run. ||

### <a name="run"></a> Step 2: Running the pipeline

Once the required installations are completed and the pipeline is configured, execute the following command:

```
python workflow/scripts/converge.py --cores [number of cores]
```

### <a name="output"></a> Step 3: Analyzing output files

After the pipeline finishes running, the final species tree estimated by ROADIES will be saved as `roadies.nwk` inside a separate folder mentioned in the `--OUT_DIR` parameter in the `config/config.yaml` file. 

Other intermediate output files for each stage of the pipeline are also saved in `--OUT_DIR`, containing the following subfolders:

- `alignments` - this folder contains the LASTZ alignment output of all individual input genomes aligned with randomly sampled gene sequences.
- `benchmarks` - this folder contains the runtime value of each of the individual jobs for each of the stages in the pipeline. These files will only be used if you want to estimate and compare the stagewise runtime of various pipeline stages and will not be used in final tree estimation. 
- `genes`- this folder contains the output files of multiple sequence alignment and tree-building stages (run by PASTA, IQTREE/FastTree) of the pipeline. 
- `genetrees` - this folder contains a file `gene_tree_merged.nwk`, which lists all gene trees together. This is used by ASTRAL-Pro to estimate the final species tree from the list of gene trees in the `gene_tree_merged.nwk` file.
- `plots` - this folder contains four plots
  - `gene_dup.png` - this histogram plot represents the count of the number of gene duplicates on the Y-axis vs. the number of genes having duplication on the X-axis.
  - `homologues.png` - this histogram plot represents the count of the number of genes on the Y-axis vs. the number of homologous species on the X-axis.
  - `num_genes.png` - this plot represents how many genes out of `--GENE_COUNT` parameter have been aligned to each of the input genomes after the LASTZ step. The X-axis represents different genomes, and the Y-axis represents the number of genes.
  - `sampling.png` - the plot shows how many genes have been sampled from each of the input genomes after the random sampling step. The X-axis represents different genomes, and the Y-axis represents the number of genes.
- `samples` - this folder contains the list of randomly sampled genes from individual input genomes. `<species_name>_temp.fa` files contain genes sampled from that input genome.`out.fa` combines all sampled genes from individual genomes into one file before the LASTZ step. 
- `statistics` - this folder contains CSV data for the plots shown in the `plots` directory mentioned above.
- `roadies_stats.nwk`- this is the final estimated species tree (same as `roadies.nwk`, along with the support branch values in the Newick tree). 
- `roadies.nwk`- this is the final estimated species tree in Newick format.
- `time_stamps.csv` - this file contains the start time, number of gene trees required for estimating species tree, end time, and total runtime (in seconds), respectively. 

## <a name="support"></a> Contributions and Support

We welcome contributions from the community to enhance the capabilities of ROADIES. If you encounter any issues or have suggestions for improvement, please open an issue on GitHub. For general inquiries and support, reach out to our team.

## <a name="citation"></a> Citing ROADIES

If you use the ROADIES pipeline for species tree inference in your research or publications, we kindly request that you cite the following paper:



