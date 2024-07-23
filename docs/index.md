 <div align="center">

<img src="images/ROADIES_logo.png" width="300" height="300"/>

</div>

# Reference-free Orthology-free Alignment-free DIscordance aware Estimation of Species tree (ROADIES)

## ROADIES Video Tutorial

<iframe width="1000" height="600" src="https://youtu.be/1sR741TvZnM" frameborder="0" allowfullscreen></iframe>

## Introduction

Welcome to the official wiki of ROADIES, a novel pipeline designed for phylogenetic tree inference of the species directly from their raw genomic assemblies. ROADIES offers a fully automated, easy-to-use, scalable solution, eliminating any error-prone manual steps and providing unique flexibility in adjusting the tradeoff between accuracy and runtime. 
<br>

## Key Features
- **Reference-free**: ROADIES ensures unbiased results by eliminating reference bias, enabling accurate species tree inference by randomly sampling genes from raw genome assemblies.
- **Orthology-free**: ROADIES automates the process of species tree inference from their raw genome assemblies without requiring any intermediate gene annotations or orthologous groups. It allows multi-copy gene trees (inferred from homologous regions) and does not require the challenging and error-prone step of orthologous detection prior to gene tree inference. 
- **Alignment-free**: ROADIES does not require any input alignment to be provided, as all alignments required for gene tree inference are constructed within the ROADIES pipeline itself.
- **Discordance-aware**: Instead of single-copy genes, ROADIES considers multi-copy genes while analyzing species tree and takes care of the possible gene discordances such as paralogs, horizonal gene transfer, incomplete lineage sorting. It uses a state-of-the-art and statistically consistent discordance-aware method to combine gene trees into a species tree.
- **Scalability**: ROADIES handles both small-scale and large-scale datasets efficiently, including diverse life forms such as mammals, flies, and birds. ROADIES also scales efficiently with multiple cores and produces faster results.
- **Flexibility**: ROADIES allows users to tune the tradeoff between accuracy and runtime by configuring the parameters and tailoring the pipeline to their specific needs.
- **Debugging options**: ROADIES provides multiple plots as output for graphical analysis, making it easier for the user to debug. 

## ROADIES Pipeline Overview
ROADIES pipeline consists of multiple stages, from raw genome assemblies to species tree estimation, with several user-configurable parameters in each stage. 

- **Stage 1: Random Sampling**: ROADIES randomly samples a configured number of subsequences from input genomic assemblies. Each of the subsequences is treated as a gene.
- **Stage 2: Pairwise alignment**: All sampled subsequences or genes are aligned with all input assemblies individually to find the homologous regions using the pairwise alignment tool [LASTZ](https://lastz.github.io/lastz/). 
- **Stage 3: Filtering of alignments**: ROADIES filters the low-quality alignments to reduce further redundant computation and limits the number of repetitive alignments from small genomic regions. 
- **Stage 4: Multiple sequence alignment**: ROADIES gathers the homologous regions for all genes across species and performs multiple sequence alignments for each of them using [PASTA](https://github.com/smirarab/pasta). 
- **Stage 5: Gene tree estimation**: ROADIES estimates gene trees from multiple sequence alignments of each of the genes using the maximum-likelihood based tree estimation tool [RAxML-NG](https://github.com/amkozlov/raxml-ng).
- **Stage 6: Species tree estimation**: After getting all the gene trees, ROADIES performs species tree estimation using a widely adopted discordance aware species tree estimation tool [ASTRAL-Pro](https://github.com/chaoszhang/A-pro). 

<div align="center">

<img src="images/drawing_github.png"width="1000" height="300" />

</div>

## Modes of operation

ROADIES supports multiple modes of operation based on various user requirements considering the tradeoff between accuracy and runtime. 

- **Accurate-Mode**: This is the default mode of operation and is preferred for accuracy-critical use cases. Here, the multiple sequence alignment stage is performed by [PASTA](https://github.com/smirarab/pasta) and the tree building stage is governed by [RAxML-NG](https://github.com/amkozlov/raxml-ng).
- **Fast-Mode**: This mode of operation is preferred for achieving faster results, for runtime-critical use cases. Here, the multiple sequence alignment and tree building stage is performed by [MashTree](https://github.com/lskatz/mashtree).
- **Balanced-Mode**: This mode of operation is preferred where the user wants an optimal runtime vs accuracy tradeoff. Here, the multiple sequence alignment stage is performed by [PASTA](https://github.com/smirarab/pasta), and the tree building stage is performed using [FastTree](http://www.microbesonline.org/fasttree/). 

!!! Note
    These modes of operation can be modified using command line argument `--mode` (details mentioned in the [Usage](index.md#other-command-line-arguments) section).

## Convergence Mechanism

The initial count of the genes is crucial to get the accurate species tree at the end. The number of genes sufficient for getting the accurate tree also varies with datasets. Hence, ROADIES incorporates an adaptive algorithm for establishing accurate trees by tracking its confident scores. It performs multiple iterations of the entire pipeline and stops if it gets the confident tree, otherwise it continues with more gene counts. The confidence of the tree is evaluated by the confidence of its branches (or local posterior probability). The tree having most of the confident branches with high posterior probability are considered to be confident and stable. 

!!! Note
    [ASTRAL-Pro2](https://github.com/chaoszhang/A-pro) provides the information of all the internal nodes in the form of quartets (and its support values such as local posterior probability) for every species tree per iteration. ROADIES gathers this information and keeps track of all the nodes with high support values. If the percentage change in the number of highly supported nodes gets minimal with a given number of iterations, then we say that the species tree is now converged.

!!! Note
    Users have the option to run ROADIES in converge mode using `--converge` argument (details mentioned in [Usage](index.md#other-command-line-arguments) section).

<img src="images/converge_manuscript.png"width="1000" height="300" />

## Ways to install ROADIES

### Using DockerHub

To run ROADIES using DockerHub, follow these steps:

1. Pull the ROADIES Docker image from DockerHub:

```
docker pull ang037/roadies:latest
```
2. Run the Docker container:

```
docker run -it ang037/roadies:latest
```

### Using Docker locally

First, clone the repository (requires `git` to be installed in the system):

```
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
```

Then build and run the Docker container:

```
docker build -t roadies_image .
docker run -it roadies_image
```

### Using installation script (requires sudo access)

First clone the repository:

```
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
```

Then, execute the installation script:

```
chmod +x roadies_env.sh
source roadies_env.sh
```

This will install and build all tools and dependencies. Once the setup is complete, it will print `Setup complete` in the terminal and activate the `roadies_env` environment with all Conda packages installed. 

!!! Note 
    ROADIES is built on [Snakemake (workflow parallelization tool)](https://snakemake.readthedocs.io/en/stable/). It also requires various tools (PASTA, LASTZ, RAxML-NG, MashTree, FastTree, ASTRAL-Pro2) to be installed before performing the analysis. To ease the process, instead of individually installing the tools, we provide `roadies_env.sh` script to automatically download all dependencies into the user system.

#### Required dependencies

To run this script, ensure the following dependencies are installed:
- Java Runtime Environment (version 1.7 or higher)
- Python (version 3 or higher)
- `wget` and `unzip` commands
- GCC (version 11.4 or higher)
- cmake (Download here: https://cmake.org/download/)
- Boost library (Download here: https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/)
- zlib (Download here: http://www.zlib.net/)

For Ubuntu, you can install these dependencies with: 

```
sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git default-jre libgomp1 libboost-all-dev cmake
```

!!! Warning
    If you encounter issues with the Boost library, add its path to `$CPLUS_LIBRARY_PATH` and save it in `~/.bashrc`.


## Quick start (with provided test dataset)

Once setup is done, you can run the ROADIES pipeline using the provided test dataset. Follow these steps for a 16-core machine:

1. Go to ROADIES repository directory if not there:

```
cd ROADIES
```

2. Create a directory for the test data and download the test datasets (using the following one line command):

```
mkdir -p test/test_data && cat test/input_genome_links.txt | xargs -I {} sh -c 'wget -O test/test_data/$(basename {}) {}'
```
3. Run the pipeline with the following command (from ROADIES directory):

```
python run_roadies.py --cores 16
```

The second command will download the 11 Drosophila genomic datasets (links provided in `test/input_genome_links.txt`) and save them in the `test/test_data` directory. The third command will run ROADIES for those 11 Drosophila genomes and save the final newick tree as `roadies.nwk` in a separate `output_files` folder upon completion.

**Running ROADIES with different modes of operation**: To run ROADIES in various other modes of operation (fast, balanced, accurate) (description of these modes are mentioned in [Modes of operation](index.md#modes-of-operation) section), try the following commands:

```
python run_roadies.py --cores 16 --mode accurate
```

```
python run_roadies.py --cores 16 --mode balanced
```

```
python run_roadies.py --cores 16 --mode fast
```
!!! Note
    Accurate mode is the default mode of operation. If you don't specify any particular mode using `--mode` argument, default mode will run.

For each modes, the output species tree will be saved as `roadies.nwk` in a separate `output_files` folder.

**Running ROADIES in converge mode**: To run ROADIES with converge mode (details mentioned in [convergence mechanism](index.md#convergence-mechanism) section), run the following command (notice the addition of `--converge` argument):

```
python run_roadies.py --cores 16 --converge
```

Try following commands for other modes:

```
python run_roadies.py --cores 16 --mode balanced --converge
```
```
python run_roadies.py --cores 16 --mode fast --converge
```

The output files for all iterations will be saved in a separate `converge_files` folder. `output_files` will save the results of the last iteration. Species tree for all iterations will be saved in `converge_files` folder with the nomenclature `iteration_<iteration_number>.nwk`.

## Detailed Usage

This section provides detailed instructions on how to configure the ROADIES pipeline further for various user requirements with your own genomic dataset. Once the required environment setup process is complete, follow the steps below.

### Step 1: Specify input genomic dataset

After installing the environment, you need to get input genomic sequences for creating the species tree. To run ROADIES with your own dataset, update the `config.yaml` file (found in the ROADIES directory - `config` folder) to include the path to your input datasets under the `GENOMES` parameter.

!!! Note 
    All input genome assemblies in the path mentioned in `GENOMES` should be in `.fa` or `.fa.gz` format. The genome assembly files should be named according to the species' names (for example, Aardvark's genome assembly is to be named `Aardvark.fa`). Each file should contain the genome assembly of one unique species. If a file contains multiple species, split it into individual genome files (fasplit can be used for this: `faSplit byname <input_dir> <output_dir>`). Moreover, the file name should not have any special characters like `.` (apart from `_`) - for example, if the file name is `Aardvark.1.fa`, rename it to `Aardvark_1.fa`.

### Step 2: Modify Other Configuration Paramters

Adjust other parameters listed in `config.yaml` as per specific user requirements. Details of the parameters are mentioned below.
 
!!! Note
    ROADIES has default values for some of the parameters that give the best results and are recommended in general. However, users can optionally modify the values specific to their needs.

| Parameters | Description | Default value |
| --- | --- | --- |
| **GENOMES** | Specify the path to your input files which includes raw genome assemblies of the species. | |
| **REFERENCE** (optional) | Specify the path for the reference tree (state-of-the-art) in Newick format to compare ROADIES' results with a state-of-the-art approach. If you don't want to specify any reference tree, set it to `NULL`. | `NULL` |
| **LENGTH** | Configure the lengths of each of the randomly sampled subsequences or genes. | 500 |
| **GENE_COUNT** | Configure the number of genes to be sampled across all input genome assemblies. In normal mode, this will be the count of the genes to be sampled. In `--converge` mode, this will be the initial count of the number of genes for the first iteration and this value will be doubled iteratively. | 250 |
| **UPPER_CASE** | Configure the lower limit threshold of upper cases for valid sampling. ROADIES samples the genes only if the percentage of upper cases in each gene is more than this value. | 0.9 (Recommended) |
| **OUT_DIR** | Specify the path for ROADIES output files (this saves the current iteration results in converge mode). | |
| **ALL_OUT_DIR** | Specify the path for ROADIES output files for all iterations in converge mode. | |
| **MIN_ALIGN** | Specify the minimum number of allowed species to exist in gene fasta files after LASTZ. This parameter is used for filtering gene fasta files which has very less species representation. It is recommended to set the value greater than or equal to 4 since ASTRAL-Pro follows a quartet-based topology for species tree inference. For larger evolutionary timescales, we recommended setting it to a much higher value. In such cases, 15 to 20 would be a good start. | 4 |
| **COVERAGE** | Set the percentage of input sequence included in the alignment for LASTZ. | 85 |
| **CONTINUITY** | Define the allowable percentage of non-gappy alignment columns for LASTZ. | 85 |
| **IDENTITY** | Set the percentage of the aligned base pairs (matches/mismatches) for LASTZ. For larger evolutionary timescales, consider lowering the identity values than default for more homologous hits to be encountered. | 65 | 
| **MAX_DUP** | Specify maximum number of allowed gene copies from one input genome in an alignment. | 10|
| **STEPS** |Specify the number of steps in the LASTZ sampling (increasing number speeds up alignment but decreases LASTZ accuracy).|1 |
| **FILTERFRAGMENTS** | Specify the portion so that sites with less than the specified portion of non-gap characters in PASTA alignments will be masked out. If it is set to 0.5, then sites with less than 50% of non-gap characters will be masked out. | 0.5 |
| **MASKSITES** | Specify the portion so that sequences with less than the specified portion of non-gap sequences will be removed in PASTA alignment. If it is set to 0.05, then sequences having less than 5% of non-gap characters (i.e., more than 95% gaps) will be masked out.| 0.02 |
| **SUPPORT_THRESHOLD** | Specify the threshold so that support values with equal to or higher than this threshold is considered as highly supported node. Such highly supported nodes crossing this threshold will be counted at every iteration to check the confidence of the tree (works in `--converge` mode). | 0.95 |
| **NUM_INSTANCES** | Specify the number of instances for PASTA, LASTZ, MashTree and RAxML-NG to run in parallel. It is recommended to set the number of instances equal to (`--cores`/4) for optimal runtime. | 4 | 

### Step 3: Run the ROADIES pipeline

Once the required installations are completed and the parameters are configured in `config.yaml` file, execute the following command (from ROADIES repo home directory):

```
python run_roadies.py --cores <number of cores>
```

This will let ROADIES run in accurate mode by default with specified number of cores. After the completion of the execution, the output species tree in Newick format will be saved as `roadies.nwk` in a separate `output_files` folder.

### Command line arguments

There are multiple command line arguments through which user can change the mode of operation, specify the custom config file path, etc.

| Argument | Description |
| --- | --- |
| `--cores` | Specify the number of cores |
| `--mode` | Specify [modes of operation](index.md#modes-of-operation) (`accurate`, `balanced` or `fast`).`accurate` mode is the default mode. | 
| `--converge` | Run ROADIES in [converge](index.md#convergence-mechanism) mode if you do not know the optimal gene count to start with |
| `--config` | Provide optional custom YAML files (in the same format as `config.yaml` provided with this repository). If not given, by default `config/config.yaml` file will be considered.|

For example:

```
python run_roadies.py --cores 16 --mode balanced --converge --config config/config.yaml
```

Use `--help` to get the list of command line arguments.

### Step 4: Analyze output files

#### Without convergence

After the pipeline finishes running, the final species tree estimated by ROADIES will be saved as `roadies.nwk` inside a separate folder mentioned in the `--OUT_DIR` parameter in the `config/config.yaml` file. 

ROADIES also provides a number of intermediate output files for extensive debugging by the user. These files will be saved in `--OUT_DIR`, containing the following subfolders:

1. `alignments` - this folder contains the LASTZ alignment output of all individual input genomes aligned with randomly sampled gene sequences.
2. `benchmarks` - this folder contains the runtime value of each of the individual jobs for each of the stages in the pipeline. These files will only be used if you want to estimate and compare the stagewise runtime of various pipeline stages and will not be used in final tree estimation. 
3. `genes` - this folder contains the output files of multiple sequence alignment and tree-building stages (run by PASTA, IQTREE/FastTree, MashTree) of the pipeline. 
4. `genetrees` - this folder contains two files as follows:
    - `gene_tree_merged.nwk` - this file lists all gene trees together generated by IQTREE/FastTree/MashTree. It is used by ASTRAL-Pro to estimate the final species tree from this list of gene trees.
    - `original_list.txt` - this file lists all gene trees together corresponding to their gene IDs. Some lines will have only gene IDs but no associated gene trees. This is because some genes will be filtered out from tree building and MSA step if it has less than four species. Hence this file also lists those gene IDs with missing gene trees for further debugging. 
5. `plots` - this folder contains four following plots:
    - `gene_dup.png` - this histogram plot represents the count of the number of gene duplicates on the Y-axis vs. the number of genes having duplication on the X-axis.
    - `homologues.png` - this histogram plot represents the count of the number of genes on the Y-axis vs. the number of homologous species on the X-axis.
    - `num_genes.png` - this plot represents how many genes out of `--GENE_COUNT` parameter have been aligned to each of the input genomes after the LASTZ step. The X-axis represents different genomes, and the Y-axis represents the number of genes.
    - `sampling.png` - the plot shows how many genes have been sampled from each of the input genomes after the random sampling step. The X-axis represents different genomes, and the Y-axis represents the number of genes.
6. `samples` - this folder contains the list of randomly sampled genes from individual input genomes. 
    - `<species_name>_temp.fa` - these files contain genes sampled from the particular input genome.
    - `out.fa` - this file contains all sampled subsequences (genes) from individual genomes combined, which is given to the the LASTZ step. 
7. `statistics` - this folder contains CSV data for the plots shown in the `plots` directory mentioned above.
    - `gene_to_species.csv` - this is an additional CSV file (corresponding plots to be added in future) which provides the information about which genes are aligned to what species after LASTZ step (`num_genes.csv` only gives the total count of the genes per species, `gene_to_species.csv` also gives the ID number of those aligned genes). Along with each gene ID number, it also provides the [score, line number in .maf file, position] of all the homologs of that particular gene. Score, position and line number information is collected from the corresponding species' .maf file (generated by LASTZ), saved in `results/alignments` folder.
8. `roadies_stats.nwk`- this is the final estimated species tree (same as `roadies.nwk`), along with the support branch values in the Newick tree. 
9. `roadies.nwk`- this is the final estimated species tree in Newick format.
10. `roadies_rerooted.nwk` (optional) - this is the final estimated species tree, re-rooted corresponding to the outgroup node from the given reference tree (provided as `REFERENCE` in `config.yaml`).
11. `time_stamps.csv` - this file contains the start time, number of gene trees required for estimating species tree, end time, and total runtime (in seconds), respectively.
12. `ref_dist.csv` - this file provides the number of gene trees and the Normalized Robinson-Foulds distance between the final estimated species tree (i.e., `roadies.nwk`) and the reference tree (i.e., REFERENCE parameter in `config.yaml`).

#### With convergence

If converge option is enabled, the results of all iterations (along with the corresponding species tree in the name `iteration_<iteration_number>.nwk`) will be saved in a separate folder mentioned in the `--ALL_OUT_DIR` parameter in the `config/config.yaml` file.

!!! Note
    With `--converge` option, `--OUT_DIR` saves the results of the current ongoing iteration (if pipeline execution is finished, then the last iteration), whereas `--ALL_OUT_DIR` saves the results of all iterations executed. 

For extensive debugging, other intermediate output files for each stage of the pipeline for each iterations are saved in `--ALL_OUT_DIR` as follows:

1. Folder with `iteration_<iteration_number>` - this folder contains results from the specific iteration corresponding to the iteration number in the folder name.
    -  Folder with name in `--OUT_DIR` - this contains the results of all stages of the pipeline (as described above in non convergence section). 
    - `gene_tree_merged.nwk` - this file lists all gene trees together generated by IQTREE/FastTree/MashTree in that particular iteration. It is concatenated with master list of gene trees from all past iterations before providing to ASTRAL-Pro to estimate the final converged species tree.
    - `iteration_<iteration_number>.log` - this file contains the log information of the corresponding iteration execution. 
    - `mapping.txt` - This file maps all gene names in the gene trees with the corresponding species name from where it originates. It is required by ASTRAL-Pro, along with the master list of gene trees from all iterations, to infer species tree. 
2. `iteration_<iteration_number>_stats.nwk` - this is the final estimated species tree for the corresponding iteration (same as `iteration_<iteration_number>.nwk`), along with the support branch values in the Newick tree. 
3. `iteration_<iteration_number>.nwk` - this is the final estimated species tree for the corresponding iteration
4. `iteration_<iteration_number>.rerooted.nwk` - (optional) - this is the final estimated species tree for the corresponding iteration, re-rooted to the outgroup node from the given reference tree (provided as `REFERENCE` in `config.yaml`).
5. `master_gt.nwk` - this is the concatenated list of all gene trees from all iterations together.
6. `master_map.txt` - this is the concatenated list of all mapping files from all iterations together. This `master_gt.nwk` and `master_map.txt` is provided to ASTRAL-Pro after every iteration to get the converged species tree. 
7. `ref_dist.csv` - this file provides the iteration number, number of gene trees and the Normalized Robinson-Foulds distance between the final estimated species tree (i.e., `roadies.nwk`) and the reference tree (i.e., REFERENCE parameter in `config.yaml`), for all iterations.
8. `time_stamps.csv`- this file contains the start time in first line, iteration number, number of gene trees required for estimating species tree, end time, and total runtime (in seconds), respectively, for all iterations in subsequent lines.

## Contributions

We welcome contributions from the community. If you encounter any issues or have suggestions for improvement, please open an issue on GitHub. For general inquiries and support, reach out to our team.

Anshu Gupta - ang037 [at] ucsd [dot] edu

Yatish Turakhia - yturakhia [at] ucsd [dot] edu

## Citing ROADIES

If you use ROADIES in your research or publications, please cite the following paper:

Gupta A, Mirarab S, Turakhia Y, (2024). Accurate, scalable, and fully automated inference of species trees from raw genome assemblies using ROADIES. _bioRxiv_. [https://www.biorxiv.org/content/10.1101/2024.05.27.596098v1](https://www.biorxiv.org/content/10.1101/2024.05.27.596098v1).

### Accessing ROADIES output files

The output files with the gene trees and species trees generated by ROADIES in the manuscript are deposited to [Dryad](https://datadryad.org/stash). To access it, please refer to the following:

Gupta, Anshu; Mirarab, Siavash; Turakhia, Yatish (2024). Accurate, scalable, and fully automated inference of species trees from raw genome assemblies using ROADIES [Dataset]. Dryad. [https://doi.org/10.5061/dryad.tht76hf73](https://doi.org/10.5061/dryad.tht76hf73).
