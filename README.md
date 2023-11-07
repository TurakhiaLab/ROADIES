 <div align="center">

![ROADIES_logo](https://github.com/TurakhiaLab/wga-phylo/assets/114828525/05cd206e-542c-4ee4-bfd6-d4c03fed5984)

# Reference free Orthology free Alignment free DIscordant Aware Estimation of Species Tree (ROADIES)

</div>

## Table of Contents
- [Introduction](#overview)
- [Using ROADIES](#usage)
- [Contributions and Support](#support)
- [Citing ROADIES](#citation)

## <a name="overview"></a> Introduction

Welcome to the official repository of ROADIES, a novel pipeline designed for phylogenetic tree inference of the species directly from their raw genomic assemblies. ROADIES pipeline offers a fully automated, easy-to-use, scalable solution, eliminating any error-prone manual steps and providing unique flexibility in adjusting the tradeoff between accuracy and runtime. 
<br>

<div align="center">

<img src="drawing_github.png">

</div>

## <a name="usage"></a> Quick start

This section provides brief overview on how to get started with the tool. To know more details about all the exisiting features and settings, please read [this documentation](https://turakhialab.github.io/wga-phylo/).

To install and build all required tools and dependencies to get started, execute the bash script `roadies_env.sh` by following the commands below:

```
chmod +x roadies_env.sh
./roadies_env.sh
```

Once setup is complete, it will print `Setup complete` in the terminal. On its completion, a snakemake environment named `roadies_env` will be activated with all conda packages installed in it. 

To run the ROADIES pipeline with 32 cores, specify the path to your input files as `--GENOMES` in `config/config.yaml`. Then, run the following command:

```
python run_roadies.py --cores 32
```
The output species tree will be saved as `roadies.nwk` in a separate results folder. 

**Note**: All input genome assemblies should be in `.fa` or `.fa.gz` format. The genome assembly files should be named according to the species' names (for example, Aardvark's genome assembly is to be named `Aardvark.fa`). Each file should contain the genome assembly of one unique species. If a file contains multiple species, split it into individual genome files (fasplit can be used for this: `faSplit byname <input_dir> <output_dir>`)

**Modes of operation**: ROADIES also supports multiple modes of operation (`fast`, `balanced`, `accurate`) by controlling the accuracy-runtime tradeoff. Try the following commands for various modes of operation (`accurate` mode is the default mode)


```
python run_roadies.py --cores 32 --mode accurate
```

```
python run_roadies.py --cores 32 --mode balanced
```

```
python run_roadies.py --cores 32 --mode fast
```

## <a name="support"></a> Contributions and Support

We welcome contributions from the community to enhance the capabilities of ROADIES. If you encounter any issues or have suggestions for improvement, please open an issue on GitHub. For general inquiries and support, reach out to our team.

## <a name="citation"></a> Citing ROADIES

If you use the ROADIES pipeline for species tree inference in your research or publications, we kindly request that you cite the following paper:



