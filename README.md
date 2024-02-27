 <div align="center">

![ROADIES_logo](https://github.com/TurakhiaLab/ROADIES/assets/114828525/05cd206e-542c-4ee4-bfd6-d4c03fed5984)

# Reference-free Orthology-free Alignment-free DIscordance aware Estimation of Species tree (ROADIES)

</div>

## Table of Contents
- [Introduction](#overview)
- [Quick Start](#usage)
- [Run ROADIES for your own dataset](#runpipeline)
- [Contributions and Support](#support)
- [Citing ROADIES](#citation)

## <a name="overview"></a> Introduction

Welcome to the official repository of ROADIES, a novel pipeline designed for phylogenetic tree inference of the species directly from their raw genomic assemblies. ROADIES pipeline offers a fully automated, easy-to-use, scalable solution, eliminating any error-prone manual steps and providing unique flexibility in adjusting the tradeoff between accuracy and runtime. 
<br>

<div align="center">

<img src="drawing_github.png">

</div>

## <a name="usage"></a> Quick start

This section provides brief overview on how to get started with the tool. To know more details about all the exisiting features and settings, please read [this documentation](https://turakhialab.github.io/ROADIES/).

### Quick install

#### Using installation script (requires sudo access)

First clone the repository, as follows:

```
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
```

Then, execute the bash script `roadies_env.sh` by following the commands below:

```
chmod +x roadies_env.sh
./roadies_env.sh
```

Once setup is complete, it will print `Setup complete` in the terminal. On its completion, a snakemake environment named `roadies_env` will be activated with all conda packages installed in it. 

#### Using docker

First clone the repository

```
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
```

Then build and run docker

```
docker build roadies_image .
docker run -it roadies_image
```

### Run ROADIES pipeline

Once setup is done, run the following command for 16-core machine:

```
python run_roadies.py --cores 16
```

This will run ROADIES for 11 Bacterial genomes, whose genomic sequences are provided in `test/test_data` folder. After the completion, the final newick tree for these 11 species will be saved as `roadies.nwk` in a separate `ROADIES/output_files` folder.

## <a name="runpipeline"></a> Run ROADIES with your own datasets

To run ROADIES with your own datasets,follow the steps below:

### Specify input genomic dataset

Specify the path of the input genomic dataset in `config.yaml` file (`GENOMES` parameter).

**Note**: All input genome assemblies should be in `.fa` or `.fa.gz` format. The genome assembly files should be named according to the species' names (for example, Aardvark's genome assembly is to be named `Aardvark.fa`). Each file should contain the genome assembly of one unique species. If a file contains multiple species, split it into individual genome files (fasplit can be used for this: `faSplit byname <input_dir> <output_dir>`)

### Configure other parameters

Configure other parameters in `config.yaml` file based on your use-case requirements. The detailed information of all parameters are mentioned [here](https://turakhialab.github.io/ROADIES/usage).

### Run the pipeline

After modifying the configurations, run the following command to execute ROADIES pipeline with 16 cores:

```
python run_roadies.py --cores 32
```

After the completion of the execution, the output species tree in Newick format will be saved as `roadies.nwk` in a separate `output_files` folder.


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



