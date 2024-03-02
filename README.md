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

#### Using installation script

First clone the repository, as follows (requires `git` to be installed in the system):

```
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
```

Then, execute the bash script `roadies_env.sh` by following the commands below (**Warning:** check the dependencies below before running this script):

```
chmod +x roadies_env.sh
source roadies_env.sh
```

This will install and build all tools and dependencies required by the user to get started. Once setup is complete, it will print `Setup complete` in the terminal. On its completion, a snakemake environment named `roadies_env` will be activated with all conda packages installed in it. Now you are ready to run our pipeline (follow [Run ROADIES pipeline](index.md#Run-ROADIES-pipeline) section).

##### Required dependencies

To run this script, user should have the following installations:
- Java Runtime Environment (version 1.7 or higher)
- Python (version 3 or higher)
- `wget` and `unzip` commands
- GCC (version 11.4 or higher)
- cmake command: https://cmake.org/download/
- Boost library: https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/ and zlib http://www.zlib.net/ are required when running cmake and make.

**Note:** The current version of ROADIES is extensively tested with Linux environment only. For Ubuntu, to install above dependencies, please run the following command OR uncomment the initial lines of `roadies_env.sh` file. 

```
sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git default-jre libgomp1 libboost-all-dev cmake
```

**Note:** As a non-root user, the `make` command won't work because these libraries hasn't configured to an environment variable. You have to add your boost library path into `$CPLUS_LIBRARY_PATH` and save it into `~/.bashrc`, then gcc will be able to find `boost/program_option.hpp`. All these requirement only work in a version of gcc which greater than 7.X (or when running `make`, it will report error: `unrecognized command line option '-std=c++17‘!` ).


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

Once setup is done, run the following commands for 16-core machine:


```
mkdir -p test/test_data && cat test/input_genome_links.txt | xargs -I {} sh -c 'wget -O test/test_data/$(basename {}) {}'

python run_roadies.py --cores 16
```

The first line will download the 11 Drosophila genomic datasets (links are provided in `test/input_genome_links.txt`) and save it in `test/test_data` directory. Second line will run ROADIES for those 11 Drosophila genomes and save the final newick tree as `roadies.nwk` in a separate `ROADIES/output_files` folder after the completion.

## <a name="runpipeline"></a> Run ROADIES with your own datasets

To run ROADIES with your own datasets,follow the steps below:

### Specify input genomic dataset

Specify the path of the input genomic dataset in `config.yaml` file (`GENOMES` parameter).

**Note**: All input genome assemblies should be in `.fa` or `.fa.gz` format. The genome assembly files should be named according to the species' names (for example, Aardvark's genome assembly is to be named `Aardvark.fa`). Each file should contain the genome assembly of one unique species. If a file contains multiple species, split it into individual genome files (fasplit can be used for this: `faSplit byname <input_dir> <output_dir>`)

### Configure other parameters

Configure other parameters in `config.yaml` file based on your use-case requirements. The detailed information of all parameters are mentioned in the `Usage` section [here](https://turakhialab.github.io/ROADIES/).

### Run the pipeline

After modifying the configurations, run the following command to execute ROADIES pipeline with 16 cores:

```
python run_roadies.py --cores 16
```

After the completion of the execution, the output species tree in Newick format will be saved as `roadies.nwk` in a separate `output_files` folder.


**Modes of operation**: ROADIES also supports multiple modes of operation (`fast`, `balanced`, `accurate`) by controlling the accuracy-runtime tradeoff. Try the following commands for various modes of operation (`accurate` mode is the default mode)


```
python run_roadies.py --cores 16 --mode accurate
```

```
python run_roadies.py --cores 16 --mode balanced
```

```
python run_roadies.py --cores 16 --mode fast
```

## <a name="support"></a> Contributions and Support

We welcome contributions from the community to enhance the capabilities of ROADIES. If you encounter any issues or have suggestions for improvement, please open an issue on GitHub. For general inquiries and support, reach out to our team.

## <a name="citation"></a> Citing ROADIES

If you use the ROADIES pipeline for species tree inference in your research or publications, we kindly request that you cite the following paper:



