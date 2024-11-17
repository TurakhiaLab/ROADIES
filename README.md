<div align="center">
    
# Reference-free Orthology-free Annotation-free DIscordance aware Estimation of Species tree (ROADIES)

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/TurakhiaLab/ROADIES/blob/main/LICENSE

[![License][license-badge]][license-link]
[![Build Status](https://github.com/TurakhiaLab/ROADIES/actions/workflows/ci.yml/badge.svg)](https://github.com/TurakhiaLab/ROADIES/actions)
[<img src="https://img.shields.io/badge/Made with-Snakemake-green.svg?logo=snakemake">](https://snakemake.readthedocs.io/en/v7.19.1/index.html)
[<img src="https://img.shields.io/badge/Install with-Biooconda-brightgreen.svg?logo=conda">](http://bioconda.github.io/recipes/roadies/README.html)
[<img src="https://img.shields.io/badge/Install with-DockerHub-informational.svg?logo=Docker">](https://hub.docker.com/r/ang037/roadies)
[<img src="https://img.shields.io/badge/Submitted to-bioRxiv-critical.svg?logo=LOGO">](https://www.biorxiv.org/content/10.1101/2024.05.27.596098v1)
[<img src="https://img.shields.io/badge/DOI-10.5061/dryad.tht76hf73-yellowgreen.svg?logo=LOGO">](https://doi.org/10.5061/dryad.tht76hf73)
[<img src="https://img.shields.io/badge/Watch it on-Youtube-FF0000.svg?logo=YouTube">](https://youtu.be/1sR741TvZnM?si=vVNAnonvzNEzrLKq)

<div align="center">

<img src="docs/images/ROADIES_logo.png" style="height: 350px; width: auto;">

</div>

</div>

## Table of Contents
- [Introduction](#overview)
- [Quick Install](#usage)
    - [Using ROADIES Bioconda package](#conda)
    - [Using DockerHub](#dockerhub)
    - [Using Docker locally](#docker)
    - [Using Installation Script](#script)
- [Quick Start](#start)
- [Run ROADIES with your own datasets](#runpipeline)
- [Citing ROADIES](#citation)

<br>

## <a name="overview"></a> Introduction

Welcome to the official repository of ROADIES, a novel pipeline designed for phylogenetic tree inference of the species directly from their raw genomic assemblies. ROADIES offers a fully automated, easy-to-use, scalable solution, eliminating any manual steps and providing unique flexibility in adjusting the tradeoff between accuracy and runtime. 

**For more detailed information on all the features and settings of ROADIES, please refer to our [Wiki](https://turakhialab.github.io/ROADIES/).**

<br>

<div align="center">

  <figure>
    <img src="docs/images/drawing_github.png" alt="ROADIES Pipeline Stages">
    <figcaption>Figure: ROADIES Pipeline Stages</figcaption>
  </figure>

</div>

<br>

## <a name="usage"></a> Quick Install

### <a name="conda"></a> Using ROADIES Bioconda package (recommended)

To run ROADIES using Bioconda package, follow these steps:

**Note:** You need to have conda installed in your system. Also make sure you have updated version of glibc in your system (`GLIBC >= 2.29`).

To install and use conda in Ubuntu machine, execute the set of commands below:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh

export PATH="$HOME/miniconda3/bin:$PATH"
source ~/.bashrc

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

After this, try running `conda` in your terminal to check if conda is properly installed. Once it is installed, follow the steps below:

1. Create and activate custom conda environment with Python version 3.9

```
conda create -n myenv python=3.9
conda activate myenv
```

2. Install ROADIES bioconda package

```
conda install roadies
```

All files of ROADIES along with dependencies will be found in `<conda_install_path>/miniconda3/envs/myenv/ROADIES`.

### <a name="dockerhub"></a> Using DockerHub

To run ROADIES using DockerHub, follow these steps:

1. Pull the ROADIES Docker image from DockerHub:

```
docker pull ang037/roadies:latest
```
2. Run the Docker container:

```
docker run -it ang037/roadies:latest
```

### <a name="docker"></a> Using Docker locally

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

### <a name="script"></a> Using installation script (requires sudo access)

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

#### Required dependencies

To run this script, ensure the following dependencies are installed:
- Java Runtime Environment (Version 1.7 or higher)
- Python (Version 3.9 or higher)
- `wget` and `unzip` commands
- GCC (Version 11.4 or higher)
- cmake (Download here: https://cmake.org/download/)
- Boost library (Download here: https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/)
- zlib (Download here: http://www.zlib.net/)
- GLIBC (Version 2.29 or higher)

For Ubuntu, you can install these dependencies with: 

```
sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git default-jre libgomp1 libboost-all-dev cmake
```

**Note:** If you encounter issues with the Boost library, add its path to `$CPLUS_LIBRARY_PATH` and save it in `~/.bashrc`.

<br>

## <a name="start"></a> Quick Start

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

The second command will download the 11 Drosophila genomic datasets (links provided in `test/input_genome_links.txt`) and save them in the `test/test_data` directory. The third command will run ROADIES pipeline for those 11 Drosophila genomes and save the final newick tree as `roadies.nwk` in a separate `output_files` folder upon completion.

<br>

## <a name="runpipeline"></a> Run ROADIES with your own datasets

To run ROADIES with your own datasets, follow these steps:

1. **Specify Input Genomic Dataset**: Update the `config.yaml` file (found in the ROADIES directory - `config` folder) to include the path to your input datasets under the `GENOMES` parameter. Ensure all input genomic assemblies are in `.fa` or `.fa.gz` format and named according to the species' name (e.g., `Aardvark.fa`). 

**Note**: Each file should contain the genome assembly of one unique species. If a file contains multiple species, split it into individual genome files (`fasplit` can be used: `faSplit byname <input_dir> <output_dir>`).

2. **Configure Other Parameters**: Adjust other parameters in `config.yaml` as needed. Detailed information on each parameter is available in the [`Usage` section](https://turakhialab.github.io/ROADIES/).

3. **Run the Pipeline**: Execute the pipeline with the following command (example for 16 cores):

```
python run_roadies.py --cores 16
```

The output species tree in Newick format will be saved as `roadies.nwk` in the `output_files` folder.

4. **Modes of operation**: ROADIES supports multiple modes of operation (`fast`, `balanced`, `accurate`) by controlling the accuracy-runtime tradeoff. Use any one of the following commands to select a mode (`accurate` mode is the default):


```
python run_roadies.py --cores 16 --mode accurate

python run_roadies.py --cores 16 --mode balanced

python run_roadies.py --cores 16 --mode fast
```

<br>

### For troubleshooting and contribution details, refer to [Wiki](https://turakhialab.github.io/ROADIES/)

<br>

## <a name="citation"></a> Citing ROADIES

If you use ROADIES in your research or publications, please cite the following paper:

Gupta A, Mirarab S, Turakhia Y, (2024). Accurate, scalable, and fully automated inference of species trees from raw genome assemblies using ROADIES. _bioRxiv_. [https://www.biorxiv.org/content/10.1101/2024.05.27.596098v1](https://www.biorxiv.org/content/10.1101/2024.05.27.596098v1).

### Accessing ROADIES output files

The output files with the gene trees and species trees generated by ROADIES in the manuscript are deposited to [Dryad](https://datadryad.org/stash). To access it, please refer to the following:

Gupta, Anshu; Mirarab, Siavash; Turakhia, Yatish (2024). Accurate, scalable, and fully automated inference of species trees from raw genome assemblies using ROADIES [Dataset]. Dryad. [https://doi.org/10.5061/dryad.tht76hf73](https://doi.org/10.5061/dryad.tht76hf73).


