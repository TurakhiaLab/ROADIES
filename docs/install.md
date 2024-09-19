# Installation Methods

## Using ROADIES Bioconda package

To run ROADIES using Bioconda package, follow these steps:

**Note:** You need to have conda installed in your system. Also make sure you have updated version of glibc in your system (`GLIBC >= 2.29`).

To install and use conda in Ubuntu machine, execute the set of commands below:

```bash
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

```bash
conda create -n myenv python=3.9
conda activate myenv
```

2. Install ROADIES bioconda package

```
conda install roadies
```

All files of ROADIES along with dependencies will be found in `<conda_install_path>/miniconda3/envs/new_env/ROADIES`.

## Using DockerHub

To run ROADIES using DockerHub, follow these steps:

1. Pull the ROADIES Docker image from DockerHub:

```bash
docker pull ang037/roadies:latest
```
2. Run the Docker container:

```bash
docker run -it ang037/roadies:latest
```

## Using Docker locally

First, clone the repository (requires `git` to be installed in the system):

```bash
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
```

Then build and run the Docker container:

```bash
docker build -t roadies_image .
docker run -it roadies_image
```

## Using installation script (requires sudo access)

First clone the repository:

```bash
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
```

Then, execute the installation script:

```bash
chmod +x roadies_env.sh
source roadies_env.sh
```

This will install and build all tools and dependencies. Once the setup is complete, it will print `Setup complete` in the terminal and activate the `roadies_env` environment with all Conda packages installed. 

!!! Note 
    ROADIES is built on [Snakemake (workflow parallelization tool)](https://snakemake.readthedocs.io/en/stable/). It also requires various tools (PASTA, LASTZ, RAxML-NG, MashTree, FastTree, ASTRAL-Pro2) to be installed before performing the analysis. To ease the process, instead of individually installing the tools, we provide `roadies_env.sh` script to automatically download all dependencies into the user system.

### Required dependencies

To run this script, ensure the following dependencies are installed:
- Java Runtime Environment (version 1.7 or higher)
- Python (version 3 or higher)
- `wget` and `unzip` commands
- GCC (version 11.4 or higher)
- cmake (Download here: https://cmake.org/download/)
- Boost library (Download here: https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/)
- zlib (Download here: http://www.zlib.net/)
- GLIBC (Version 2.29 or higher)

For Ubuntu, you can install these dependencies with: 

```bash
sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git default-jre libgomp1 libboost-all-dev cmake
```

!!! Warning
    If you encounter issues with the Boost library, add its path to `$CPLUS_LIBRARY_PATH` and save it in `~/.bashrc`.
