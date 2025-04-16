#!/bin/bash

# Required installations: (uncomment next 2 lines if you have sudo access, otherwise make sure following tools are installed before proceeding)
# sudo apt-get update
# sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git vim screen default-jre libgomp1 libboost-all-dev cmake

# Define installation paths and check for directory
CONDA_PATH="${HOME}/conda"
ROADIES_ENV_SETUP="roadies_env.sh"

# Download and install Mambaforge if not already installed
if [ ! -d "${CONDA_PATH}" ]; then
    wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-Linux-x86_64.sh"
    bash Miniforge3.sh -b -p "${CONDA_PATH}"
fi

# Create and setup the Conda environment if it doesn't exist
source "${CONDA_PATH}/etc/profile.d/conda.sh"  # Temporarily source conda for this script
source "${CONDA_PATH}/etc/profile.d/mamba.sh"  # Temporarily source mamba for this script

if ! conda env list | grep -q "roadies_env"; then
    mamba create -y -c conda-forge -c bioconda --name roadies_env snakemake alive-progress biopython iqtree=2.2.0.3 numpy lastz mashtree matplotlib seaborn treeswift=1.1.28 fasttree=2.1.11 python=3.11 raxml-ng ete3 lastz=1.04.52 aster=1.19 pyyaml seaborn
fi
conda activate roadies_env

# Clone PASTA if not already done
if [ ! -d "pasta" ]; then
    git clone https://github.com/smirarab/pasta.git
fi

# Clone Sate-Tools if not already done
if [ ! -d "sate-tools-linux" ]; then
    git clone https://github.com/smirarab/sate-tools-linux.git
fi

# Setup PASTA if not already done
if [ -d "pasta" ]; then
    mafft_file="pasta/bin/mafft"
    if [ ! -f "$mafft_file" ]; then
        cd pasta
        python3 setup.py develop --user
        cd ..
    fi
fi

# Build sampling code
if [ ! -d "workflow/scripts/sampling/build" ]; then
    cd workflow/scripts/sampling
    mkdir build
    cd build
    cmake ..
    make
    cd ../../../..
fi

# Install ete3 for the user, only if it is not already installed
python3 -m pip show ete3 &>/dev/null || python3 -m pip install --user ete3

echo "Setup complete. Remember to source '${ROADIES_ENV_SETUP}' before running your projects."
