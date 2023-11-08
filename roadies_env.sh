#!/bin/bash

# Download and install Mambaforge
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh -b -p "${HOME}/conda"

# Source Conda and Mamba scripts
source ${HOME}/conda/etc/profile.d/conda.sh
source ${HOME}/conda/etc/profile.d/mamba.sh
conda activate base
mamba create -y -c conda-forge -c bioconda --name roadies_env snakemake alive-progress biopython iqtree=2.2.0.3 numpy lastz mashtree matplotlib seaborn treeswift=1.1.28 fasttree=2.1.11 python=3.11 ete3
echo "conda activate roadies_env" >> ~/.bashrc

# Download ASTER repository
wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip
unzip Linux.zip
cd ASTER-Linux
make
cd ..

# Clone PASTA repository and install
git clone https://github.com/smirarab/pasta.git

# Clone Sate-Tools repository
git clone https://github.com/smirarab/sate-tools-linux.git
    
cd pasta
python3 setup.py develop --user
cd ..

#build sampling code
cd workflow/scripts/sampling
mkdir build
cd build
cmake ..
make
cd ../../../..

# Source the updated .bashrc to activate Conda environment
source ~/.bashrc

echo "Setup complete"
