#!/bin/bash

# Install necessary tools
sudo apt-get update
sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git vim screen default-jre libgomp1 libboost-all-dev cmake

# Download and install Mambaforge
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh -b -p "${HOME}/conda"

# Source Conda and Mamba scripts
echo "source ${HOME}/conda/etc/profile.d/conda.sh" >> ~/.bashrc
echo "source ${HOME}/conda/etc/profile.d/mamba.sh" >> ~/.bashrc
echo "conda activate base" >> ~/.bashrc
echo "mamba create -y -c conda-forge -c bioconda --name roadies_env snakemake alive-progress biopython iqtree=2.2.0.3 numpy lastz mashtree matplotlib seaborn treeswift=1.1.28" >> ~/.bashrc
echo "conda activate roadies_env" >> ~/.bashrc

pip3 install ete3

# Download ASTER repository
wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip
unzip Linux.zip
cd ASTER-Linux
make
cd ..

# Clone PASTA repository and install
RUN git clone https://github.com/smirarab/pasta.git

# Clone Sate-Tools repository
RUN git clone https://github.com/smirarab/sate-tools-linux.git
    
cd pasta
python setup.py develop --user
cd ..

# Source the updated .bashrc to activate Conda environment
source ~/.bashrc

echo "Setup complete"
