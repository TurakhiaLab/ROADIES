#!/bin/bash

sudo apt-get update
sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git vim screen default-jre libgomp1 libboost-all-dev cmake

wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-Linux-x86_64.sh"
bash Miniforge3.sh -b -p "${HOME}/conda"

source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"
mamba create -y -c conda-forge -c bioconda --name roadies_env_test snakemake alive-progress biopython iqtree=2.2.0.3 numpy lastz mashtree matplotlib seaborn treeswift=1.1.28 fasttree=2.1.11 python=3.11 raxml-ng ete3 lastz=1.04.52 aster=1.19 pyyaml seaborn
conda activate roadies_env_test

git clone https://github.com/smirarab/pasta.git

git clone https://github.com/smirarab/sate-tools-linux.git

cd pasta
python3 setup.py develop --user
cd ..

# Build sampling code
cd workflow/scripts/sampling
mkdir build
cd build
cmake ..
make
cd ../../../..

echo "Setup complete. Remember to source before running your projects."

