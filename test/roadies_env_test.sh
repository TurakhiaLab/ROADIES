#!/bin/bash

sudo apt-get update
sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git vim screen default-jre libgomp1 libboost-all-dev cmake

wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-Linux-x86_64.sh"
bash Miniforge3.sh -b -p "${HOME}/conda"

source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"
conda activate roadies_env_test

wget -q https://github.com/chaoszhang/ASTER/archive/refs/tags/v1.19.zip -O Linux.zip
unzip -q Linux.zip
mv ASTER-1.19 ASTER-Linux
cd ASTER-Linux
make
g++ -D CASTLES -std=gnu++11 -march=native -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro3
cd ..
cp ASTER-Linux/bin/astral-pro3 .

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

# Build LASTZ
cd workflow/scripts
wget https://github.com/lastz/lastz/archive/refs/tags/1.04.45.zip
unzip 1.04.45.zip 
cd lastz-1.04.45/src/
make lastz_40
make install_40
cp lastz_40 ../../../../
cd ../../../../

echo "Setup complete. Remember to source before running your projects."

