#!/bin/bash

sudo apt-get update
sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git vim screen default-jre libgomp1 libboost-all-dev cmake

source /usr/share/miniconda3/etc/profile.d/conda.sh
source /usr/share/miniconda3/etc/profile.d/mamba.sh
conda activate roadies_env_test

wget -q https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip -O Linux.zip
unzip -q Linux.zip
cd ASTER-Linux
make
g++ -D CASTLES -std=gnu++11 -march=native -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro3
cd ..

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
wget https://github.com/lastz/lastz/archive/refs/heads/master.zip
unzip master.zip
cd lastz-master/src/
make lastz_32 flagsFor32="-Dmax_sequence_index=63 -Dmax_malloc_index=40 -Ddiag_hash_size=4194304"
make install_32
cp lastz_32 ../../
cd ../../../../

echo "Setup complete. Remember to source before running your projects."

