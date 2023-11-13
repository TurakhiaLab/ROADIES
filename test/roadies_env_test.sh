#!/bin/bash

sudo apt-get update
sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git vim screen default-jre libgomp1 libboost-all-dev cmake

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

echo "Setup complete"
