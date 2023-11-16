# Use Ubuntu 22.09 as the base image
FROM ubuntu:22.04

USER root
# Install necessary tools
RUN apt-get update && \
    apt-get install -y openssh-client git

WORKDIR /root

ARG ssh_prv_key
RUN mkdir /root/.ssh && \
    chmod -R 700 /root/.ssh && \
    echo "${ssh_prv_key}" >> /root/.ssh/id_rsa

RUN echo ${ssh_prv_key}
RUN chmod 0400 /root/.ssh/id_rsa && echo "StrictHostKeyChecking no" > /root/.ssh/config
RUN cat /root/.ssh/id_rsa
RUN ssh-keyscan github.com >> /root/.ssh/known_hosts

RUN git clone -v git@github.com:TurakhiaLab/wga-phylo.git && \
    rm /root/.ssh/id_rsa*

WORKDIR wga-phylo

# Set environment variables
ENV HOME=/root

RUN apt-get update && \
    apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git vim screen default-jre libgomp1 libboost-all-dev cmake

# Download and install Mambaforge
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash Mambaforge-Linux-x86_64.sh -b -p "${HOME}/conda"

# Source Conda and Mamba scripts
RUN echo "source ${HOME}/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "source ${HOME}/conda/etc/profile.d/mamba.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    echo "mamba create -y -c conda-forge -c bioconda --name roadies_env snakemake alive-progress biopython iqtree=2.2.0.3 numpy lastz mashtree matplotlib seaborn treeswift=1.1.28 fasttree=2.1.11 python=3.11 ete3" >> ~/.bashrc && \
    echo "conda activate roadies_env" >> ~/.bashrc

# Download ASTER repository
RUN wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip && \
    unzip Linux.zip && \
    cd ASTER-Linux && \
    make && \
    cd ..

# Clone PASTA repository and install
RUN git clone https://github.com/smirarab/pasta.git

# Clone Sate-Tools repository
RUN git clone https://github.com/smirarab/sate-tools-linux.git
    
RUN echo "cd pasta" >> ~/.bashrc && \
    echo "python3 setup.py develop --user" >> ~/.bashrc && \
    echo "cd .." >> ~/.bashrc

#build sampling code
RUN echo "cd workflow/scripts/sampling" >> ~/.bashrc && \ 
    echo "mkdir build" >> ~/.bashrc && \
    echo "cd build" >> ~/.bashrc && \
    echo "cmake .." >> ~/.bashrc && \
    echo "make" >> ~/.bashrc && \
    echo "cd ../../../.." >> ~/.bashrc

# Source the updated .bashrc to activate Conda environment
CMD ["/bin/bash"]

