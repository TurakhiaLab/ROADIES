# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

USER root

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y openssh-client git wget unzip make g++ python3 python3-pip \
        python3-setuptools vim screen default-jre libgomp1 libboost-all-dev cmake && \
    rm -rf /var/lib/apt/lists/*

# Clone the ROADIES repository
RUN git clone https://github.com/TurakhiaLab/ROADIES.git /ROADIES

WORKDIR /ROADIES

# Run the environment setup at build time so the image is ready to use immediately
SHELL ["/bin/bash", "-lc"]
RUN chmod +x roadies_env.sh && bash roadies_env.sh

# Ensure conda is sourced and the environment is activated for all subsequent RUN steps
ENV PATH="/root/conda/envs/roadies_env/bin:/root/conda/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=roadies_env

ENTRYPOINT ["/bin/bash", "-lc", "source /root/conda/etc/profile.d/conda.sh && conda activate roadies_env && exec bash"]
