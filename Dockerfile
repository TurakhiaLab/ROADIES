# Use Ubuntu 22.09 as the base image
FROM ubuntu:22.04

USER root

RUN apt-get update && \
    apt-get install -y openssh-client git wget unzip make g++ python3 python3-pip python3-setuptools git vim screen default-jre libgomp1 libboost-all-dev cmake

# Clone the ROADIES repository
RUN git clone https://github.com/TurakhiaLab/ROADIES.git

WORKDIR ROADIES

RUN chmod +x roadies_env.sh

# Conda environment activation and source the setup script on container start
RUN echo "source roadies_env.sh" >> ~/.bashrc && \
    echo "source ${HOME}/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate roadies_env" >> ~/.bashrc

# Start the container with bash
CMD ["/bin/bash"]
