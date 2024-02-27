# Installation

This section provides detailed instructions on how to install and set up the environment to run ROADIES in your system.

ROADIES is built on [Snakemake (workflow parallelization tool)](https://snakemake.readthedocs.io/en/stable/). It also requires various tools (PASTA, LASTZ, IQTREE, MashTree, FastTree, ASTRAL-Pro) to be installed before performing the analysis. To ease the process, instead of individually installing the tools, we provide various ways to automatically download all dependencies into the user system, as follows:

## Using Conda

Coming soon

## Using Docker

### Using DockerHub

Coming soon

### Running locally

First clone the repository

```
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
```

Then build and run docker

```
docker build roadies_image .
docker run -it roadies_image
```


## Using Installation scripts (requires sudo access)

First clone the repository

```
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
```
Then, execute bash script `roadies_env.sh` by following the commands below:

```
chmod +x roadies_env.sh
./roadies_env.sh
```
This will install and build all required tools and dependencies required by the user to get started. Once setup is complete, it will print `Setup complete` in the terminal. On its completion, a snakemake environment named `roadies_env` will be activated with all conda packages installed in it. 

Once this is done, follow the steps mentioned in [Usage](usage.md) tab for detailed configuration tailored to the specific requirements.