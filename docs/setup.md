# Installation

This section provides detailed instructions on how to install and set up the environment to run ROADIES in your system.

ROADIES is built on [Snakemake (workflow parallelization tool)](https://snakemake.readthedocs.io/en/stable/). It also requires various tools (PASTA, LASTZ, IQTREE, MashTree, FastTree, ASTRAL-Pro) to be installed before performing the analysis. To ease the process, instead of individually installing the tools, we provide a script to automatically download all dependencies into the user system. 

## Using Linux

Execute bash script `roadies_env.sh` by following the commands below:

```
chmod +x roadies_env.sh
./roadies_env.sh
```
This will install and build all required tools and dependencies required by the user to get started. Once setup is complete, it will print `Setup complete` in the terminal. On its completion, a snakemake environment named `roadies_env` will be activated with all conda packages installed in it. 

Once this is done, follow the steps mentioned in [Getting Started](start.md) page for quick run. For detailed configuration tailored to the specific requirements, follow the steps mentioned in [Usage](usage.md).

## Using Docker