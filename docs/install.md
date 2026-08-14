# Installation Methods

Please follow any of the options below to install ROADIES on your system.

!!! Note
    ROADIES requires Linux (tested on Ubuntu 20.04/22.04). GPU features (ROADIES_XP, optional) additionally need an NVIDIA GPU with CUDA.

## Option 1: Install via Bioconda (Recommended)

1. Install Conda, if you don't already have it:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
export PATH="$HOME/miniconda3/bin:$PATH" && source ~/.bashrc
```

2. Create an environment and install ROADIES into it:

```bash
conda create -n roadies_env -c bioconda -c conda-forge python=3.9 ete3 seaborn
conda activate roadies_env
conda install roadies=0.1.10
```

3. `conda install` puts the full repository contents (Snakemake rules, scripts, `config.yaml`, `run_roadies.py`, etc.) under `$CONDA_PREFIX/ROADIES` — that's your working directory from now on:

```bash
cd $CONDA_PREFIX/ROADIES
```

You're now ready for Quick Start — run it from this `$CONDA_PREFIX/ROADIES` directory (`cd ROADIES` if you've since moved elsewhere and need to get back).

!!! Note
    Bioconda's `pasta` package (the default multiple-sequence aligner) already installs `run_pasta.py`/`run_seqtools.py` correctly — no manual PASTA build should be needed. If the pipeline does fail on the PASTA step, see [Troubleshooting: Issues with PASTA](troubleshooting.md#error-1-issues-with-pasta).

## Option 2: Install via DockerHub

1. Pull and run the prebuilt image:

```bash
docker pull ang037/roadies:latest
docker run -it ang037/roadies:latest
```

This launches an interactive container with the `roadies_env` conda environment already active and the working directory set to the ROADIES repository. Proceed to Quick Start.

## Option 3: Install via Local Docker Build

1. Clone the repository and build the image:

```bash
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
docker build -t roadies_image .
docker run -it roadies_image
```

Proceed to Quick Start once you're inside the container.

## Option 4: Install via Source Script

1. Install the system dependencies (**requires sudo access**): Java Runtime Environment (1.7+), Python (3.9+), `wget`/`unzip`, GCC (11.4+), [cmake](https://cmake.org/download/), [Boost](https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/), [zlib](http://www.zlib.net/). On Ubuntu:

```bash
sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git default-jre libgomp1 libboost-all-dev cmake
```

2. Clone the repository and run the setup script:

```bash
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
source roadies_env.sh
```

`roadies_env.sh` creates and activates the `roadies_env` conda environment, then builds everything ROADIES needs from source — PASTA, and (for ROADIES_XP) TWILIGHT plus, best-effort, MLIPPER if CUDA and libpll are already present. A `Setup complete` message means you're ready for Quick Start.

!!! Note
    If you encounter issues with the Boost library, add its path to `$CPLUS_LIBRARY_PATH` and save it in `~/.bashrc`.

!!! Note
    No extra steps are needed to use `--mode placement`; GPU placement (`--gpu`) does require a CUDA-capable GPU. See the [ROADIES_XP guide](roadies_xp.md#gpu-placement) for details, and re-run `roadies_env.sh` (or `bash MLIPPER/install/setup_host.sh`) once CUDA/libpll are available if MLIPPER was skipped at first setup.
