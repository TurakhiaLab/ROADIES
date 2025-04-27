# Installation Methods

Please follow any of the options below to install ROADIES in your system. 

## Option 1: Install via Bioconda (Recommended)

1. Install Conda (if not installed):

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh

export PATH="$HOME/miniconda3/bin:$PATH"
source ~/.bashrc
```

2. Configure Conda channels:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Verify the installation by running `conda` in your terminal

3. Create and activate a custom environment:

```bash
conda create -n roadies_env python=3.9 ete3 seaborn
conda activate roadies_env
```

4. Install ROADIES:

```bash
conda install roadies
```

5. Locate the installed files:

```bash
cd $HOME/miniconda3/envs/roadies_env/ROADIES

```

Now you are ready to follow the Quick Start section to run the pipeline. 

## Option 2: Install via DockerHub

If you would like to install ROADIES using DockerHub, follow these steps:

1. Pull the ROADIES image from DockerHub:

```bash
docker pull ang037/roadies:latest
```
2. Launch a container:

```bash
docker run -it ang037/roadies:latest
```

Once you are able to access the ROADIES repository, refer to the Quick Start section to run the pipeline. 

## Option 3: Install via Local Docker Build

1. Clone the ROADIES repository:

```bash
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
```

2. Build and run the Docker container:

```bash
docker build -t roadies_image .
docker run -it roadies_image
```

Once you are able to access the ROADIES repository, refer to Quick Start instructions to run the pipeline. 

## Option 4: Install via Source Script

1. Install the following dependencies (**requires sudo access**):

- Java Runtime Environment (Version 1.7 or higher)
- Python (Version 3.9 or higher)
- `wget` and `unzip` commands
- GCC (Version 11.4 or higher)
- cmake (Download here: https://cmake.org/download/)
- Boost library (Download here: https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/)
- zlib (Download here: http://www.zlib.net/)

For Ubuntu, you can install these dependencies with: 

```bash
sudo apt-get install -y wget unzip make g++ python3 python3-pip python3-setuptools git default-jre libgomp1 libboost-all-dev cmake
```

2. Clone the repository:

```bash
git clone https://github.com/TurakhiaLab/ROADIES.git
cd ROADIES
```

3. Run the installation script:

```bash
chmod +x roadies_env.sh
source roadies_env.sh
```

After successful setup (Setup complete message), your environment roadies_env will be activated. Proceed to Quick Start.

**Note:** If you encounter issues with the Boost library, add its path to `$CPLUS_LIBRARY_PATH` and save it in `~/.bashrc`.