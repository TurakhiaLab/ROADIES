# wga-phylo
Snakemake and Snakedeploy are best installed via the Mamba package manager (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via Mambaforge. 

Given that Mamba is installed, run

mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
to install both Snakemake and Snakedeploy in an isolated environment. For all following commands ensure that this environment is activated via
conda activate snakemake
SnakeMake workflow for sequence selection to lastz, change your configuration data in configuration file

`snakemake --core [number of cores] --use-conda`

Snakemake will automatically detect the main Snakefile in the workflow subfolder and execute the workflow module that has been defined by the deployment in step 2.

In order to run the pipeline, a manual installation of ASTRAL-PRO(latest version) is required. Replace the ASTER-Linux directory with the ASTRAL-PRO git repository. 

rmdir ASTER-Linux\
wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip \
unzip Linux.zip
