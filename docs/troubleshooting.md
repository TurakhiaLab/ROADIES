# Troubleshooting Steps

## 1. Mamba not found in the shell

When running the following command:
```bash
$ python ROADIES-main/run_roadies.py --cores 1
```
You may encounter this error:

```bash
rm: cannot remove 'output_files': No such file or directory
Unlocking working directory.
snakemake --cores 1 --config mode=accurate config_path=config/config.yaml num_threads=0 --use-conda --rerun-incomplete
Config file config/config.yaml is extended by additional config specified via the command line.
Building DAG of jobs...
CreateCondaEnvironmentException:
The 'mamba' command is not available in the shell /usr/bin/bash that will be used by Snakemake. You have to ensure that it is in your PATH, e.g., first activating the conda base environment with `conda activate base`.The mamba package manager (https://github.com/mamba-org/mamba) is a fast and robust conda replacement. It is the recommended way of using Snakemake's conda integration. It can be installed with `conda install -n base -c conda-forge mamba`. If you still prefer to use conda, you can enforce that by setting `--conda-frontend conda`.
```

### Cause

The `mamba` package manager is missing or not available in the environment.

### Solution

Install mamba:

```
conda install -n base -c conda-forge mamba
```

If you prefer using `conda`, you can enforce it by adding the `--conda-frontend` conda argument.

**Step 1:** In the downloaded ROADIES repository, open the file `noconverge.py` inside the `workflow` folder (`ROADIES/workflow/noconverge.py`).

**Step 2:** At line 31, add the argument `--conda-frontend conda` to the `cmd` command, as shown below:

```python
cmd = [
    "snakemake",
    "--cores",
    str(cores),
    "--config",
    "mode=" + str(mode),
    "config_path=" + str(config_path),
    "num_threads=" + str(num_threads),
    "--use-conda",
    "--rerun-incomplete",
    "--conda-frontend", "conda"
]
```
**Step 3:** Rerun the pipeline as follows:

```
python run_roadies.py --cores 16
```

## 2. Conda not recognized

### Cause

Conda is not added to your system's PATH.

### Solution

Ensure conda is added to the PATH by running the following commands:

```bash
export PATH="$HOME/miniconda3/bin:$PATH"
source ~/.bashrc
```

## 3. Singularity issues

### Cause

Problems arise when trying to run the pipeline with Singularity.

### Solution

We recommend using Docker instead of Singularity. Ensure Docker is installed and running on your system. We have also provided Bioconda support for users who face issues with Singularity.

## 4. Handling dependencies (glibc)

Ensure that the glibc version on your system is updated to 2.29 or higher. Update your system libraries if necessary. Otherwise you may encounter this error:

```bash
workflow/scripts/lastz_32: /lib64/libm.so.6: version 'GLIBC_2.29' not found
```

## 5. PASTA fails with insufficient core count

### Cause

Pasta fails when the number of cores is insufficient for the number of instances.

The pipeline provides `NUM_INSTANCES` as a configuration parameter in `config.yaml` to run multiple instances in parallel. Each instance can also be parallelized using threads. The number of threads per instance is calculated as:

```makefile
num_threads = number_of_cores / num_instances
```
If `num_instances > number_of_cores`, then `num_threads` will be 0 and the process (e.g., `pasta`) will fail.

### Solution

Ensure that the number of cores is greater than or equal to the number of instances. By default, `NUM_INSTANCES` is set to 4, so the number of cores (`--cores` in command line argument) must be at least 4. To run the pipeline with fewer cores, modify the `NUM_INSTANCES` parameter in the config file:

```bash
python run_roadies.py --cores <available_cores> --config_path config/config.yaml
```