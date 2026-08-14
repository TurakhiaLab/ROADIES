# Troubleshooting Steps

## Error 1. Issues with PASTA

### Solution

Bioconda's `pasta` package (>=1.9.0) normally installs `run_pasta.py`/`run_seqtools.py` correctly on its own - the rules in `workflow/rules/multi_align.smk` already call those names directly, so this shouldn't come up in a clean install. If it still does (e.g. `run_pasta.py`/`run_seqtools.py` not found, or resolving to the wrong install), build PASTA from source instead. Run the following from the main ROADIES repository directory (after doing `cd ROADIES`), within the activated Conda environment:

```bash
git clone https://github.com/smirarab/pasta.git
git clone https://github.com/smirarab/sate-tools-linux.git
cd pasta
python3 setup.py develop --user
```

This installs `run_pasta.py`/`run_seqtools.py` as scripts under `~/.local/bin`, which takes priority on `PATH` over any conda environment's own copy - including in *other* conda environments on the same machine later on. If PASTA behaves oddly after switching environments or machines, check for a stray `~/.local/bin/run_pasta.py`/`run_seqtools.py` from a past run of this workaround before assuming something else is wrong.

After doing this change, please re-run the ROADIES pipeline.

## Error 2. Environment conflict

### Solution 

If you encounter the following error message - `"ls: relocation error: /lib64/libacl.so.1: symbol getxattr, version ATTR_1.0 not defined in file libattr.so.1 with link time reference"`, please run the following command to resolve it: 

```bash
export LD_LIBRARY_PATH=/usr/lib64/:${LD_LIBRARY_PATH}
```

## Error 3. Mamba not found in the shell

When running the following command:
```bash
$ python ROADIES/run_roadies.py --cores 1
```
You may encounter this error:

```bash
rm: cannot remove '<OUT_DIR>': No such file or directory
Unlocking working directory.
snakemake --cores 1 --config mode=accurate config_path=config/config.yaml num_threads=0 --use-conda --rerun-incomplete
Config file config/config.yaml is extended by additional config specified via the command line.
Building DAG of jobs...
CreateCondaEnvironmentException:
The 'mamba' command is not available in the shell /usr/bin/bash that will be used by Snakemake. You have to ensure that it is in your PATH, e.g., first activating the conda base environment with `conda activate base`.The mamba package manager (https://github.com/mamba-org/mamba) is a fast and robust conda replacement. It is the recommended way of using Snakemake's conda integration. It can be installed with `conda install -n base -c conda-forge mamba`. If you still prefer to use conda, you can enforce that by setting `--conda-frontend conda`.
```
This means `mamba` package manager is missing or not available in the environment.

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

## Error 4. Conda not recognized

This can happen if conda is not added to your system's PATH.

### Solution

To resolve this, please ensure conda is added to the PATH by running the following commands:

```bash
export PATH="$HOME/miniconda3/bin:$PATH"
source ~/.bashrc
```

## Error 5. Handling dependencies (glibc)

### Solution

Ensure that the glibc version on your system is updated to 2.29 or higher. Update your system libraries if necessary. Otherwise you may encounter this error:

```bash
workflow/scripts/lastz_32: /lib64/libm.so.6: version 'GLIBC_2.29' not found
```

## Error 6. PASTA fails with insufficient core count

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

## Error 7. MLIPPER fails with `undefined symbol: ATL_dGetNB` (GPU placement mode)

You may see the following when running `--mode placement --gpu`:

```bash
MLIPPER/MLIPPER: symbol lookup error: /lib/x86_64-linux-gnu/liblapack.so.3: undefined symbol: ATL_dGetNB
```

This means your system's `liblapack.so.3`/`libblas.so.3` (via Debian/Ubuntu's `update-alternatives`) currently resolves to an ATLAS build that's missing its own `libatlas.so.3` dependency - a broken/incomplete ATLAS package install, unrelated to MLIPPER or ROADIES itself.

### Solution

Check which LAPACK/BLAS variant is active:

```bash
update-alternatives --display liblapack.so.3-x86_64-linux-gnu
update-alternatives --display libblas.so.3-x86_64-linux-gnu
```

If a non-ATLAS alternative is listed (commonly under `/usr/lib/x86_64-linux-gnu/lapack/` and `/usr/lib/x86_64-linux-gnu/blas/`), you can point the dynamic linker at it for your ROADIES session without changing the system-wide default:

```bash
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/lapack:/usr/lib/x86_64-linux-gnu/blas:$LD_LIBRARY_PATH"
```

Then re-run `run_roadies.py`. If no non-ATLAS alternative exists on your system, reinstalling `libatlas3-base` (or your distro's equivalent) should restore the missing `libatlas.so.3` dependency.