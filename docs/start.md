# Getting Started

This section provides quick steps to run ROADIES in your system once you have setup the environment, mentioned in [Installations](setup.md) page.

## Without convergence

After setting up the environment (explained in [setup](setup.md)), to run the ROADIES pipeline with 32 cores, specify the path to your input files as `GENOMES` in `config/config.yaml`. Then, run the following command:

```
python run_roadies.py --cores 32
```
!!! Note
    All input genome assemblies should be in `.fa` or `.fa.gz` format. The genome assembly files should be named according to the species' names (for example, Aardvark's genome assembly is to be named `Aardvark.fa`). Each file should contain the genome assembly of one unique species. If a file contains multiple species, split it into individual genome files (fasplit can be used for this: `faSplit byname <input_dir> <output_dir>`)

The output species tree will be saved as `roadies.nwk` in a separate folder provided by `--OUT_DIR` parameter.

To run ROADIES in various modes of operation (fast, balanced, accurate) (description of these modes are mentioned in [Home](index.md#modes-of-operation)), try the following commands:


```
python run_roadies.py --cores 32 --mode accurate
```

```
python run_roadies.py --cores 32 --mode balanced
```

```
python run_roadies.py --cores 32 --mode fast
```
!!! Note
    Accurate mode is the default mode of operation. If you don't specify any particular mode using `--mode` argument, default mode will run.

## With convergence

To run ROADIES with convergence (details mentioned in [Home](index.md#convergence-mechanism)) converge mode, specify the input file in the same way mentioned above, then specify the number of iterations you want to run the convergence of ROADIES with `ITERATIONS` parameter in `config.yaml`. Then, run the following command (notice the addition of --converge argument).

```
python run_roadies.py --cores 32 --converge
```
To try other modes, run as follows:

```
python run_roadies.py --cores 32 --mode balanced --converge
```
```
python run_roadies.py --cores 32 --mode fast --converge
```
The output files for all iterations will be saved in a separate folder specified in `--ALL_OUT_DIR`. Species tree for all iterations will be saved in this folder with the nomenclature `iteration_<iteration_number>.nwk`.