# Quick Start

This section provides quick steps to get acquainted with ROADIES. 

For more details about various possible ways of installations, refer to [Installations](setup.md) tab. For more details about the usage and further configurability, refer to [Usage](usage.md) tab.

## Quick Install (using docker)

First clone the repository

```
git clone https://github.com/TurakhiaLab/wga-phylo.git
cd wga-phylo
```

Then build and run docker

```
docker build roadies_image .
docker run -it roadies_image
```

## Get input genomic data
After installing the environment, we need to get input genomic sequences for creating the species tree. To start with this, we have provided few test genomes, present in the repository in `test/test_data` folder,

OR, download few genomes by executing the following command:

## Configure the config file

To run ROADIES with test data, modify the path for `GENOMES` in `config/config.yaml` as `test/test_data`.

To run ROADIES with downloaded genomes using `wget` commands mentioned above, provide the path of the downloaded genomes to `GENOMES` argument.

## Running the pipeline 

After modifying the config file, run the following command:

```
python run_roadies.py --cores 32
```

After the completion of the execution, the output species tree in Newick format will be saved as `roadies.nwk` in a separate `output_files` folder.

To run ROADIES in various other modes of operation (fast, balanced, accurate) (description of these modes are mentioned in [Home](index.md#modes-of-operation)), try the following commands:


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

For each modes, the output species tree will be saved as `roadies.nwk` in a separate `output_files` folder.

#### With convergence

To run ROADIES with converge mode (details mentioned in [Home](index.md#convergence-mechanism)), specify the input file in the same way mentioned above, then specify the number of iterations you want to run the convergence of ROADIES with `ITERATIONS` parameter in `config.yaml`. Then, run the following command (notice the addition of --converge argument).

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
The output files for all iterations will be saved in a separate `converge_files` folder. `output_files` will save the results of the last iteration. 

Species tree for all iterations will be saved in `converge_files` folder with the nomenclature `iteration_<iteration_number>.nwk`.