# Quick start (with provided test dataset)

Once setup is done, you can run the ROADIES pipeline using the provided test dataset. Follow these steps for a 16-core machine:

**Step 1:** Go to ROADIES repository directory if not there:

```bash
cd ROADIES
```

**Step 2:** Create a directory for the test data and download the test datasets (using the following one line command):

```bash
mkdir -p test/test_data && cat test/input_genome_links.txt | xargs -I {} sh -c 'wget -O test/test_data/$(basename {}) {}'
```
**Step 3:** Run the pipeline with the following command (from ROADIES directory):

```bash
python run_roadies.py --cores 16
```

Step 2 will download the 11 Drosophila genomic datasets (links provided in `test/input_genome_links.txt`) and save them in the `test/test_data` directory. Step 3 will run ROADIES for those 11 Drosophila genomes and save the final newick tree as `roadies.nwk` in a separate `output_files` folder for the current iteration. The final output files for all iterations will be saved in `converge_files` folder upon completion. 

## Running ROADIES with different modes of operation

To run ROADIES in various other modes of operation (fast, balanced, accurate) (description of these modes are mentioned in [Modes of operation](index.md#modes-of-operation) section), try the following commands:

```bash
python run_roadies.py --cores 16 --mode accurate
```

```bash
python run_roadies.py --cores 16 --mode balanced
```

```bash
python run_roadies.py --cores 16 --mode fast
```
!!! Note
    Accurate mode is the default mode of operation. If you don't specify any particular mode using `--mode` argument, default mode will run.

For each modes, the output files for all iterations will be saved in a separate `converge_files` folder. `output_files` will save the results of the last iteration. Species tree for all iterations will be saved in `converge_files` folder with the nomenclature `iteration_<iteration_number>.nwk`.

## Running ROADIES in non converge mode (single iteration mode)

By default, ROADIES will run for multiple iteration until it gets a stable tree at the end (details mentioned in [convergence mechanism](index.md#convergence-mechanism) section). To run ROADIES with non converge mode (only for one iteration), execute the following command (notice the addition of `--noconverge` argument):

```bash
python run_roadies.py --cores 16 --noconverge
```

Try following commands for other modes:

```bash
python run_roadies.py --cores 16 --mode balanced --noconverge
```
```bash
python run_roadies.py --cores 16 --mode fast --noconverge
```