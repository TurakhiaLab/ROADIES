# Quick start

After installing using one of the options mentioned in Quick Install, you're ready to run ROADIES! To get started:

## Step 1: Download the test dataset (11 Drosophila genomes):**

Make sure to run this step from the main ROADIES repository directory. 

```bash
mkdir -p test/test_data && cat test/input_genome_links.txt | xargs -I {} sh -c 'wget -O test/test_data/$(basename {}) {}'
```

This will save the datasets on a separate `test/test_data` folder within the repository

## Step 2: Run the pipeline

!!! Note
    ROADIES by default runs multiple iterations for generating highly accurate trees. For quick testing, use `--noconverge` to run a single iteration.

**Full run (multiple iterations)**
```bash
python run_roadies.py --cores 16
```
**OR**

**Quick test run (one iteration)**
```bash
python run_roadies.py --cores 16 --noconverge 
```

## Step 3: Analyze Output:

 - Final **UNROOTED** newick tree saved as `roadies.nwk` in the folder set by `OUT_DIR` in `config.yaml`, kept up to date after every iteration.
 - Intermediate files (if `--noconverge` not used) saved in the folder set by `ALL_OUT_DIR` in `config.yaml`.


!!! Note
    ROADIES outputs unrooted trees by default. You can reroot trees on your own or use the provided `reroot.py` script in `workflow/scripts/` (given a reference rooted species tree as input). 

!!! Note
    Already have a species tree and want to add more genomes to it (rather than rebuild from scratch), or want GPU acceleration? See the [ROADIES_XP guide](roadies_xp.md).

