# Detailed Usage

This section provides detailed instructions on how to configure the ROADIES pipeline further for various user requirements with your own genomic dataset. Once the required environment setup process is complete, follow the steps below.

## Step 1: Specify input genomic dataset

After installing the environment, you need to get input genomic sequences for creating the species tree. To run ROADIES with your own dataset, update the `config.yaml` file (found in the ROADIES directory - `config` folder) to include the path to your input datasets under the `GENOMES` parameter.

!!! Note 
    All input genome assemblies in the path mentioned in `GENOMES` should be in `.fa` or `.fa.gz` format. The genome assembly files should be named according to the species' names (for example, Aardvark's genome assembly is to be named `Aardvark.fa`). Each file should contain the genome assembly of one unique species. If a file contains multiple species, split it into individual genome files (fasplit can be used for this: `faSplit byname <input_dir> <output_dir>`). Moreover, the file name should not have any special characters like `.` (apart from `_`) - for example, if the file name is `Aardvark.1.fa`, rename it to `Aardvark_1.fa`.

## Step 2: Modify Other Configuration Paramters

Adjust other parameters listed in `config.yaml` as per specific user requirements. Details of the parameters are mentioned below.
 
!!! Note
    ROADIES has default values for some of the parameters that give the best results and are recommended in general. However, users can optionally modify the values specific to their needs.

| Parameters | Description | Default value |
| --- | --- | --- |
| **GENOMES** | Specify the path to your input files which includes raw genome assemblies of the species. | |
| **REFERENCE** (optional) | Specify the path for the reference tree (state-of-the-art) in Newick format to compare ROADIES' results with a state-of-the-art approach. If you don't want to specify any reference tree, set it to `NULL`. | `NULL` |
| **LENGTH** | Configure the lengths of each of the randomly sampled subsequences or genes. | 500 |
| **GENE_COUNT** | Configure the number of genes to be sampled across all input genome assemblies. In normal mode, this will be the count of the genes to be sampled. In `--converge` mode, this will be the initial count of the number of genes for the first iteration and this value will be doubled iteratively. | 250 |
| **UPPER_CASE** | Configure the lower limit threshold of upper cases for valid sampling. ROADIES samples the genes only if the percentage of upper cases in each gene is more than this value. | 0.9 (Recommended) |
| **OUT_DIR** | Specify the path for ROADIES output files (this saves the current iteration results in converge mode). | |
| **ALL_OUT_DIR** | Specify the path for ROADIES output files for all iterations in converge mode. | |
| **MIN_ALIGN (deprecated)** | Specify the minimum number of allowed species to exist in gene fasta files after LASTZ. This parameter is used for filtering gene fasta files which has very less species representation. It is recommended to set the value greater than or equal to 4 since ASTRAL-Pro follows a quartet-based topology for species tree inference. For larger evolutionary timescales, we recommended setting it to a much higher value. In such cases, 15 to 20 would be a good start. | 4 |
| **COVERAGE** | Set the percentage of input sequence included in the alignment for LASTZ. | 85 |
| **CONTINUITY** | Define the allowable percentage of non-gappy alignment columns for LASTZ. | 85 |
| **IDENTITY** | Set the percentage of the aligned base pairs (matches/mismatches) for LASTZ. For larger evolutionary timescales, consider lowering the identity values than default for more homologous hits to be encountered. | 65 | 
| **IDENTITY_DEEP** | Set the percentage of the aligned base pairs (matches/mismatches) for LASTZ, for larger evolutionary timescales (needed when `--deep True` - details below). | 40 | 
| **MAX_DUP** | Specify maximum number of allowed gene copies from one input genome in an alignment. | 10 |
| **STEPS** |Specify the number of steps in the LASTZ sampling (increasing number speeds up alignment but decreases LASTZ accuracy).|1 |
| **FILTERFRAGMENTS** | Specify the portion so that sites with less than the specified portion of non-gap characters in PASTA alignments will be masked out. If it is set to 0.5, then sites with less than 50% of non-gap characters will be masked out. | 0.5 |
| **MASKSITES** | Specify the portion so that sequences with less than the specified portion of non-gap sequences will be removed in PASTA alignment. If it is set to 0.05, then sequences having less than 5% of non-gap characters (i.e., more than 95% gaps) will be masked out.| 0.02 |
| **SUPPORT_THRESHOLD** | Specify the threshold so that support values with equal to or higher than this threshold is considered as highly supported node. Such highly supported nodes crossing this threshold will be counted at every iteration to check the confidence of the tree (works in `--converge` mode). | 0.95 |
| **NUM_INSTANCES** | Specify the number of instances for PASTA, LASTZ, MashTree and RAxML-NG to run in parallel. It is recommended to set the number of instances equal to (`--cores`/4) for optimal runtime. | 4 | 
| **SCORES** | Set the alignment scores for LASTZ (needed when `--deep True` - details below). | `HOXD55.q` - file is provided along with the ROADIES package |

## Step 3: Run the ROADIES pipeline

Once the required installations are completed and the parameters are configured in `config.yaml` file, execute the following command (from ROADIES repo home directory):

```bash
python run_roadies.py --cores <number of cores>
```

This will let ROADIES run in accurate mode by default with specified number of cores. After the completion of the execution, the output species tree in Newick format will be saved as `roadies.nwk` in a separate `output_files` folder.

## Command line arguments

There are multiple command line arguments through which user can change the mode of operation, specify the custom config file path, etc.

| Argument | Description |
| --- | --- |
| `--cores` | Specify the number of cores |
| `--mode` | Specify [modes of operation](index.md#modes-of-operation) (`accurate`, `balanced` or `fast`).`accurate` mode is the default mode. | 
| `--noconverge` | Run ROADIES in non converge mode (for single iteration) if you know the optimal gene count to start with |
| `--config` | Provide optional custom YAML files (in the same format as `config.yaml` provided with this repository). If not given, by default `config/config.yaml` file will be considered.|
| `--deep` | Specify if ROADIES will evaluate deeper phylogeny. Set it to `True` or `False`. By default, its set to `False`. |

For example:

```
python run_roadies.py --cores 16 --mode balanced --noconverge --config config/config.yaml --deep True
```

Use `--help` to get the list of command line arguments.

## Step 4: Analyze output files

### Output from current iteration

After the pipeline finishes running, the final species tree of current iteration estimated by ROADIES will be saved as `roadies.nwk` inside a separate folder mentioned in the `--OUT_DIR` parameter in the `config/config.yaml` file. 

ROADIES also provides a number of intermediate output files for extensive debugging by the user, described below:

1. `alignments` - this folder contains the LASTZ alignment output of all individual input genomes aligned with randomly sampled gene sequences.
2. `benchmarks` - this folder contains the runtime value of each of the individual jobs for each of the stages in the pipeline. These files will only be used if you want to estimate and compare the stagewise runtime of various pipeline stages and will not be used in final tree estimation. 
3. `genes` - this folder contains the output files of multiple sequence alignment and tree-building stages (run by PASTA, IQTREE/FastTree, MashTree) of the pipeline. 
4. `genetrees` - this folder contains two files as follows:
    - `gene_tree_merged.nwk` - this file lists all gene trees together generated by IQTREE/FastTree/MashTree. It is used by ASTRAL-Pro to estimate the final species tree from this list of gene trees.
    - `original_list.txt` - this file lists all gene trees together corresponding to their gene IDs. Some lines will have only gene IDs but no associated gene trees. This is because some genes will be filtered out from tree building and MSA step if it has less than four species. Hence this file also lists those gene IDs with missing gene trees for further debugging. 
5. `plots` - this folder contains four following plots:
    - `gene_dup.png` - this histogram plot represents the count of the number of gene duplicates on the Y-axis vs. the number of genes having duplication on the X-axis.
    - `homologues.png` - this histogram plot represents the count of the number of genes on the Y-axis vs. the number of homologous species on the X-axis.
    - `num_genes.png` - this plot represents how many genes out of `--GENE_COUNT` parameter have been aligned to each of the input genomes after the LASTZ step. The X-axis represents different genomes, and the Y-axis represents the number of genes.
    - `sampling.png` - the plot shows how many genes have been sampled from each of the input genomes after the random sampling step. The X-axis represents different genomes, and the Y-axis represents the number of genes.
6. `samples` - this folder contains the list of randomly sampled genes from individual input genomes. 
    - `<species_name>_temp.fa` - these files contain genes sampled from the particular input genome.
    - `out.fa` - this file contains all sampled subsequences (genes) from individual genomes combined, which is given to the the LASTZ step. 
7. `statistics` - this folder contains CSV data for the plots shown in the `plots` directory mentioned above.
    - `gene_to_species.csv` - this is an additional CSV file (corresponding plots to be added in future) which provides the information about which genes are aligned to what species after LASTZ step (`num_genes.csv` only gives the total count of the genes per species, `gene_to_species.csv` also gives the ID number of those aligned genes). Along with each gene ID number, it also provides the [score, line number in .maf file, position] of all the homologs of that particular gene. Score, position and line number information is collected from the corresponding species' .maf file (generated by LASTZ), saved in `results/alignments` folder.
8. `roadies_stats.nwk`- this is the final estimated species tree (same as `roadies.nwk`), along with the support branch values in the Newick tree. 
9. `roadies.nwk`- this is the final estimated species tree in Newick format.
10. `roadies_rerooted.nwk` (optional) - this is the final estimated species tree, re-rooted corresponding to the outgroup node from the given reference tree (provided as `REFERENCE` in `config.yaml`).
11. `time_stamps.csv` - this file contains the start time, number of gene trees required for estimating species tree, end time, and total runtime (in seconds), respectively.
12. `ref_dist.csv` - this file provides the number of gene trees and the Normalized Robinson-Foulds distance between the final estimated species tree (i.e., `roadies.nwk`) and the reference tree (i.e., REFERENCE parameter in `config.yaml`).

### Extra output files from all iterations (these are not generated if --noconverge is used)

By default, results of all iterations (along with the corresponding species tree in the name `iteration_<iteration_number>.nwk`) will be saved in a separate folder mentioned in the `--ALL_OUT_DIR` parameter in the `config/config.yaml` file.

!!! Note
    With `--noconverge` option, ROADIES only saves the results of the current ongoing iteration in the folder specified by `--OUT_DIR` and the files below won't be generated.

For extensive debugging, other intermediate output files for each stage of the pipeline for each iterations are saved in `--ALL_OUT_DIR` as follows:

1. Folder with `iteration_<iteration_number>` - this folder contains results from the specific iteration corresponding to the iteration number in the folder name.
    -  Folder with name in `--OUT_DIR` - this contains the results of all stages of the pipeline (as described above in non convergence section). 
    - `gene_tree_merged.nwk` - this file lists all gene trees together generated by IQTREE/FastTree/MashTree in that particular iteration. It is concatenated with master list of gene trees from all past iterations before providing to ASTRAL-Pro to estimate the final converged species tree.
    - `iteration_<iteration_number>.log` - this file contains the log information of the corresponding iteration execution. 
    - `mapping.txt` - This file maps all gene names in the gene trees with the corresponding species name from where it originates. It is required by ASTRAL-Pro, along with the master list of gene trees from all iterations, to infer species tree. 
2. `iteration_<iteration_number>_stats.nwk` - this is the final estimated species tree for the corresponding iteration (same as `iteration_<iteration_number>.nwk`), along with the support branch values in the Newick tree. 
3. `iteration_<iteration_number>.nwk` - this is the final estimated species tree for the corresponding iteration
4. `iteration_<iteration_number>.rerooted.nwk` - (optional) - this is the final estimated species tree for the corresponding iteration, re-rooted to the outgroup node from the given reference tree (provided as `REFERENCE` in `config.yaml`).
5. `master_gt.nwk` - this is the concatenated list of all gene trees from all iterations together.
6. `master_map.txt` - this is the concatenated list of all mapping files from all iterations together. This `master_gt.nwk` and `master_map.txt` is provided to ASTRAL-Pro after every iteration to get the converged species tree. 
7. `ref_dist.csv` - this file provides the iteration number, number of gene trees and the Normalized Robinson-Foulds distance between the final estimated species tree (i.e., `roadies.nwk`) and the reference tree (i.e., REFERENCE parameter in `config.yaml`), for all iterations.
8. `time_stamps.csv`- this file contains the start time in first line, iteration number, number of gene trees required for estimating species tree, end time, and total runtime (in seconds), respectively, for all iterations in subsequent lines.

# Run ROADIES in a multi-node cluster (using SLURM) (currently being tested)

To run ROADIES in a multi-node cluster, make the following changes in the file `workflow/scripts/converge.py` (for `--noconverge` mode - make changes in `workflow/scripts/noconverge.py`)

Replace below lines:

```
    cmd = [
        "snakemake",
         "--cores",
         str(cores),
         "--config",
         "mode=" + str(mode),
         "config_path=" + str(config_path),
         "num_threads=" + str(num_threads),
         "deep_mode=" + str(deep_mode),
         "MIN_ALIGN=" + str(MIN_ALIGN),
         "--use-conda",
         "--rerun-incomplete",
    ]
```

With below lines (you can change the value of `--jobs` and other account details based on your cluster configuration):

```
    cmd = [
    "snakemake",
    "--jobs",
    "4",
    "--groups",
    "lastz=group0",
    "--group-components",
    "group0=8",
    "--config",
    "mode=" + str(mode),
    "config_path=" + str(config_path),
    "num_threads=" + str(num_threads),
    "--use-conda",
    "--rerun-incomplete",
    "--cluster",
    (
        "sbatch "
        "--job-name=XXX "
        "--partition=XXX "
        "--account=XXX "
        "--nodes=1 "
        "--ntasks-per-node=4 "
        "--cpus-per-task=8 "
        "--time=8-0 "
        "--mem-per-cpu=11G "
        "--output=%x_%j.out "
        "--error=%x_%j.err "
        "--mail-user=XXX "
        "--mail-type=ALL"
    )
]
```

After the above changes, save the following lines of code as separate file called `roadies.slurm` and run `sbatch roadies.slurm`.
```
#! /bin/bash
#SBATCH -J ROADIES_XXX
#SBATCH -p XXX
#SBATCH --account=XXX
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=8-0
#SBATCH --mem-per-cpu=11G
#SBATCH -o %x_%j.out 
#SBATCH -e %x_%j.err
#SBATCH --mail-user=XXX
#SBATCH --mail-type=ALL

echo Starting at `date`
echo This is job $SLURM_JOB_ID
echo Running on `hostname`

source /<PATH>/miniconda3/etc/profile.d/conda.sh
conda activate myenv
cd /<PATH>/miniconda3/envs/myenv/ROADIES
srun --nodes=1 python run_roadies.py --cores 128

echo Exiting at `date`
srun sleep 30
```

