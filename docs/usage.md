# Usage

This section provides detailed instructions on how to configure the ROADIES pipeline further for various user requirements. Once the required environment setup process is complete, follow the steps below.

## Configuration paramters

For specific user requirements, ROADIES provides multiple parameters to be configured before running the pipeline. These input parameters are listed in `config/config.yaml`. 
 
!!! Note
    ROADIES has default values for some of the parameters that give the best results, users can optionally modify the values specific to their needs.

!!! Note 
    All input genome assemblies in the path mentioned in `GENOMES` should be in `.fa` or `.fa.gz` format. The genome assembly files should be named according to the species' names (for example, Aardvark's genome assembly is to be named `Aardvark.fa`). Each file should contain the genome assembly of one unique species. If a file contains multiple species, split it into individual genome files (fasplit can be used for this: `faSplit byname <input_dir> <output_dir>`)

| Parameters | Description | Default value |
| --- | --- | --- |
| **GENOMES** | Specify the path to your input files which includes raw genome assemblies of the species. | |
| **REFERENCE** (optional) | Specify the path for the reference tree (state-of-the-art) in Newick format to compare ROADIES' results with a state-of-the-art approach. If you don't want to specify any reference tree, set it to `null`. | `null` |
| **LENGTH** | Configure the lengths of each of the randomly sampled subsequences or genes. | 500 |
| **GENE_COUNT** | Configure the number of genes to be sampled across all input genome assemblies. | 750 |
| **UPPER_CASE** | Configure the lower limit threshold of upper cases for valid sampling. ROADIES samples the genes only if the percentage of upper cases in each gene is more than this value. | 0.9 (Recommended) |
| **OUT_DIR** | Specify the path for ROADIES output files (this saves the current iteration results in converge mode). | |
| **ALL_OUT_DIR** | Specify the path for ROADIES output files for all iterations in converge mode. | |
| **MIN_ALIGN** | Specify the minimum number of allowed species to exist in gene fasta files after LASTZ. This parameter is used for filtering gene fasta files which has very less species representation. It is recommended to set the value more than the default value since ASTRAL-Pro follows a quartet-based topology for species tree inference. | 4 (Recommended) |
| **COVERAGE** | Set the percentage of input sequence included in the alignment for LASTZ. | 85 (Recommended) |
| **CONTINUITY** | Define the allowable percentage of non-gappy alignment columns for LASTZ. | 85 (Recommended) |
| **IDENTITY** | Set the percentage of the aligned base pairs for LASTZ. | 65 (Recommended) | 
| **MAX_DUP** | Specify maximum number of allowed gene copies from one input genome in an alignment. | 10|
| **STEPS** |Specify the number of steps in the LASTZ sampling (increasing number speeds up alignment but decreases LASTZ accuracy).|1 (Recommended)|
| **FILTERFRAGMENTS** | Specify the portion so that sites with less than the specified portion of non-gap characters in PASTA alignments will be masked out. If it is set to 0.5, then sites with less than 50% of non-gap characters will be masked out. | 0.5 (Recommended)|
| **MASKSITES** | Specify the portion so that sequences with less than the specified portion of non-gap sequences will be removed in PASTA alignment. If it is set to 0.05, then sequences having less than 5% of non-gap characters (i.e., more than 95% gaps) will be masked out.| 0.02 (Recommended)|
| **SUPPORT_THRESHOLD** | Specify the threshold so that support values with equal to or higher than this threshold is considered as highly supported node. Such highly supported nodes crossing this threshold will be counted at every iteration to check the confidence of the tree (works in --converge mode). | 0.7 (Recommended) |
| **ITERATIONS** | Specify the number of iterations to be run in --converge mode. | | 


## Running the pipeline

Once the required installations are completed and the parameters are configured in `config/config.yaml` file, execute the following command:

```
python run_roadies.py --cores <number of cores>
```

This will let ROADIES run in accurate mode by default with specified number of cores in non converge mode. 

### Various command line arguments

There are multiple command line arguments through which user can change the mode of operation, specify the custom config file path, etc.

To modify the modes of operation, add the `--mode` command line argument as follows:

```
python run_roadies.py --cores <number of cores> --mode <`fast` OR `balanced` OR `accurate`>
```

This will run ROADIES with specified mode of operation in non converge mode. To run this in converge mode, add the `--converge` argument as follows:

```
python run_roadies.py --cores <number of cores> --mode <`fast` OR `balanced` OR `accurate`> --converge
```

Additionally, user can have custom YAML files (in the same format as config.yaml provided with this repository) with their own parameterizable values. Provide the custom YAML file using `--config` argument as follows (if not given, by default `config/config.yaml` file will be considered):

```
python run_roadies.py --cores <number of cores> --mode <`fast` OR `balanced` OR `accurate`> --config <add own config path>
```
Also, with the converge mode and the custom yaml file together, run the command as follows:

```
python run_roadies.py --cores <number of cores> --mode <`fast` OR `balanced` OR `accurate`> --config <add own config path> --converge
```

Use `--help` to get the list of command line arguments.


## Analyzing output files

### Without convergence

After the pipeline finishes running, the final species tree estimated by ROADIES will be saved as `roadies.nwk` inside a separate folder mentioned in the `--OUT_DIR` parameter in the `config/config.yaml` file. 

ROADIES also provides a number of intermediate output files for extensive debugging by the user. These files will be saved in `--OUT_DIR`, containing the following subfolders:

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

### With convergence

If converge option is enabled, the results of all iterations (along with the corresponding species tree in the name `iteration_<iteration_number>.nwk`) will be saved in a separate folder mentioned in the `--ALL_OUT_DIR` parameter in the `config/config.yaml` file.

!!! Note
    With `--converge` option, `--OUT_DIR` saves the results of the current ongoing iteration (if pipeline execution is finished, then the last iteration), whereas `--ALL_OUT_DIR` saves the results of all iterations executed. 

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