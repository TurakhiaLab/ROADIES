# wga-phylo

SnakeMake workflow for sequence selection to lastz

`snakemake --core [number of cores] --config PATH=["path directory containing all genomic fasta files"] LENGTH=[Length of region] KREG=[Number of regions] OUT=["output directory"]`

If you want to rerun please use following code to get rid of previous files

`snakemake clean --core [number of cores] --config PATH=["path directory containing all genomic fasta files"] LENGTH=[Length of region] KREG=[Number of regions] OUT=["output directory"]`
