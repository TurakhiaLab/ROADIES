Use the following code to download NCBI bird assemblies:

To generate a ranking of NCBI bird genome assemblies by Contig N50 use:
`python Scrape.py' in order to generate a text file Datarank.txt that is a csv file with species name, N50, and link to assembly stats.

To download ALL genomes use:
`./get_fasta.sh [output directory]`

To download the top n genomes by contig N50 use:
`./get_fasta_small.sh [number of genomes] [output directory]
default value is 30`

