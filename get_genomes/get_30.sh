#!/bin/bash
input="fasta_links_30.txt"
mkdir fasta_genomes_30
while IFS= read -r line
do
	wget ${line} -P ./fasta_genomes
done < "$input"
