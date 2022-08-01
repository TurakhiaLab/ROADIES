#!/bin/bash
input="fasta_links.txt"
mkdir fasta_genomes
while IFS= read -r line
do
	wget ${line} -P ./fasta_genomes
done < "$input"
