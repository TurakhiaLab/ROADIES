#!/bin/bash
python get_fasta_small.py $1
input="fasta_links_small.txt"
mkdir fasta_genomes_small
while IFS= read -r line
do
	wget ${line} -P $2
done < "$input"
