#!/bin/bash
python get_seq.py
input="fasta_links.txt"
mkdir fasta_genomes
while IFS= read -r line
do
	wget ${line} -P $1
done < "$input"
