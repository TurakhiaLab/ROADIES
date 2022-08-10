#!/bin/bash
home = `pwd`
python get_seq.py
input="fasta_links.txt"
mkdir fasta_genomes
while IFS= read -r line
do
	wget ${line} -P $1
done < "$input"
cd $1
gunzip *.gz
cd ${home}
