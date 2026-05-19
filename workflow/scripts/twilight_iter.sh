#!/bin/bash

seqFile=$1
threads=$2
tempDir=$3
output=$4

export OMP_NUM_THREADS=$threads
        
# Initial Guide Tree: MAFFT 
mkdir -p $tempDir
mafft --retree 0 --treeout --reorder --thread $threads $seqFile > mafft.out
python3 workflow/scripts/mafft2nwk.py  $seqFile.tree $seqFile $tempDir/tree_iter0.nwk
rm mafft.out && rm $seqFile.tree

# MSA Iter 1: TWILIGHT
# twilight -i $seqFile -t $tempDir/tree_iter0.nwk -o $tempDir/msa_iter1.fa -C $threads --no-filtering --check
/home/ang037@AD.UCSD.EDU/TWILIGHT/bin/twilight -i $seqFile -t $tempDir/tree_iter0.nwk -o $tempDir/msa_iter1.fa -C $threads --match 40 --mismatch -7 --transition 17 --gap-open -140 --gap-extend -10 --overwrite

cp $tempDir/msa_iter1.fa $output 