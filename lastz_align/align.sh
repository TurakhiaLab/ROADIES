#!/bin/bash
mkdir alignments2
filepath=`pwd`
echo "$1" | parallel /home/AD.UCSD.EDU/ttl074/lastz-distrib/bin/lastz_32 [multiple] ../Sequence_Select/out.fasta[multiple] --notransition --step=20  --format=maf --noytrim --filter=coverage:70 --filter=identity:70 >${FILE}.maf 
cd $1
mv ./*.maf ${filepath}/alignments2
cd ${filepath}
