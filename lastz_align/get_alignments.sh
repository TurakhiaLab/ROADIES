#!/bin/bash
mkdir test
filepath=`pwd`
cat ../Sequence_Select/index.csv | parallel --colsep=',' /home/AD.UCSD.EDU/ttl074/lastz-distrib/bin/lastz_32 $1/{1}[multiple] ../Sequence_Select/out.fasta[multiple] --notransition --step=20  --format=maf --noytrim --filter=coverage:70 --filter=identity:70 >> alignments.maf 
