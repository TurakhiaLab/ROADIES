#!/bin/bash
mkdir alignments
filepath=`pwd`
cat ../sequence_select/index.csv | parallel --colsep=',' /home/AD.UCSD.EDU/ttl074/lastz-distrib/bin/lastz_32 $1/{1}[multiple] ../sequence_select/out.fasta[multiple] --notransition  --format=maf --noytrim --filter=coverage:70 --filter=identity:70 --output={1}.maf
mv *.maf alignments
