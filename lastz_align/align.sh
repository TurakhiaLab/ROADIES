#!/bin/bash
mkdir alignments
filepath=`pwd`
for FILE in $1/*; do
	echo "aligning ${FILE} to sequences"
	/home/AD.UCSD.EDU/ttl074/lastz-distrib/bin/lastz_32 ${FILE}[multiple] ../Sequence_Select/out.fasta[multiple] --notransition --step=20  --format=maf --noytrim --filter=coverage:70 --filter=identity:70 >${FILE}.maf 
done
cd $1
mv ./*.maf ${filepath}/alignments
cd ${filepath}
