#bootstrap.sh [gene tree file] [resampling size] [resampling number]

mkdir -p results/bootstrap
rm -f results/bootstrap/*

for i in $(seq 1 $3);
do
    python3 workflow/scripts/resampleLines.py $1 results/bootstrap/resampled_$i.newick $2
    ASTER-Linux/bin/astral-pro -i results/bootstrap/resampled_$i.newick -o results/bootstrap/speciesTree_$i.newick -a results/genes/mapping.txt
done
# Rscript workflow/scripts/rfcalc.r results/bootstrap speciesTree
Rscript workflow/scripts/rfcalc.r results/bootstrap speciesTree | grep mean | cut -d "=" -f 2 | cut -d " " -f 2

