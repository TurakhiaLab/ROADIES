#get_run.sh copies necessary ROADIES output files to the corresponding converge directory
#USAGE: ./workflow/scripts/get_run.sh [converge directory] [run id]
mkdir -p $1/$2
latest_log=`ls -Art .snakemake/log/ | tail -n 1`
cp .snakemake/log/${latest_log} $1/$2/$2.log
cp -r results $1/$2
cp results/genetrees/gene_tree_merged.nwk $1/$2
cp results/genes/mapping.txt $1/$2
cp freqQuad.csv $1/$2
cp -r results/benchmarks $1/$2
rm -r ~/.pasta