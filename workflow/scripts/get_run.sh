#get_run.sh copies necessary ROADIES output files to the corresponding converge directory
#USAGE: ./workflow/scripts/get_run.sh [converge directory] [run id]
mkdir -p $1/$2
latest_log=`ls -Art .snakemake/log/ | tail -n 1`
echo "Latest log is .snakemake/log/${latest_log} adding to $1/$2/$2.log"
python workflow/scripts/logparser.py ${latest_log} $1/$2
cp .snakemake/log/${latest_log} $1/$2/$2.log
echo "adding results to $1/$2"
cp results/roadies.nwk $1/$2
cp results/roadies_stats.nwk $1/$2
cp results/genetrees/gene_tree_merged.nwk $1/$2
cp results/genes/mapping.txt $1/$2
cp freqQuad.csv $1/$2
cp $3 $1/$2
cp $4 $1/$2
cp $5 $1/$2
rm -r ~/.pasta

