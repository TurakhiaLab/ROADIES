#get_run.sh copies necessary ROADIES output files to the corresponding converge directory
#USAGE: ./workflow/scripts/get_run.sh [converge directory] [run id]
mkdir $1/$2
latest_log=`ls -Art .snakemake/log/ | tail -n 1`
echo "Latest log is ${latest_log}"
python workflow/scripts/logparser.py ${latest_log} $1/$2
cp -a results/plots $1/$2
cp -a results/statistics $1/$2
cp results/roadies.nwk $1/$2
cp results/geneTree/gene_tree_merged.nwk $1/$2
cp results/genes/mapping.txt $1/$2
