#get_run.sh copies necessary ROADIES output files to the corresponding converge directory
#USAGE: ./workflow/scripts/get_run.sh [converge directory] [run id]
latest_log=`ls -Art .snakemake/log/ | tail -n 1`
snakemake --report report.html 
rm -r ~/.pasta