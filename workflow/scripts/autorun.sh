mkdir -p runs
for L in 200 500 1000 1500
do
    for K in 200 400 600 800 1000 1200 1400 1600 2000
    do
        rm -r results
        mkdir results
        snakemake --core $1 --config LENGTH=$L KREG=$K --use-conda --rerun-incomplete
        ./workflow/scripts/get_run.sh runs run_L$L.K$K
        log=`ls .snakemake/log | tail -n 1`
        echo "$log"
        python3 workflow/scripts/logparser.py .snakemake/log/$log runs/run_L$L.K$K
        snakemake --report runs/run_L$L.K$K/report.html
    done
done