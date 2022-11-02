#usage:: ./autorun.sh [outputdir] [cores] [statdir] [referencetree]
mkdir -p $1
for L in 200 500 1000 1500
do
    for K in 200 400 600 800 1000 1200 1400 1600 2000
    do
        rm -r results
        mkdir results
        mkdir -p $1/length_$L
        snakemake --core $2 --config LENGTH=$L KREG=$K --use-conda --rerun-incomplete
        ./workflow/scripts/get_run.sh $1/length_$L run_L$L.K$K
        log=`ls .snakemake/log | tail -n 1`
        echo "$log"
        python3 workflow/scripts/logparser.py .snakemake/log/$log $1/length_$L/run_L$L.K$K
        snakemake --report $1/length_$L/run_L$L.K$K/report.html
    done
    ./workflow/scripts/get_stats.sh $1/length_$L $3 $4
done
mkdir -p $3/combined
python workflow/scripts/combine_plots.py $4 $4/combined autorun