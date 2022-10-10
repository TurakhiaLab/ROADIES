for L in 200 500 1000
do
    for n in 1 2 3
    do
        rm -r test/*
        snakemake -c16 --config LENGTH=$L
        ./workflow/scripts/get_run.sh run_L$L.$n
        snakemake --report runs/run_L$L.$n/report.html
    done
done