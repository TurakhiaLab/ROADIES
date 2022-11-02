#usage ./converge.sh [output dir] [cores] [number of iterations] [number of genes per iter] [statdir]
mkdir -p $1 
END=$3
mkdir -p $5
for L in 500 600 700
do
    mkdir -p $1/length_$L
    touch $1/length_$L/combined_gt.tre
    touch $1/length_$L/combined_mapping.txt
    for i in $(seq 1 $END);
    do
        if [ ${i} -lt 10 ];then
            i="0$i"
        fi
        echo $i
        rm -r results
        mkdir results
        snakemake --core $2 --config LENGTH=$L KREG=$4 --use-conda --rerun-incomplete
        ./workflow/scripts/get_run.sh $1/length_$L run_$i
        log=`ls .snakemake/log | tail -n 1`
        echo "$log"
        python3 workflow/scripts/logparser.py .snakemake/log/$log $1/length_$L/run_$i
        snakemake --report $1/length_$L/run_$i report.html
        cat $1/length_$L/run_$i/gene_tree_merged.newick >> $1/length_$L/combined_gt.tre
        cat $1/length_$L/run_$i/mapping.txt >> $1/length_$L/combined_mapping.txt
        ASTER-LINUX/bin/astral-pro -i {input} -o {output} -a {params.genes}/mapping.txt --polylimit {params.p}
    ./workflow/scripts/c_stats.sh $1/length_$L $5/length$_L $L $4
    done
done
mkdir -p $5/combined
python workflow/scripts/combine_plots.py $5 $5/combined converge