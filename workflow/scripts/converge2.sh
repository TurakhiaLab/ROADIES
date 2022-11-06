#usage ./converge.sh [output dir] [cores] [number of consecutive] [threshold] [number of genes per iter] [statdir] [max runs] [ref nwk]
#example: ./converge.sh converge 16 3 0.05 200 cstats 20 trees/cn48.nwk
mkdir -p $1 
mkdir -p $6
for L in 500 600 700
do
    mkdir -p $1/length_$L
    touch $1/length_$L/combined_gt.tre
    touch $1/length_$L/combined_mapping.txt
    counter=0
    i=1
    while [ $counter -le $3 ] && [ $i -le $7 ];
    do
        if [ ${i} -lt 10 ];then
            i="0$i"
        fi
        echo $i
        rm -r results
        mkdir results
        snakemake --core $2 --config LENGTH=$L KREG=$5 --use-conda --rerun-incomplete
        ./workflow/scripts/get_run.sh $1/length_$L run_$i
        log=`ls .snakemake/log | tail -n 1`
        echo "$log"
        python3 workflow/scripts/logparser.py .snakemake/log/$log $1/length_$L/run_$i
        snakemake --report $1/length_$L/run_$i report.html
        cat $1/length_$L/run_$i/gene_tree_merged.newick >> $1/length_$L/combined_gt.tre
        cat $1/length_$L/run_$i/mapping.txt >> $1/length_$L/combined_mapping.txt
        ASTER-Linux/bin/astral-pro -i $1/length_$L/combined_gt.tre -o $1/length_$L/run_$i.nwk -a $1/length_$L/combined_mapping.txt
        dist=($(Rscript --vanilla workflow/scripts/catdist.R $1 $2))
        echo $dist
        if [[ ${dist} -lt $4 ]]
        then
            counter=$((counter+1))
        else
            counter=0
        fi
        i=$((i+1))
        echo $i
        echo $counter
    ./workflow/scripts/c_stats.sh $1/length_$L $6/length$_L $L $8
    done
done
mkdir -p $6/combined
python workflow/scripts/combine_plots.py $6 $6/combined converge