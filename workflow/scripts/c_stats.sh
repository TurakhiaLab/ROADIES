#USAGE: ./c_stats.sh [convergence dir] [outputdir] [length] [k] [ref newick]
#make sure that directories do not end with '/'
rm -r $2
mkdir $2
arr=()
touch $2/refdist.txt
for d in $1/*.nwk; do
    echo "$d"
    echo "Rscript --vanilla workflow/scripts/dist.R ${5} ${d} $2/refdist.txt"
    Rscript --vanilla workflow/scripts/dist.R ${5} ${d} $2/refdist.txt
    arr+=($d)
done 
echo ${arr[@]}
ITER=0
len=$(expr ${#arr[@]} - 1)
echo "${len}"
for (( i=0; i<$len;i++));
do 
    ITER2=$(expr $i + 1)
    echo "comparing ${arr[$i]} to ${arr[${ITER2}]}" 
    Rscript --vanilla workflow/scripts/dist.R ${arr[${i}]} ${arr[${ITER2}]} $2/iterdist.txt
done
python3 workflow/scripts/c_stats.py $2 $3 $4
