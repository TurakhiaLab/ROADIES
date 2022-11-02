#USAGE: ./c_stats.sh [convergence dir] [outputdir] [length] [ref newick]
#make sure that directories do not end with '/'
mkdir $2
arr=()
for d in ${1}/*.nwk ; do
    echo "$d"
    Rscript --vanilla workflow/scripts/dist.R $5 ${d} $2/refdist.txt
    Rscript --vanilla workflow/scripts/dist2.R $5 ${d} $2/refdist2.txt
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
    echo "alternative"
    Rscript --vanilla workflow/scripts/dist2.R ${arr[${i}]} ${arr[${ITER2}]} $2/iterdist2.txt
done
python3 workflow/scripts/c_stats.py $2 $3 $4
