#USAGE: ./c_stats.sh [convergence dir] [outputdir] [length] [ref newick]
#make sure that directories do not end with '/'
rm -r $2
mkdir $2
arr=()
for d in ${1}/*.nwk ; do
    echo "$d"
    Rscript --vanilla dist.R $4 ${d} $2/refdist.txt
    Rscript --vanilla dist2.R $4 ${d} $2/refdist2.txt
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
    Rscript --vanilla dist.R ${arr[${i}]} ${arr[${ITER2}]} $2/iterdist.txt
    echo "alternative"
    Rscript --vanilla dist2.R ${arr[${i}]} ${arr[${ITER2}]} $2/iterdist2.txt
done
python3 c_stats.py $2 $3
