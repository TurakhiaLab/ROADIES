#!/bin/bash
rm -r $2
mkdir $2
arr=()
for d in ${1}/*/ ; do
    echo "$d"
    Rscript --vanilla dist.R $3 ${d}speciesTree.newick $2/refdist.txt
    Rscript --vanilla dist2.R $3 ${d}speciesTree.newick $2/refdist2.txt

    arr+=($d)
done 
python get_stats.py $1 $2
echo ${arr[@]}
ITER=0
len=$(expr ${#arr[@]} - 1)
echo "${len}"
for (( i=0; i<$len;i++))
do 
    ITER2=$(expr $i + 1)
    echo "comparing ${arr[$i]} to ${arr[${ITER2}]}" 
    Rscript --vanilla dist.R ${arr[${i}]}speciesTree.newick ${arr[${ITER2}]}speciesTree.newick $2/iterdist.txt
    Rscript --vanilla dist2.R ${arr[${i}]}speciesTree.newick ${arr[${ITER2}]}speciesTree.newick $2/iterdist2.txt

done
python rf.py $2