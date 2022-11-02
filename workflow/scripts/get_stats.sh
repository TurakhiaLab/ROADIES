#!/bin/bash
# usage ./get_stats.sh [input dir] [output dir] [reference tree]
mkdir -p $2
arr=()
for d in ${1}/*/ ; do
    echo "$d"
    Rscript --vanilla workflow/scripts/dist.R $3 ${d}speciesTree.newick $2/refdist.txt
    Rscript --vanilla workflow/scripts/dist2.R $3 ${d}speciesTree.newick $2/refdist2.txt

    arr+=($d)
done 
python workflow/scripts/get_stats.py $1 $2
echo ${arr[@]}
ITER=0
len=$(expr ${#arr[@]} - 1)
echo "${len}"
for (( i=0; i<$len;i++))
do 
    ITER2=$(expr $i + 1)
    echo "comparing ${arr[$i]} to ${arr[${ITER2}]}" 
    Rscript --vanilla workflow/scripts/dist.R ${arr[${i}]}speciesTree.newick ${arr[${ITER2}]}speciesTree.newick $2/iterdist.txt
    Rscript --vanilla workflow/scripts/dist2.R ${arr[${i}]}speciesTree.newick ${arr[${ITER2}]}speciesTree.newick $2/iterdist2.txt

done
python workflow/scripts/rf.py $2