num=`grep -n '>' $1 | wc -l`
if [[ ${num} -gt $4 ]]
then
	mafft --localpair --maxiterate 1000 $1 > $2
else
	touch $2
fi
