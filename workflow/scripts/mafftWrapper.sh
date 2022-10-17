num=`grep -n '>' $1 | wc -l`
if [[ ${num} -gt $4 ]]
then
	mafft --auto $1 > $2
else
	touch $2
fi
