num=`grep -n '>' $1 | wc -l`
if [[ ${num} > 2 ]]
then
	mafft --auto $1 > $2
else
	cp $1 $2
fi
