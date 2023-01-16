num=`grep -n '>' $1 | wc -l`
if [[ ${num} -gt $4 ]]
then
	mafft --anysymbol --localpair --retree 2 --maxiterate 10 $1 > $2
	#mafft --auto $1 > $2
else
	touch $2
fi
