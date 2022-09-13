num=`grep -n '>' $1 | wc -l`
if [[ ${num} > 2 ]]
then
	iqtree -s $1
else
	cp $1 $2
fi
