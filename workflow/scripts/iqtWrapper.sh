num=`grep -n '>' $1 | wc -l`
if [[ ${num} -gt $4 ]]
then
	iqtree -s $1 -redo
	if [[ $? -eq 0 ]]
	then
        mv "${1}.treefile" $2
        mv "${1}."* "${3}"
    else
		touch $2
	fi
else
	touch $2
fi
