num=`grep -n '>' $1 | wc -l`
if [[ ${num} -gt 2 ]]
then
	FastTree -nt $1 > $2
	if [[ $? -ne 0 ]]
	then
		touch $2
	fi
else
	touch $2
fi
