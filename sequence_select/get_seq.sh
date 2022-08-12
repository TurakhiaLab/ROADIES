cat $1 | parallel --colsep=',' python new_select.py $2/{1} -k {2} -l $4 
cat *.fa >> $3
