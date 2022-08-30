cat $1 | parallel --colsep=',' python my_select.py $2/{1} -k {2} -l $4 
cat *.fa >> $3
