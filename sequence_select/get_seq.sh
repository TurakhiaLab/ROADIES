cat $1 | parallel --colsep=',' python my_select.py $2/{1} -k {2} -l 1000
cat *.fa >> $3
