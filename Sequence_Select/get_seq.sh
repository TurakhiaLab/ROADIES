cat index.csv | parallel --colsep=',' python my_select.py $1/{1} -k {2} -l 1000
cat *.fa >> out.fasta
