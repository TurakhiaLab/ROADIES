touch $1
for fileName in $@
do
    cat $fileName >> $1
done