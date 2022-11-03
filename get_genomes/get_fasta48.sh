#!/bin/bash
home=`pwd`
input="links48.txt"
while IFS= read -r line
do
	wget ${line} -P $1
done < "$input"
cd $1
gunzip *.gz
cd ${home}
