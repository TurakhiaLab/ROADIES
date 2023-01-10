size=`stat --printf="%s" $1`
echo $size
if [[ ${size} -gt 3000000000 ]]
then
    echo $1
    echo "file greater than lastz limit"
    python3 workflow/scripts/get_names.py $1 $6 4
    echo "aligning $6/$7.0.subset"
    lastz_32 $1[subset=$6/$7.0.subset,multiple] $2[multiple] --filter=coverage:$3 --filter=identity:$4 --format=maf --output=$6/$7.0.maf --ambiguous=iupac
    echo "aligning $6/$7.1.subset"
    lastz_32 $1[subset=$6/$7.1.subset,multiple] $2[multiple] --filter=coverage:$3 --filter=identity:$4 --format=maf --output=$6/$7.1.maf --ambiguous=iupac
    cat $6/$7.0.maf >> $5
    cat $6/$7.1.maf >> $5
else
    echo $1
    echo "aligning normally"
    lastz_32 $1[multiple] $2[multiple] --filter=coverage:$3 --filter=identity:$4 --format=maf --output=$5 --ambiguous=iupac
fi