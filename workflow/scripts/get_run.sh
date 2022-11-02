
mkdir $1/$2
#mkdir $1/$2/plots
#mkdir $1/$2/statistics
cp -a results/plots $1/$2
cp -a results/statistics $1/$2
cp results/speciesTree.newick $1/$2
cp results/geneTree/gene_tree_merged.newick $1/$2
cp results/genes/mapping.txt $1/$2
