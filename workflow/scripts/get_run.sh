mkdir runs/$1
mkdir runs/$1/plots
mkdir runs/$1/statistics
cp -a results/plots runs/$1/plots
cp -a results/statistics runs/$1/statistics
cp results/speciesTree.newick runs/$1
cp results.speciesTree.newick runs/trees 
