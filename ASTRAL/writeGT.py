from treeswift import *

gt = []
for filePath in snakemake.input:
    gt.append(read_tree_newick(filePath))

gtStr = ''
for t in gt:
    gtStr = gtStr + t.newick() + '\n'
gtFile = open(snakemake.output[0], 'w')
gtFile.write(gtStr)
gtFile.close()