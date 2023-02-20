#usage: python3 reroot.py [reference tree] [tree to be rerooted] [output]
import sys
import numpy as np
from ete3 import Tree

#by ETE
def rerootTree(refTr, rertTr):
    refSub = refTr.children[0]
    refSubTips = set()
    for l in refSub:
        refSubTips.add(l.name)

    maxScore = 0
    for n in rertTr.traverse('preorder'):
        cladeTips = set()
        for l in n:
            cladeTips.add(l.name)

        common = refSubTips.intersection(cladeTips)
        score = len(refTr) - len(refSubTips) - len(cladeTips) + 2*len(common)
       
        if score > maxScore:
            newRoot = n
            maxScore = score

    rertTr.set_outgroup(newRoot)



# refTree = read_tree_newick(sys.argv[1])
# rerootedTree = read_tree_newick(sys.argv[2])
# reroot_TS(refTree, rerootedTree)
# rerootedTree.write_tree_newick(sys.argv[3])

# refTree = Phylo.read(sys.argv[1], 'newick')
# rerootedTree = Phylo.read(sys.argv[2], 'newick')
# reroot_BP(refTree, rerootedTree)
# Phylo.write(rerootedTree, sys.argv[3], 'newick')
if __name__=="__main__":
    refTree = Tree(sys.argv[1])
    rerootedTree = Tree(sys.argv[2])
    rerootTree(refTree, rerootedTree)
    print(rerootedTree)
    rerootedTree.write(outfile=sys.argv[3])