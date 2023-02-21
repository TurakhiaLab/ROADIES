#reroot.py takes in a reference and a input tree and reroots the input with respect to the reference
#REQUIREMENTS: numpy and ete3
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

if __name__=="__main__":
    refTree = Tree(sys.argv[1])
    rerootedTree = Tree(sys.argv[2])
    rerootTree(refTree, rerootedTree)
    print(rerootedTree)
    rerootedTree.write(outfile=sys.argv[3])
