# reroot.py takes in a reference and a input tree and reroots the input with respect to the reference
# REQUIREMENTS: numpy and ete3
# usage: python3 reroot.py [reference tree] [tree to be rerooted] [output]

import sys
import numpy as np
from ete3 import Tree


def rerootTree(refTr, rertTr):
    refSub1 = refTr.children[0]
    refSubTips1 = set()
    for l in refSub1:
        refSubTips1.add(l.name)

    length1 = len(refSubTips1)

    refSub2 = refTr.children[1]
    refSubTips2 = set()
    for l in refSub2:
        refSubTips2.add(l.name)

    length2 = len(refSubTips2)
    refSubTips = set()

    if length1 > length2:
        refSubTips = refSubTips2
    else:
        refSubTips = refSubTips1

    maxScore = 0
    rertTr_temp = rertTr

    nodes = []
    for n in rertTr.traverse("preorder"):
        nodes.append(n)

    for n in nodes:
        if n == rertTr_temp:
            continue
        cladeTips = set()
        for l in n:
            cladeTips.add(l.name)

        common = refSubTips.intersection(cladeTips)
        score = len(refTr) - len(refSubTips) - len(cladeTips) + 2 * len(common)

        if score > maxScore:
            newRoot = n
            maxScore = score

    rertTr.set_outgroup(newRoot)


if __name__ == "__main__":
    refTree = Tree(sys.argv[1])
    rerootedTree = Tree(sys.argv[2])
    rerootTree(refTree, rerootedTree)
    rerootedTree.write(outfile=sys.argv[3])
