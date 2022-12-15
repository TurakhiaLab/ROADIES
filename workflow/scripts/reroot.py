#usage: python3 reroot.py [reference tree] [tree to be rerooted] [output]
import sys
import numpy as np
from scipy.linalg import expm
from treeswift import *

def reroot(refTr, rertTr):
    #find root
    for n in refTr.traverse_preorder():
        if n.is_root():
            root = n

    #find bipartitions by root in the reference tree
    rootChildren = root.child_nodes()
    refSubtree = refTr.extract_subtree(rootChildren[0])
    refStLeaves = [] # leaves of ref subtree
    leaves = [] # all leaves
    for n in refTr.traverse_leaves():
        leaves.append(n.get_label())
    for n in refSubtree.traverse_leaves():
        refStLeaves.append(n.get_label())

    #find the new root for the other tree
    maxScore = 0
    for n in rertTr.traverse_preorder():
        # print(n.get_label())
        subtree = rertTr.extract_subtree(n)
        stLeaves = []
        for l in subtree.traverse_leaves():
            stLeaves.append(l.get_label())
            # print(l.get_label())
        common = [l for l in stLeaves if l in refStLeaves]
        score = len(leaves) - len(stLeaves) - len(refStLeaves) + 2*len(common)
        # print(score)
        # print('-----')
        if score > maxScore:
            newRoot = n
            maxScore = score

    #rerooting
    rertTr.reroot(newRoot)

refTree = read_tree_newick(sys.argv[1])
rerootedTree = read_tree_newick(sys.argv[2])
reroot(refTree, rerootedTree)
rerootedTree.write_tree_newick(sys.argv[3])