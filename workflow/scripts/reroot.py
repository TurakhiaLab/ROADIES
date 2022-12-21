#usage: python3 reroot.py [reference tree] [tree to be rerooted] [output]
import sys
import numpy as np
from scipy.linalg import expm
from treeswift import *
from Bio import Phylo
from ete3 import Tree

#by tree swift
#not used
def reroot_TS(refTr, rertTr):
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

#by Bio.Phylo
#notused
def reroot_BP(refTr, rertTr):
    # Phylo.draw_ascii(rertTr)
    terminals = []
    for t in refTr.get_terminals(order='preorder'):
        terminals.append(t)

    refSub = refTr.get_path(terminals[0])[0]
    refSubTips = set()
    for t in refSub.get_terminals():
        refSubTips.add(str(t))

    maxScore = 0
    for c in rertTr.find_clades():
        cladeTips = set()
        for t in c.get_terminals():
            cladeTips.add(str(t))

        common = refSubTips.intersection(cladeTips)
        score = len(terminals) - len(refSubTips) - len(cladeTips) + 2*len(common)
        
        # tips=[]
        # if 'Ostrich' in cladeTips and 'Great_Tinamou' in cladeTips:
        #     for t in cladeTips:
        #         print(t)
        #     print()
        #     for t in refSubTips:
        #         print(t)
        #     print()
        #     for t in common:
        #         print(t)

        #     print(len(terminals))
        #     print(len(refSubTips))
        #     print(len(cladeTips))
        #     print(len(common))
        #     print(score)
        #     print('------')

        if score > maxScore:
            newRoot = c
            maxScore = score

    rertTr.root_with_outgroup(newRoot)

    # Phylo.draw_ascii(rertTr)

    # terminals = []
    # for t in rertTr.get_terminals(order='preorder'):
    #     terminals.append(t)
    # rtsub = rertTr.get_path(terminals[-2])[0]
    # for t in rtsub.get_terminals():
    #     print(t)

#by ETE
def reroot(refTr, rertTr):
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

refTree = Tree(sys.argv[1])
rerootedTree = Tree(sys.argv[2])
reroot(refTree, rerootedTree)
print(rerootedTree)
rerootedTree.write(outfile=sys.argv[3])