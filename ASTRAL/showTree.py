from treeswift import *
import sys

tree = read_tree_newick(sys.argv[1])
tree.draw(show_labels=True)