from treeswift import *
from os import listdir
from os.path import isfile, getsize

# for i in range(1,301):
#     if isfile('test/geneTree/gene_tree_' + str(i) + '.newick') == False:
#         print(i)

gtList = [f for f in listdir('test/geneTree/') if f.find('.newick') != -1 and f.find('merged') == -1 and getsize('test/geneTree/' + f) != 0]
# print(len(gtList))

# tree = read_tree_newick('test/geneTree/gene_tree_4.newick')
# tree.draw()
for gtFile in gtList:
    # print(gtFile)
    flag = 0
    tree = read_tree_newick('test/geneTree/' + gtFile)
    for n in tree.traverse_preorder(leaves=False):
        if(len(n.child_nodes()) > 2):
            print(len(n.child_nodes()))
            flag = 1
            # break
    # if flag == 0:
    #     print('aaaaaaaaaa  ' + gtFile)