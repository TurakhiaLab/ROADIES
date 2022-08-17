import numpy as np
from scipy.linalg import expm
from treeswift import *
import random

def random_seq(k, gtr_probs):
    '''
    This function generates a random sequence of length k using GTR stationary probabilities
    :param k: The length of the output sequence
    :param gtr_probs: The GTR stationary probabilities as a list [prob_A, prob_C, prob_G, prob_T]
    :return: A random string of length k generated using gtr_probs
    '''
    # prob_A, prob_C, prob_G, prob_T = gtr_probs # You can use these if it's more convenient
    arr = np.random.choice(['A', 'C', 'G', 'T'], size=k, p=gtr_probs)
    seq = ''.join(list(arr))
    return seq

def evolve(tree, root_seq, gtr_probs, gtr_rates):
    '''
    This function simulates the evolution of a root sequence down a given tree under the GTR model
    :param tree: A TreeSwift Tree object representing the phylogenetic tree
    :param root_seq: The root sequence to evolve down the tree
    :param gtr_probs: The GTR stationary probabilities as a list [prob_A, prob_C, prob_G, prob_T]
    :param gtr_rates: The GTR transition rates as a list [rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG]
    :return: A dictionary where keys are the labels of the given tree and values are evolved sequences
    '''
    prob_A, prob_C, prob_G, prob_T = gtr_probs # You can use these if it's more convenient
    rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG = gtr_rates # You can use these if it's more convenient
    initProb = {'A':prob_A, 'C':prob_C, 'G':prob_G, 'T':prob_T}

    #initialize R
    #-------------------------------------------------------
    R = {
        'A': {'A':0             , 'C':rate_AC*prob_C, 'G':rate_AG*prob_G, 'T':rate_AT*prob_T},
        'C': {'A':rate_AC*prob_A, 'C':0             , 'G':rate_CG*prob_G, 'T':rate_CT*prob_T},
        'G': {'A':rate_AG*prob_A, 'C':rate_CG*prob_C, 'G':0             , 'T':rate_GT*prob_T},
        'T': {'A':rate_AT*prob_A, 'C':rate_CT*prob_C, 'G':rate_GT*prob_G, 'T':0             }
    }

    for letter1 in R:
        v = 0
        for letter2 in R[letter1]:
            v += R[letter1][letter2]
        R[letter1][letter1] = -v

    # print('\tA\tC\tG\tT')
    # for letter1 in R:
    #     line = letter1
    #     for letter2 in R[letter1]:
    #         line += '\t' + str(round(R[letter1][letter2], 3))
    #     print(line)
    #-------------------------------------------------------

    #normalize R
    #-------------------------------------------------------
    transRate = 0
    for letter in R:
        transRate -= R[letter][letter] * initProb[letter]

    for letter1 in R:
        for letter2 in R[letter1]:
            R[letter1][letter2] /= transRate

    # print('transition rate:' + str(transRate))

    # print('\tA\tC\tG\tT')
    # for letter1 in R:
    #     line = letter1
    #     for letter2 in R[letter1]:
    #         line += '\t' + str(round(R[letter1][letter2], 3))
    #     print(line)
    #-------------------------------------------------------

    #convert to numpy array
    #-------------------------------------------------------
    letters = ['A', 'C', 'G', 'T']
    Rmat = np.zeros((4,4),dtype=float)
    for (letter1, i) in zip(letters,range(len(letters))):
        for (letter2, j) in zip(letters,range(len(letters))):
            Rmat[i][j] = R[letter1][letter2]

    # print(Rmat)
    #-------------------------------------------------------
    
    #generate all sequences from root to leaves
    #-------------------------------------------------------
    seqs = dict() # This will be your output (keys = leaf labels (str) and values = sequences (str))
    for node in tree.traverse_preorder():
        if node.is_root():
            seqs[node] = root_seq

        for child in node.child_nodes():
            P = expm(child.get_edge_length()*Rmat)
            subs = dict()
            ind = dict()
            for (letter, i) in zip(letters,range(len(letters))):
                cnt = seqs[node].count(letter)
                subs[letter] = np.random.choice(letters, size=cnt, p=P[i])
                ind[letter] = 0

            newSeq = ''
            for letter in seqs[node]:
                newSeq += subs[letter][ind[letter]]
                ind[letter] += 1

            seqs[child] = newSeq
    #-------------------------------------------------------

    #extract leaf sequences
    #-------------------------------------------------------
    leafSeqs = dict()
    for node in tree.traverse_leaves():
        leafSeqs[node.get_label()] = seqs[node]
    #-------------------------------------------------------
    
    return leafSeqs

def genTree(depth):
    tree = read_tree_newick('root;')
    for n in tree.traverse_leaves():
        root = n

    for i in range(depth):
        leaves = []
        for n in tree.traverse_leaves():
            leaves.append(n)

        for n in leaves:
            n.add_child(Node())
            n.add_child(Node())

    for n in tree.traverse_preorder():
        if not n.is_root():
            n.set_edge_length(random.random())

    i = 1
    for n in tree.traverse_leaves():
        n.set_label(str(i))
        i = i + 1

    # tree.draw(show_labels=True)

    return tree

def genGen(tree, genLen):
    gtr_probs = [0.292187, 0.211771, 0.229896, 0.266146]
    gtr_rates = [1.066667, 0.2, 0.3333334, 0.3333334, 0.4, 1]
    root_seq = random_seq(genLen, gtr_probs)
    seqs = evolve(tree, root_seq, gtr_probs, gtr_rates)

    return seqs