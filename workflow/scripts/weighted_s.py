#USAGE: `python workflow/scripts/weighted_s.py [arguments described in argparse help]
import sys,glob,os
import argparse
import random
import numpy as np
from collections import Counter
#getting arguments needed
parser = argparse.ArgumentParser(description='choosing weighted sampling')
#parser.add_argument('-t',type=float,default = 0.6,help="")
parser.add_argument('-k',type=int,default=1000,help="number of genes to sample total (weighed + unweighted)")
parser.add_argument('--quartets',default='freqQuad.csv',help="location of ASTRAL-PRO quartet file")
parser.add_argument('--weighted',type=float,default=0,help="percentage of genes using weighted sampling; default is 0 which is totally unweighted")
parser.add_argument('--species_out',default="species_id.csv")
parser.add_argument('--genomes',default='/home/roadies-datasets/birds',help="genome/assembly/fasta directory")
parser.add_argument('--id_out',default="gene_id.csv",help="output file for gene ids")
parser.add_argument('--align_out',default="species_lists.csv",help="output file for list of species in weighted sampling")
parser.add_argument('--to_align',type=float,default=0.5,help="percentage of total species to run lastz with")
args = parser.parse_args()
#t = args.t
KREG = args.k
WEIGHTED = args.weighted
QUARTET_FILE = args.quartets
GENE_ID = args.id_out
SPECIES_LIST = args.align_out
TO_ALIGN = args.to_align
GENOMES = args.genomes
SPECIES_ID = args.species_out
counts ={}
#getting list of species names
for filename in glob.glob(os.path.join(GENOMES,'*.fa')):
    s = filename.split('/')
    name = s[len(s)-1]
    species = name.replace(".fa","")
    counts[species] = 0
#MASTER LIST OF SPECIES NAMES
MASTER_SPECIES = list(counts.keys())
#number of species we are aligning to each gene in lastz
num_align = int(TO_ALIGN*len(MASTER_SPECIES))
# if we are doing weighted sampling
if WEIGHTED != 0:
    #get number of weighted 'genes'
    num_weighted = int(KREG * WEIGHTED)
    print("Doing weighted sampling on {0} of {1} genes".format(num_weighted,KREG))
    #getting scores of quartet to assign weights
    lines = open(QUARTET_FILE,'r').readlines()
    scores = {}
    #go through freqQuad.csv and get pp scores of t1 of each quartet
    for i in range(len(lines)):
        if i % 3 == 0:
            s = lines[i].split('\t')
            scores[i] = float(s[3])
    #sort list of quartets by score
    sorted_scores = sorted(scores.items(),key=lambda x:x[1])
    #only use quartets with a score lesser than t
    num_sampled = 0
    idx = 0
    #assign weighted scheme of (1/c_i)/summation_i of (1/c_i)
    weights = []
    scores = [1/x[1] for x in sorted_scores]
    m = np.sum(scores)
    for i in range(len(scores)):
        weights.append(scores[i]/m)
    #print(weights)
    #now we are going to get the list of species for each gene
    species_list = []
    #keep on sampling from low scoring quartet branches until we sample required number of weighted genes
    while num_sampled < num_weighted:
        #list of species to align gene to
        to_align = []
        #keep on sampling species until we get to the number we want to align using lastz
        while len(to_align) < num_align:
            idx = np.random.choice(len(weights),1,p=weights)[0]
            line = lines[sorted_scores[idx][0]]
            quartet = line.split('\t')[2].strip()
            s2 = quartet.split('#')
            s3 = s2[0].split('|')
            s4 = s2[1].split('|')
            R = s3[0]
            L = s3[1]
            S = s4[0]
            O = s4[1]
            nodes = [R,L,S,O]
            nodes_list = []
            #just reformatting for ease
            for node in nodes:
                n = node.replace('{','')
                n = n.replace('}','')
                s = n.split(',')
                nodes_list.append(s)
            itr = 0
            #iterate through each branch of quartet and randomly select one species from that branch
            while True:
                node = nodes_list[itr]
                #pick random species on branch
                s = random.randint(0,len(node)-1)
                #only add if species not found already
                if node[s] not in to_align:
                    to_align.append(node[s])
                    #check if number of species is less than limit 
                    if len(to_align) == num_align:
                        break
                #remove species from the node list
                del node[s]
                if len(node) == 0:
                    break
                # increment node counter or reset if its at O
                else:
                    if itr == 3:
                        itr= 0
                    else:
                        itr += 1

        #print(to_align)
        species_list.append(to_align)
        num_sampled += 1
    #done with weighted sampling
    #now do the rest of unweighted sampling
    if num_sampled < KREG:
        for i in range(KREG-num_sampled):
            r = random.randint(0,len(MASTER_SPECIES)-1)
            species = MASTER_SPECIES[r]
            rand_list = random.sample(MASTER_SPECIES,num_align)
           # print(species,rand_list)
           # counts[species] += 1
            species_list.append(rand_list)
    #we should have a list of k species of size TO_ALIGN at this point with (WEIGHTED * k) of them being from weighted sampling
    #print(species_list)
    ctr = 0
    #a list of pairs where each pair is (species,index on species_list)
    mapping = []
    #picking 1 random species from each species list
    for i in range(len(species_list)):
        r = random.randint(0,len(species_list[i])-1)
        #if ".fa" not in species_list[i][r]:
          #  species = species_list[i][r]+".fa"
        #else:
        species = species_list[i][r]
        #print(species)
        counts[species] += 1
        ctr += 1
        mapping.append((species,i))
    #sort the tuples by species
    mapping_sorted = sorted(mapping)
    #iterating through mapping to get index where species changes
    curr = mapping_sorted[0][0]
    species_stops=[curr]
    id_stops = [1]
    for i in range(len(mapping_sorted)):
        species = mapping_sorted[i][0]
        if species != curr:
            stop = i+1
            id_stops.append(stop)
            curr = species
            species_stops.append(species)
        idx = mapping_sorted[i][1]
        list_of_species = species_list[idx]
    id_stops.append(KREG)
    #writing gene ids
    with open(GENE_ID,'w') as w:
        for i in range(len(species_stops)):
            species = species_stops[i]
            start = id_stops[i]
            end=id_stops[i+1]
            #print(species,start,end)
            w.write(species+','+str(start)+','+str(end)+'\n')
    #print list of species lists
    with open(SPECIES_LIST,'w') as w2:
        for i in range(len(species_list)):
            li = species_list[i]
            for j in range(len(li)):
                if j == len(li)-1:
                    w2.write(li[j]+'\n')
                else:
                    w2.write(li[j]+',')
    species_ids = {}
    for i in range(KREG):
        mapping = mapping_sorted[i]
        idx = mapping[1]
        li = species_list[idx]
        for species in li:
           # print(species)
            if species not in species_ids:
                species_ids[species] = [i+1]
            else:
                species_ids[species].append(i+1)
    #for species in species_ids:
        #print(species,len(species_ids[species]))
    with open(SPECIES_ID,'w') as w3:
        for species in species_ids:
            #print(species)
            w3.write(species+',')
            #print(species_ids[species])
            for i in range(len(species_ids[species])):
                if i == len(species_ids[species])-1:
                    w3.write(str(species_ids[species][i])+'\n')
                else:
                    w3.write(str(species_ids[species][i])+',')

else:
    print("weight must be greater than 0")
    exit()
