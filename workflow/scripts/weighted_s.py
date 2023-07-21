# USAGE: `python workflow/scripts/weighted_s.py [arguments described in argparse help]
import sys, glob, os
import argparse
import random
import numpy as np
from collections import Counter

# getting arguments needed
parser = argparse.ArgumentParser(description="choosing weighted sampling")
# parser.add_argument('-t',type=float,default = 0.6,help="")
parser.add_argument(
    "-k",
    type=int,
    default=1000,
    help="number of genes to sample total (weighed + unweighted)",
)
parser.add_argument(
    "--weighted",
    type=float,
    default=0,
    help="percentage of genes using weighted sampling; default is 0 which is totally unweighted",
)
parser.add_argument(
    "--genomes",
    default="/home/roadies-datasets/birds",
    help="genome/assembly/fasta directory",
)

parser.add_argument(
    "--to_align", type=int, required=True, help="# species to align per gene"
)
parser.add_argument(
    "--num_species", type=int, required=True, help="# species to align per gene"
)
args = parser.parse_args()
# t = args.t
KREG = args.k
WEIGHTED = args.weighted
TO_ALIGN = args.to_align
GENOMES = args.genomes
NUM_SPECIES = args.num_species
counts = {}
# getting list of species names
for filename in glob.glob(os.path.join(GENOMES, "*.fa")):
    s = filename.split("/")
    name = s[len(s) - 1]
    species = name.replace(".fa", "")
    counts[species] = 0
# MASTER LIST OF SPECIES NAMES
MASTER_SPECIES = list(counts.keys())
# number of species we are aligning to each gene in lastz
# if we are doing weighted sampling
if NUM_SPECIES != TO_ALIGN:
    # get number of weighted 'genes'
    num_weighted = int(KREG * WEIGHTED)
    num_sampled = 0
    species_list = []
    if num_weighted != 0:
        print("Doing weighted sampling on {0} of {1} genes".format(num_weighted, KREG))
    # getting scores of quartet to assign weights
        lines = open("freqQuad.csv","r").readlines()
        scores = {}
    # go through freqQuad.csv and get pp scores of t1 of each quartet
        for i in range(len(lines)):
            if i % 3 == 0:
                s = lines[i].split("\t")
                #print(s[4],s[5])
                support = float(s[4])
                teb = float(s[5])
                if teb == 0:
                    continue
                scores[i] = support / teb
    # sort list of quartets by score
        sorted_scores = sorted(scores.items(), key=lambda x: x[1])
    # only use quartets with a score lesser than t
        
        idx = 0
    # assign weighted scheme of (1/c_i)/summation_i of (1/c_i)
        weights = []
        scores = [1 / x[1] for x in sorted_scores]
        m = np.sum(scores)
        for i in range(len(scores)):
            weights.append(scores[i] / m)
    # print(weights)
    # now we are going to get the list of species for each gene
       
        # keep on sampling from low scoring quartet branches until we sample required number of weighted genes
        while num_sampled < num_weighted:
            # list of species to align gene to
            species_to_align = []
            # keep on sampling species until we get to the number we want to align using lastz
            while len(species_to_align) < TO_ALIGN:
                idx = np.random.choice(len(weights), 1, p=weights)[0]
                line = lines[sorted_scores[idx][0]]
                quartet = line.split("\t")[2].strip()
                s2 = quartet.split("#")
                s3 = s2[0].split("|")
                s4 = s2[1].split("|")
                R = s3[0]
                L = s3[1]
                S = s4[0]
                O = s4[1]
                nodes = [R, L, S, O]
                nodes_list = []
                # just reformatting for ease
                for node in nodes:
                    n = node.replace("{", "")
                    n = n.replace("}", "")
                    s = n.split(",")
                    nodes_list.append(s)
                itr = 0
                # iterate through each branch of quartet and randomly select one species from that branch
                while True:
                    node = nodes_list[itr]
                    # pick random species on branch
                    s = random.randint(0, len(node) - 1)
                    # only add if species not found already
                    if node[s] not in species_to_align:
                        species_to_align.append(node[s])
                        # check if number of species is less than limit
                        if len(species_to_align) == TO_ALIGN:
                            break
                    # remove species from the node list
                    del node[s]
                    if len(node) == 0:
                        break
                    # increment node counter or reset if its at O
                    else:
                        if itr == 3:
                            itr = 0
                        else:
                            itr += 1

            # print(to_align)
            species_list.append(species_to_align)
            num_sampled += 1
    # done with weighted sampling
    # now do the rest of unweighted sampling
    if num_sampled < KREG:
        for i in range(KREG - num_sampled):
            r = random.randint(0, len(MASTER_SPECIES) - 1)
            species = MASTER_SPECIES[r]
            rand_list = random.sample(MASTER_SPECIES, TO_ALIGN)
            # print(species,rand_list)
            # counts[species] += 1
            species_list.append(rand_list)
    # we should have a list of k species of size TO_ALIGN at this point with (WEIGHTED * k) of them being from weighted sampling
    # print(species_list)
    ctr = 0
    # a list of pairs where each pair is (species,index on species_list)
    mapping = []
    # picking 1 random species from each species list
    for i in range(len(species_list)):
        r = random.randint(0, len(species_list[i]) - 1)
        # if ".fa" not in species_list[i][r]:
        #  species = species_list[i][r]+".fa"
        # else:
        species = species_list[i][r]
        # print(species)
        counts[species] += 1
        ctr += 1
        mapping.append((species, i))
    # in case a species doesn't get sampled
    for species in counts:
        if counts[species] == 0:
            mapping.append((species, len(species_list)))
            species_list.append(random.sample(MASTER_SPECIES, TO_ALIGN))
    # sort the tuples by species
    mapping_sorted = sorted(mapping)
    # iterating through mapping to get index where species changes
    curr = mapping_sorted[0][0]
    species_stops = [curr]
    id_stops = [1]
    for i in range(len(mapping_sorted)):
        species = mapping_sorted[i][0]
        if species != curr:
            stop = i + 1
            id_stops.append(stop)
            curr = species
            species_stops.append(species)
        idx = mapping_sorted[i][1]
        list_of_species = species_list[idx]
    id_stops.append(len(mapping_sorted))
    # writing gene ids
    with open("gene_ids.csv", "w") as w:
        for i in range(len(species_stops)):
            species = species_stops[i]
            start = id_stops[i]
            end = id_stops[i + 1]
            # print(species,start,end)
            w.write(species + "," + str(start) + "," + str(end) + "\n")
    # print list of species lists
    with open("species_lists.csv", "w") as w2:
        for i in range(len(species_list)):
            li = species_list[i]
            for j in range(len(li)):
                if j == len(li) - 1:
                    w2.write(li[j] + "\n")
                else:
                    w2.write(li[j] + ",")
    species_ids = {}
    for i in range(KREG):
        mapping = mapping_sorted[i]
        idx = mapping[1]
        li = species_list[idx]
        for species in li:
            # print(species)
            if species not in species_ids:
                species_ids[species] = [i + 1]
            else:
                species_ids[species].append(i + 1)
    # for species in species_ids:
    # print(species,len(species_ids[species]))
    with open("species_ids.csv", "w") as w3:
        for species in species_ids:
            # print(species)
            w3.write(species + ",")
            # print(species_ids[species])
            for i in range(len(species_ids[species])):
                if i == len(species_ids[species]) - 1:
                    w3.write(str(species_ids[species][i]) + "\n")
                else:
                    w3.write(str(species_ids[species][i]) + ",")

else:
    print("NO SAMPLING")
    exit()
