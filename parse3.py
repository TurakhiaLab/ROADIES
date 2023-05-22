import sys,glob,os
import argparse
import random
import numpy as np
#getting arguments needed
parser = argparse.ArgumentParser(description='choosing weighted sampling')
parser.add_argument('-t',type=float,default = 0.6)
parser.add_argument('-k',type=int,default=1000)
parser.add_argument('--quartets',default='freqQuad.csv')
parser.add_argument('--weighted',type=float,default=0)
#parser.add_argument('')
parser.add_argument('--genomes',default='/home/roadies-datasets/birds')
parser.add_argument('--output',default="results/samples/counts.csv")
parser.add_argument('--num_align',type=float,default=0.5)
args = parser.parse_args()
t = args.t
k = args.k
w = args.weighted
q = args.quartets
o = args.output
a = args.num_align
genomes = args.genomes
print(t,k,w,q)
counts ={}
#getting list of species names
for filename in glob.glob(os.path.join(genomes,'*.fa')):
    s = filename.split('/')
    name = s[len(s)-1]
    species = name.replace(".maf","")
    counts[species] = 0
master_species = list(counts.keys())
num_align = int(a*len(master_species))
# if we are doing weighted sampling
if w != 0:
    print("Doing weighted sampling")
    #get number of weighted 'genes'
    num_weighted = int(k * w)
    lines = open(q,'r').readlines()
    scores = {}
    #go through freqQuad.csv and get scores of t1
    for i in range(len(lines)):
        if i % 3 == 0:
            s = lines[i].split('\t')
            scores[i] = float(s[3])
    #sort list of quartets by score
    sorted_scores = sorted(scores.items(),key=lambda x:x[1])
    #only use quartets with a score lesser than t
   
 
    num_sampled = 0
    idx = 0
    weights = []
    scores = [1/x[1] for x in sorted_scores]
    m = np.sum(scores)
    for i in range(len(scores)):
        weights.append(scores[i]/m)
    print(weights)
    print(np.sum(weights))
    species_list = []
    #keep on sampling from low scoring quartet branches until we sample required number of weighted genes
    while num_sampled < num_weighted:
        to_align = []
        while len(to_align) < num_align:
          #  print(len(species),len(weights))
            idx = np.random.choice(len(weights),1,p=weights)[0]
         #   print(idx)
      #  print(num_sampled,num_weighted)
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
            #iterate through each branch of quartet and randomly select one species from that branch
            for node in nodes:
                n = node.replace('{','')
                n = n.replace('}','')
                s = n.split(',')
                nodes_list.append(s)
            itr = 0
            while True:
                node = nodes_list[itr]
                #pick random species on branch
                s = random.randint(0,len(node)-1)
                #only add if species not found already
                if node[s] not in to_align:
                    to_align.append(node[s])
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
    if num_sampled < k:
        for i in range(k-num_sampled):
            r = random.randint(0,len(master_species)-1)
            species = master_species[r]
            rand_list = random.sample(master_species,num_align)
            print(species,rand_list)
           # counts[species] += 1
            species_list.append(rand_list)
    print(species_list)
    ctr = 0
    mapping = []
    for i in range(len(species_list)):
        r = random.randint(0,len(species_list[i])-1)
        if ".fa" not in species_list[i][r]:
            species = species_list[i][r]+".fa"
        else:
            species = species_list[i][r]
        #print(species)
        counts[species] += 1
        ctr += 1
        mapping.append((species,i))
    mapping_sorted = sorted(mapping)
    #print(counts)
   # print(ctr)
   # print(mapping_sorted)
    w = open(o,'w')
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
      #  print("gene id",i,species,list_of_species)
    id_stops.append(k)
   # print(id_stops)
   # print(species_stops)

    for i in range(len(species_stops)):
        species = species_stops[i]
        start = id_stops[i]
        end=id_stops[i+1]
        print(species,start,end)
        w.write(species+','+str(start)+','+str(end)+'\n')
        print(species_list[end-1])
    #print(counts)
    #print(len(counts))
    left_to_sample = k - num_sampled
    #the rest of sampling is completely random
    #for i in range(left_to_sample):
       
   # print(counts)
else:
    print("Doing unbiased sampling")
    for i in range(k):
        r = random.randint(0,len(species)-1)
        print(species[r])
        counts[species[r]] += 1
    print(counts)
       # print(quartet)
with open(o,'w') as w:
    for species in counts:
        w.write(species +','+str(counts[species])+'\n')