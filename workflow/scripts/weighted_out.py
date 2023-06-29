# sequence_merge.py concats all the genes together and gathers stats
# REQUIRES: Seaborn and Matplotlib
# import shutil
import glob, os, sys

# import seaborn as sns
# import matplotlib.pyplot as plt
genes = sys.argv[1]
species_id = sys.argv[2]
species = sys.argv[3]
out = sys.argv[4]
# stats = sys.argv[3]
species_ids = {}
id_dict = {}
g = open(genes, "r").readlines()
curr_id = g[0].strip().replace(">gene_", "")
curr_seq = ""
# print(curr_id)
for i in range(1, len(g)):
    # print(g[i])
    # start of new seq
    if ">gene_" in g[i]:
        #   print(g[i])
        id_dict[curr_id] = curr_seq
        # print("found")
        # print(curr_id,curr_seq)
        curr_id = g[i].strip().replace(">gene_", "")
        # print(curr_id)
        curr_seq = ""
    else:
        if i == len(g) - 1:
            curr_seq += g[i].strip()
            id_dict[curr_id] = curr_seq
        else:
            curr_seq += g[i].strip()
# print(curr_seq)
# print(id_dict)
# print(len(id_dict))

lines = open(species_id, "r").readlines()
# print(len(lines))
for i in range(len(lines)):
    s = lines[i].strip().split(",")
    if s[0] == species:
        with open(out, "w") as w:
            ids = s[1:]
            for id in ids:
                w.write(">gene_{0}\n".format(id))
                w.write(id_dict[id] + "\n")

        # species_list[]
# for species in species_ids:
