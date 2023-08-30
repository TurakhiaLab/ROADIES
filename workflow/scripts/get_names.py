# get_names.py reads a fasta file and outputs n subsets of that file by chromosome name
# REQUIREMENTS: Biopython
# USAGE: `python workflow/scripts/get_names.py [path] [output_dir] [number of subsets]`
import os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

path = sys.argv[1]
n = int(sys.argv[3])
outdir = sys.argv[2]
s = path.split("/")
file = s[len(s) - 1]
filename = file.replace(".fa", "")
print(filename)
records = list(SeqIO.parse(path, "fasta"))
num = len(records)
print(filename + " has " + str(num) + " records")
total = 0
for record in records:
    total += len(record.seq)
print(filename + " has " + str(total) + " number of base pairs")
indices = []
subsets = [[]]
for i in range(n - 1):
    index = int(((i + 1) / n) * total)
    indices.append(index)
    subsets.append([])
print(indices)
counter = 0
j = 0
for i in range(len(records)):
    counter += len(records[i].seq)
    if j < n - 1:
        if counter >= indices[j]:
            j += 1
        else:
            subsets[j].append(records[i].name)
    else:
        subsets[j].append(records[i].name)
for i in range(len(subsets)):
    fs = "/" + filename + "." + str(i) + ".subset"
    print(fs)
    print("has " + str(len(subsets[i])) + " records")
    with open(outdir + fs, "w") as w:
        for record in subsets[i]:
            # print(name)
            w.write(record + "\n")
