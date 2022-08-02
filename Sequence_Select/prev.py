from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import random
import numpy as np

# import OS module
import os
import time

start_time = time.time()

# Get the list of all files and directories
path = "./../genomes"
dir_list = os.listdir(path)

# Arguments passed

k= random.randint(0,10) 
l= random.randint(0,300) 


for i in range(1, len(sys.argv)):
    if sys.argv[i]=="-l":
        l=int(sys.argv[i+1])
    elif sys.argv[i]=="-k":
        k=int(sys.argv[i+1])
     
print("Number of regions:", k)
print("Region Length :", l)
print("\n")

print(len(dir_list))
filenum = np.random.randint(0, len(dir_list), k)

filenum.sort()
sequence=[]
for i in range(k):
    if i==0 or filenum[i]!=filenum[i-1]:
        f=filenum[i]
        SeqFile=path+"/"+dir_list[f]
        records = list(SeqIO.parse(SeqFile, "fasta"))
    for j in range(5):
        c_seq=random.randint(0,len(records)-1)
        seq_len=len(records[c_seq].seq)
        start=random.randint(0,seq_len-l)
        frag=records[c_seq].seq[start:(start+l)]
        if "a" not in frag or "N" not in frag or "t" not in frag or "g" not in frag or "c" not in frag:
            break
    record = SeqRecord(frag, "fragment_%i" % (i + 1), "", "")
    sequence.append(record)

SeqIO.write(sequence, "out.fasta", "fasta")
print("--- %s seconds ---" % (time.time() - start_time))
