from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import time
import random
# import OS module
import numpy as np
import os
from alive_progress import alive_bar 
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
with alive_bar(k) as bar:
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
            if j ==4:
                print("no good samples defaulting: "+frag)
                record = SeqRecord(frag, "fragment_%i" % (i + 1), "", "")
                sequence.append(record)
                break
            if 'a' not in frag and "N" not in frag and 't' not in frag and 'g' not in frag and 'c' not in frag:
                print(frag)
                record = SeqRecord(frag, "fragment_%i" % (i + 1), "", "")
                sequence.append(record)
                break
            

        time.sleep(0.01)
        bar()

SeqIO.write(sequence, "out.fasta", "fasta")

