from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import random
import argparse
# import OS module
import os
import time
import numpy as np

start_time = time.time()

 
# Get the list of all files and directories
# Arguments passed
k= random.randint(0,10) 
l= random.randint(0,300) 

def isPath(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

parser = argparse.ArgumentParser()
parser.add_argument('-path', type=isPath, required=True)
parser.add(argument('-k', type=int, default=k))
parser.add(argument('-l', type=int, default=l))
args = parser.parse_args()
print("Input Directory Path:"+ args.path)
print("Number of regions:" +str(args.k))
print("Region Length :"+ str(args.l))
print("\n")
path = args.path
dir_list = os.listdir(path)

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
main()
print("--- %s seconds ---" % (time.time() - start_time))
