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
# Arguments passed
k= random.randint(0,10) 
l= random.randint(0,300) 
for i in range(1, len(sys.argv)):
    if sys.argv[i]=="-l":
        l=int(sys.argv[i+1])
    elif sys.argv[i]=="-k":
        k=int(sys.argv[i+1])
path = sys.argv[1]
print("Number of regions:", k)
print("Region Length :", l)
print("Input Directory:" +path)
print("\n")
records = list(SeqIO.parse(path,"fasta"))
sequence=[]
name = ''
split = path.split('/')
for s in split:
    if 'fa' in s:
        name = s
with alive_bar(k) as bar:
    for i in range(k):
        for j in range(5):
            c_seq=random.randint(0,len(records)-1)
            seq_len=len(records[c_seq].seq)
            start = random.randint(0,seq_len-l)
            frag = records[c_seq].seq[start:(start+l)]
            if j ==4:
                record = SeqRecord(frag, name+':'+str(start)+':'+str(i),"","" )
                if "N" not in record:
                    sequence.append(record)
                break
            if 'a' not in frag and "N" not in frag and 't' not in frag and 'g' not in frag and 'c' not in frag:
                record = SeqRecord(frag, name+str(i), "", "")
                sequence.append(record)
                break
            
        time.sleep(0.01)
        bar()
frag = ''
split = path.split('/')
for s in split:
    if 'fa' in s:
        frag = s
SeqIO.write(sequence, "frag_"+ s, "fasta")

