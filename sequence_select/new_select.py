from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
import time
import random
# import OS module
import numpy as np
import os
from alive_progress import alive_bar 
import re
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
""" print("Number of regions:", k)
print("Region Length :", l)
print("Input Directory:" +path)
print("\n") """
def CheckLength(letter):
    return True if len(letter)>=l else False


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
            if l>seq_len:
                continue
            new_list=list(filter(CheckLength,re.split("[a,c,g,t,N]", str(records[c_seq].seq))))
            if len(new_list)==0:
                continue
            chosen=random.randint(0,len(new_list)-1)
            start = random.randint(0,len(new_list[chosen])-l)
            frag = Seq(new_list[chosen][start:(start+l)])
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