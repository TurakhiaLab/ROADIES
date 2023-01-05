from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import time
import random
# import OS module
import numpy as np
import os
from alive_progress import alive_bar 
import argparse, os
# Get the list of all files and directories
# Arguments passed
parser = argparse.ArgumentParser(description='genome selection')
parser.add_argument('-l',type=int,default = 1000)
parser.add_argument('-s',type=int,required=True)
parser.add_argument('-e',type=int,required=True)
parser.add_argument('-t',type=float,default=0.7)
parser.add_argument('--input',default="./../genomes")
parser.add_argument('--output',default="./index.csv")
args = parser.parse_args()
s= args.s
e= args.e
l= args.l
t = args.t
path = args.input
output=args.output
k = e-s+1
if k < 1:
    print("Invalid indices")
    exit(1)
print("Number of regions:", k)
print("ID START: "+(str(s))+", ID END: "+str(e))
print("Region Length:", l)
print("Input File:" +path)
records = list(SeqIO.parse(path,"fasta"))
total_length = 0
indices = []
for record in records:
    total_length += len(record.seq)
    indices.append(total_length)
#print(indices)
sequence=[]
index = s
threshold = int(t*l)
#keeps track of time
with alive_bar(k) as bar:
    #for specified number of genes
    for i in range(k):
        notFound =True
        count = 0
        #until a passing sample is found
        while notFound:
            #print(k)
            k = random.randint(0,total_length-l)
            loc = 0
            idx = 0
            start = 0 
            for index in indices:
                if k > index:
                    loc += 1 
                else:
                    start = index - k
                    break
            if start < l:
                count += 1
                continue
            frag = records[loc].seq[start:start+l]
            print(len(frag))
            print(frag)
            upper = 0
            uppercase = ['A','T','C','G']
            accepted = ['a','t','c','g']
            for f in frag:
                if f not in accepted and f not in uppercase:
                    break
                elif f in uppercase:
                    upper += 1
            if upper >= threshold:
                record = SeqRecord(frag, "gene_"+str(index),"","")
                sequence.append(record)
                index = index+1
                break
            count = count+1
        print('# of resampling: '+str(count))
        time.sleep(0.01)
        bar()

SeqIO.write(sequence, output, "fasta")

