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
k = e-s
if k < 1:
    print("Invalid indices")
    exit(1)
print("Number of regions:", k)
print("ID START: "+(str(s))+", ID END: "+str(e))
print("Region Length:", l)
print("Input File:" +path)
records = list(SeqIO.parse(path,"fasta"))
sequence=[]
index = s
threshold = int(t*l)
with alive_bar(k) as bar:
    for i in range(k):
        notFound =True
        while notFound:
            c_seq=random.randint(0,len(records)-1)
            seq_len=len(records[c_seq].seq)
            if l>seq_len:
                continue
            start = random.randint(0,seq_len-l)
            frag = records[c_seq].seq[start:(start+l)]
            upper = 0
            for i in range(len(frag)):
                if frag[i] =='N':
                    break
                if frag[i]=='A' or frag[i] =='T' or frag[i]=='G' or frag[i]=='T':
                    upper = upper +1
            if upper >= threshold:
                record = SeqRecord(frag, "gene_"+str(index),"","")
                sequence.append(record)
                index = index+1
                break
        time.sleep(0.01)
        bar()

SeqIO.write(sequence, output, "fasta")

