#This sampling script samples specific number of genes, each of a given specific length, from input files and saves it in fasta format
#REQUIRES: Biopython
#USAGE: `python workflow/scripts/sampling.py {args}`
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
import argparse
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
    #for specified number of genes
count = 0
for i in range(k-1):
    notFound =True
    sample = ''
    #until a passing sample is found
    while notFound:
        #print(k)
        k = random.randint(0,total_length-l)
        loc = 0
        idx = 0
        start = 0 
        for i in indices:
            if k > i:
                loc += 1 
            else:
                start = i - k
                break
        if start < l:
            count += 1
            continue
        frag = records[loc].seq[start:start+l]
        upper = 0
        uppercase = ['A','T','C','G']
        accepted = ['a','t','c','g']
        accept = True
        if 'N' in frag:
            continue
        for f in frag:
            if f not in uppercase:
                if f not in accepted:
                    accept = False
                    break
            else:
                upper += 1
        if accept == False:
            continue
        if upper >= threshold:
            sample = SeqRecord(frag, "gene_"+str(index),"","")
            break
        else:
            count = count+1
    sequence.append(sample)
    index = index+1
print('Total # of resampling: '+str(count))
SeqIO.write(sequence, output, "fasta")

