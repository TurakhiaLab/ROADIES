import re
import os, glob
import sys
import argparse
parser = argparse.ArgumentParser(description='get gene fasta')
parser.add_argument('-k',type=int,default=400)
parser.add_argument('--path')
parser.add_argument('--outdir')
parser.add_argument('-m',type=int,default=4)
args = parser.parse_args()
path = args.path
outdir = args.outdir
k = args.k
m = args.m

counts = {}
isin=[]
for filename in glob.glob(os.path.join(path,'*.maf')):
    with open(os.path.join(os.getcwd(),filename),'r') as f:
        s = filename.split('/')
        name = s[len(s)-1]
        species = name.replace(".maf","")
        lines = f.readlines()
        print(species)
        for l in range(15,len(lines)):
            if (l-15)%4 == 0:
                split = lines[l].split()
                seq = split[len(split)-1]
                seq = seq.replace("-","")
                n = lines[l+1]
                split2 = n.split()
                gene = split2[1]
                s = gene.split('_')
                num = int(s[1])
                if num not in isin:
                    isin.append(num)
                if gene in counts:
                    counts[gene]= counts[gene]+1
                else:
                    counts[gene] = 1
                with open(outdir+'/'+gene+'.fa','a') as w:
                    #print("adding ",gene)
                    w.write('>'+species+"_"+split[2]+'\n')
                    w.write(seq+'\n')
                #print(species,gene,counts[gene])
print(counts)
for i in range(1,int(k)):
    if i not in isin:
        print("no matches for gene: ",i)
for filename in glob.glob(os.path.join(path,'*.fa')):
    records = list(Seq.IO.parse(filename),"fasta")
    if len(records)< m:
        os.remove(path)
