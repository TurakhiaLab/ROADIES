import re
import os, glob
import sys
import argparse
from Bio import SeqIO
from operator import itemgetter
parser = argparse.ArgumentParser(description='get gene fasta')
parser.add_argument('-k',type=int,default=400)
parser.add_argument('--path')
parser.add_argument('--outdir')
parser.add_argument('-m',type=int,default=4)
parser.add_argument('-d',type=int,default=5)
args = parser.parse_args()
path = args.path
outdir = args.outdir
k = args.k
d= args.d
m = args.m
for filename in glob.glob(os.path.join(path,'*.maf')):
    with open(os.path.join(os.getcwd(),filename),'r') as f:
        s = filename.split('/')
        name = s[len(s)-1]
        species = name.replace(".maf","")
        lines = f.readlines()
        found = []
        genes = {}
        for l in range(15,len(lines)):
            if (l-15)%4 == 0:
                gene_line = lines[l+1].split()
                gene = gene_line[1]
                gene_s = gene.split('_')
                gene_id = gene_s[1]
                score_line = lines[l-1].split()
                score_expr = score_line[1].split('=')
                score = int(score_expr[1])
                if gene_id not in genes:
                    genes[gene_id] = [(score,l)]
                else:
                    genes[gene_id].append((score,l))
        for gene in genes:
            gene_list = genes[gene]
            if len(gene_list) < d:
                for i in range(len(gene_list)):
                    l = gene_list[i][1]
                    seq_line = lines[l].split()
                    seq = seq_line[len(seq_line)-1]
                    index = species+'_'+seq_line[2]
                    with open(outdir+'/gene_'+gene+'.fa','a') as w:
                        #print("adding ",gene)
                            w.write('>'+index+'\n')
                            w.write(seq+'\n')
                    w.close()
                    with open(outdir+"/mapping.txt",'a') as w2:
                        w2.write(index+' ' +species+'\n')
                    w2.close()
            else:
                max_scores = list(range(0,d))
                #print(max_scores)
                m = 0
                ms = gene_list[0][0]
                for i in range(1,d):
                    score = gene_list[i][0]
                    if score < ms:
                        ms = score
                        m = i
                for g in range(d,len(gene_list)):
                    score = gene_list[g][0]
                    if score > ms:
                        max_scores[m] = g
                        m = 0
                        ms= gene_list[max_scores[0]][0]
                        for i in range(d):
                            idx = max_scores[i]
                            score = gene_list[idx][0]
                            if score < ms:
                                ms = score
                                m = i
                for i in max_scores:
                    l = gene_list[i][1]
                    seq_line = lines[l].split()
                    seq = seq_line[len(seq_line)-1]
                    index = species+'_'+seq_line[2]
                    with open(outdir+'/gene_'+gene+'.fa','a') as w:
                        #print("adding ",gene)
                            w.write('>'+index+'\n')
                            w.write(seq+'\n')
                    w.close()
                    with open(outdir+"/mapping.txt",'a') as w2:
                        w2.write(index+' ' +species+'\n')
                    w2.close()
                #print(species,gene,counts[gene])

#for filename in glob.glob(os.path.join(outdir,'*.fa')):
 #   print(filename)
  #  records = list(SeqIO.parse(filename,"fasta"))
   # print(len(records))
    #if len(records)< m:
     #   print(filename)
      #  os.remove(filename)
