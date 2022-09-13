import sys, glob, os
from Bio import SeqIO
import shutil
from Bio.Align.Applications import MafftCommandline
path = sys.argv[1]
out = sys.argv[2]
n = int(sys.argv[3])
r = list(SeqIO.parse(path,"fasta"))
if len(r)> n:
    mafft_cline = MafftCommandline(input=path)
    print(mafft_cline)
    stdout,stderr = mafft_cline()
    with open(out,'w') as handle:
        handle.write(stdout)
else:    
    shutil.copyfile(path,out)
    
    

