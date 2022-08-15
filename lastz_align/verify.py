import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
f = open("alignments/GCA_008658365.1.fa.maf",'r')
lines = f.readlines()
records = list(SeqIO.parse("../../genomes_top40/GCA_008658365.1.fa","fasta"))
for l in range(15,len(lines)):
    if (l-15)%4 == 0:
        split = lines[l].split()
        seq = split[len(split)-1]
        seq = re.sub('-','',seq)
        chrnum = split[1]
        print(chrnum)
        chrnum = (re.sub('chr','',chrnum))
        if chrnum != "Z":
            chrnum = int(chrnum)
        else:
            chrnum = len(records)
        ch = records[chrnum-1].seq[0:len(records[0])-1]
        print(ch.find(seq))
