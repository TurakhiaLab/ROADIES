import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
if len(sys.argv) != 3:
    print("Needs path to file:\nUsage: python verify.py [genome] [alignments]")
    exit(0)
if sys.argv[1].endswith((".fasta",".fa")):
    genome = sys.argv[1]
    alignment = sys.argv[2]
    records = list(SeqIO.parse(genome,"fasta"))
    f = open(alignment,'r')
    lines = f.readlines()
    for l in range(15,len(lines)):
        if (l-15)%4 == 0:
            split = lines[l].split()
            seq = split[len(split)-1]
            seq = re.sub('-','',seq)
            chrnum = split[1]
            chrnum = (re.sub('chr','',chrnum))
            if chrnum.isnumeric():
                chrnum = int(chrnum)
                ch = records[chrnum-1].seq[0:len(records[0])-1]
                if ch.find(seq) == -1:
                    print("Unverified sequence alignment is faulty")
                    exit(0)
    print("All alignments are verified for chromosomes with integers names (i.e. chr1)")
else:
    print("Genome is not valid fasta file")
