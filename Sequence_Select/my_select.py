from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import random
 
# Arguments passed
k= random.randint(0,10) 
l= random.randint(0,300) 
SeqFile="sample.fa"

if len(sys.argv)>1 and (sys.argv[1][-3:]==".fa" or sys.argv[1][-6:]==".fasta"):
    SeqFile=sys.argv[1]
print("\nName of File:", SeqFile)

for i in range(2, len(sys.argv)):
    if sys.argv[i]=="-l":
        l=int(sys.argv[i+1])
    elif sys.argv[i]=="-k":
        k=int(sys.argv[i+1])
     
print("Number of regions:", k)
print("Region Length :", l)
print("\n")

records = list(SeqIO.parse(SeqFile, "fasta"))
sequence=[]
for i in range(k):
    for j in range(5):
        c_seq=random.randint(0,len(records)-1)
        seq_len=len(records[c_seq].seq)
        start=random.randint(0,seq_len-l)
        frag=records[c_seq].seq[start:(start+l)]
        if "a" not in frag or "N" not in frag or "t" not in frag or "g" not in frag or "c" not in frag:
            break
    record = SeqRecord(frag, "fragment_%i" % (i + 1), "", "")
    sequence.append(record)

SeqIO.write(sequence, "out.fasta", "fasta")

