from Bio import SeqIO as sio
from Bio.SeqRecord import SeqRecord as seqrec
from Bio.Seq import Seq
import os
import time

from test_gen import *

gtNum = 10
depth = 5
genLen = 100
threadNum = 16

class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        now = time.perf_counter()
        elapsedTime = now - self._start_time
        self._start_time = now
        return elapsedTime

tim = Timer()
tim.start()
st = genTree(depth)
print(f"tree generation : {tim.stop()}")
gt = []
seqs = []

for i in range(gtNum):
    seqs.append(genGen(st, genLen))
print(f"data generation : {tim.stop()}")

for i in range(gtNum):
    fastaGenes = []
    for k in seqs[i]:
        fastaGenes.append(seqrec(Seq(seqs[i][k]),id=k,description=k))
    sio.write(fastaGenes, 'tmp/genes' + str(i) + '.fasta', 'fasta')
print(f"writing fasta files : {tim.stop()}")

# for i in range(gtNum):
#     os.system('mafft tmp/genes' + str(i) + '.fasta > tmp/msa' + str(i) + '.fasta 2> /dev/null')
for i in range(gtNum):
    os.system('echo "mafft tmp/genes' + str(i) + '.fasta > tmp/msa' + str(i) + '.fasta 2> /dev/null" >> tmp/mafft.bash')
os.system('cat tmp/mafft.bash | parallel -j ' + str(threadNum))
print(f"mafft : {tim.stop()}")

# for i in range(gtNum):
#     os.system('FastTree -nt tmp/msa' + str(i) + '.fasta > tree.newick 2> /dev/null')
#     gt.append(read_tree_newick('tree.newick'))
for i in range(gtNum):
    os.system('echo "FastTree -nt tmp/msa' + str(i) + '.fasta > tmp/tree' + str(i) + '.newick 2> /dev/null" >> tmp/fasttree.bash')
os.system('cat tmp/fasttree.bash | parallel -j ' + str(threadNum))
for i in range(gtNum):
    gt.append(read_tree_newick('tmp/tree' + str(i) + '.newick'))
print(f"fast tree : {tim.stop()}")

print(f"gene tree generation : {tim.stop()}")

gtStr = ''
for t in gt:
    gtStr = gtStr + t.newick() + '\n'

gtFile = open('geneTrees.newick', 'w')
gtFile.write(gtStr)
gtFile.close()

print(f"save gene trees : {tim.stop()}")

os.system('astral -i geneTrees.newick -o speciesTree.newick 2> /dev/null')
print(f"ASTRAL : {tim.stop()}")
os.system('java -Djava.library.path=A-pro/ASTRAL-MP/lib -jar A-pro/ASTRAL-MP/astral.1.1.6.jar -i geneTrees.newick -o speciesTree.newick 2> /dev/null')
print(f"ASTRAL-Pro : {tim.stop()}")