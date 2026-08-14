# lastz2fasta.py takes a directory of .maf alignments for each species, parses, and outputs the alignmnents into k gene fastas
# filtering and stat-gathering steps are concurrently done that are discussed further in README.md
# REQUIREMENTS: Biopython, seaborn, matplotlib
# USAGE: `python workflow/scripts/lastz2fasta.py {args}`
import re
import os, glob
import sys
import argparse
from Bio import SeqIO
from operator import itemgetter
import seaborn as sns
import matplotlib.pyplot as plt
from collections import OrderedDict
from Bio.Seq import Seq
import matplotlib.pyplot as plt

# get arguments
parser = argparse.ArgumentParser(
    description="lastz2fasta.py takes a directory of .maf alignments for each species, parses, and outputs the alignmnents into k gene fastas"
)
parser.add_argument("-k", type=int, default=200, help="Number of genes")
parser.add_argument("--path", default="results/alignments")
parser.add_argument("--outdir", default="results/genes")
parser.add_argument("-m", type=int, default=4)
parser.add_argument("--plotdir", default="results/plots")
parser.add_argument("--statdir", default="results/statistics")
parser.add_argument("-d", type=int, default=100)
parser.add_argument("--tool", default="accurate")
parser.add_argument("--batch", required=True, help="Batch identifier to process (e.g., batch_1)")
parser.add_argument("--batch_size", type=int, default=50, help="Batch size")
args = parser.parse_args()
path = args.path
outdir = args.outdir
plotdir = args.plotdir
statdir = args.statdir
tool = args.tool
batch = args.batch 
batch_size = args.batch_size 
# m = args.m
if tool == "placement":
    m = 1 #args.m
else:
    m = args.m
k = args.k
d = args.d
num_genes = {}
num_homologues = {}
# for i in range(1, k + 1):
#     os.system("touch {0}/gene_{1}.fa".format(outdir, i))

# for batch_start in range(1, k + 1, batch_size):
#     batch_end = min(batch_start + batch_size - 1, k)
#     for i in range(batch_start, batch_end + 1):
#         os.system("touch {0}/gene_{1}.fa".format(outdir, i))

# Determine batch-specific gene range
batch_start = (int(batch) - 1) * batch_size + 1
batch_end = min(batch_start + batch_size - 1, k)

# Create only the files for this batch
for i in range(batch_start, batch_end + 1):
    os.makedirs(outdir, exist_ok=True)
    os.system(f"touch {outdir}/gene_{i}.fa")

# open all lastz alignment outputs
for filename in glob.glob(os.path.join(path, f"*batch_{batch}.maf")):
    with open(os.path.join(os.getcwd(), filename), "r") as f:
        # get species name
        s = filename.split("/")
        name = s[-1]
        species = re.sub(r'_batch_\d+\.maf$', '', name)
        lines = f.readlines()
        # make dict of genes for each species
        genes = {}
        # go through every 4th line due to maf format
        for l in range(15, len(lines)):
            if (l - 15) % 4 == 0:
                # get gene id
                gene_line = lines[l + 1].split()
                gene = gene_line[1]
                gene_s = gene.split("_")
                gene_id = gene_s[1]
                # get score of that alignment
                score_line = lines[l - 1].split()
                score_expr = score_line[1].split("=")
                score = int(score_expr[1])
                # get position in species fasta
                seq_line = lines[l].split()
                position = int(seq_line[2])
                # add to dict of genes
                if gene_id not in genes:
                    genes[gene_id] = [(score, l, position)]

                else:
                    genes[gene_id].append((score, l, position))
        # get number of genes for that species
        num_genes[species] = len(genes)
        aa = list(genes.keys())
        bb = list(genes.values())
        with open(statdir + "/gene_to_species_" + batch + ".csv", "a") as w1:
            for i in range(len(aa)):
                w1.write(species + "\t" + str(aa[i]) + "," + str(bb[i]) + "\n")
        # get number of homologues for that gene (#species)
        for gene in genes:
            if gene in num_homologues:
                num_homologues[gene] += 1
            else:
                num_homologues[gene] = 1
            # make list of genes
            gene_list = genes[gene]
            # skip if no homologue
            # sort homologues by score
            gene_list.sort(reverse=True)
            # initialize with highest score
            max_scores = [0]
            positions = [gene_list[0][2]]
            idx = 1
            # go through gene list
            while idx < len(gene_list):
                # get position of that alignment
                pos = gene_list[idx][2]
                # if it's within 2000 bp of another higher scoring alignment skip that alignment
                tooClose = False
                for p in positions:
                    if abs(p - pos) < (2 * 1000):
                        tooClose = True
                        break
                if not tooClose:
                    positions.append(pos)
                    max_scores.append(idx)
                    tooClose = False
                idx = idx + 1
            # limits the number of alignments
            n = len(max_scores)
            if n > d:
                n = d
            # for highest scoring alignments
            for i in range(n):
                # get line number
                l = gene_list[i][1]
                # get sequence
                orient_line = lines[l + 1].split()
                orientation = orient_line[4]
                seq = ""
                seq_line = lines[l].split()
                seq = seq_line[len(seq_line) - 1]
                seq = seq.replace("-", "")
                genome_pos = str(i)
                if orientation == "-":
                    seq = str(Seq(seq).reverse_complement())
                index = species + "_" + str(i)
                # output to gene fasta
                allowed = ["a", "t", "c", "g", "A", "T", "C", "G"]
                good_seq = True
                for i in range(len(seq)):
                    if seq[i] not in allowed:
                        good_seq = False
                        break
                if not good_seq:
                    continue
                with open(outdir + "/gene_" + gene + ".fa", "a") as w:
                    w.write(">" + index + "\n")
                    w.write(seq + "\n")
                w.close()
                # Per-batch mapping file (species for each placed sequence,
                # needed by astral-pro3's -a flag). Batch-scoped rather than
                # a single shared mapping.txt because multiple batches'
                # lastz2fasta_batch checkpoints can run concurrently as
                # separate Snakemake/sbatch jobs - appending to one shared
                # file from parallel processes risks interleaved/torn
                # writes. mergeTrees concatenates all mapping_batch_*.txt
                # into the final mapping.txt once, after every batch is done.
                with open(outdir + "/mapping_batch_" + batch + ".txt", "a") as w2:
                    w2.write(index + " " + species + "\n")
        w1.close()

if tool == "fast":
    batch_start = (batch - 1) * batch_size + 1
    batch_end = min(batch * batch_size, k)

    for i in range(batch_start, batch_end + 1):
        filename = os.path.join(outdir, f"gene_{i}.fa")
        if not os.path.isfile(filename):
            continue

        with open(filename, "r") as f:
            sequences = f.read().strip().split(">")[1:]

        gene_directory = os.path.join(outdir, f"gene_{i}")
        os.makedirs(gene_directory, exist_ok=True)

        for i, sequence in enumerate(sequences):
            header, *lines = sequence.split("\n")
            sequence_text = "\n".join(lines)
            # Extract the sequence name from the header
            seq_name = header.split()[0]
            # Generate the output filename using the sequence name
            output_filename = os.path.join(gene_directory, f"{seq_name}.fa")
            with open(output_filename, "w") as out:
                out.write(f">{header}\n{sequence_text}\n")
