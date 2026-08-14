import glob
from collections import OrderedDict, defaultdict
import random,os
from pathlib import Path
import subprocess
import csv
import argparse
import snakemake.io

# Set up argument parser
parser = argparse.ArgumentParser(description="Precompute sampling values.")
parser.add_argument("--group_csv", help="Path to the CSV file with group information", default=None)
parser.add_argument("num", type=int, help="Number of samples.")
parser.add_argument("genomes", help="Path to genome directory")
args = parser.parse_args()

genome_dir = Path(args.genomes)
SAMPLES = []
EXTENSIONS = set()

for f in genome_dir.glob("*.fa*"):  # Matches .fa, .fa.gz, .fa.bz2, etc.
    name = f.name
    if name.endswith(".fa"):
        sample = name[:-3]
        ext = ""
    elif ".fa." in name:
        sample, ext = name.split(".fa.", 1)
    else:
        continue  # Skip weirdly named files
    SAMPLES.append(sample)
    EXTENSIONS.add(ext)

    # Result
SAMPLES = sorted(set(SAMPLES))
EXTENSION = sorted(EXTENSIONS)  # This will be [''] if only .fa is present

# Read input from the provided group CSV file
group_map = {}
grouped_species = defaultdict(list)
ungrouped_species = []

if args.group_csv and os.path.exists(args.group_csv):
    with open(args.group_csv, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            species = row["species"].strip().replace(".fa","").replace(".fasta","").replace(".fa.gz","")
            group = row["group"].strip()
            group = group if group else None
            group_map[species] = group
            if group:
                grouped_species[group].append(species)
            else:
                ungrouped_species.append(species)

# SELECTED_SAMPLES = list(group_map.keys())

    SELECTED_SAMPLES = [s for s in group_map.keys() if s in SAMPLES]
else:
    print("No group CSV provided or file missing. Selecting all species in genome directory.")
    SELECTED_SAMPLES = SAMPLES.copy()
    group_map = {s: None for s in SELECTED_SAMPLES}
    grouped_species = {}

# SAMPLES = glob_wildcards(config["GENOMES"] + "/{samples}.fa.{extension}").samples
n = len(SELECTED_SAMPLES)
c = len(grouped_species)
m = sum(len(splist) for splist in grouped_species.values()) if grouped_species else 0

prob_dict = {}
for species in SELECTED_SAMPLES:
    group = group_map[species]
    if group and grouped_species:
        m_i = len(grouped_species[group])
        prob = (m / n) * (1 / c) * (1 / m_i)
    else:
        prob = 1 / n
    prob_dict[species] = prob

# Normalize probabilities
total = sum(prob_dict.values())
for s in prob_dict:
    prob_dict[s] /= total

# Sampling logic (adjust num as necessary)
# od = OrderedDict([(s, 0) for s in SELECTED_SAMPLES])
# sampled_species = random.choices(SELECTED_SAMPLES, weights=[prob_dict[s] for s in SELECTED_SAMPLES], k=args.num)
# for s in sampled_species:
#     od[s] += 1

use_group_sampling = bool(args.group_csv and os.path.exists(args.group_csv))

# Sampling logic
od = OrderedDict([(s, 0) for s in SELECTED_SAMPLES])

if use_group_sampling:
    # Probabilistic sampling from group CSV
    sampled_species = random.choices(
        SELECTED_SAMPLES,
        weights=[prob_dict[s] for s in SELECTED_SAMPLES],
        k=args.num
    )
    for s in sampled_species:
        od[s] += 1
else:
    # No group CSV: ensure each species is sampled at least once
    od = OrderedDict([(s, 1) for s in SELECTED_SAMPLES])  # everyone gets 1
    remaining = args.num - len(SELECTED_SAMPLES)
    if remaining > 0:
        extra_samples = random.choices(
            SELECTED_SAMPLES,
            weights=[prob_dict[s] for s in SELECTED_SAMPLES],
            k=remaining
        )
        for s in extra_samples:
            od[s] += 1

# Define od_e
temInt = 1
od_e = OrderedDict([(key, 0) for key in SELECTED_SAMPLES])
for species in SELECTED_SAMPLES:
    count = od[species]  # Number of genes sampled from this species
    if count > 0:
        od_e[species] = temInt + count - 1  # End ID
        od[species] = temInt                # Start ID
        temInt = od_e[species] + 1          # Update for next block
    else:
        od[species] = 0
        od_e[species] = 0

# Save values to text file
with open("sampling_output.txt", "w") as f:
    # Fingerprint of the inputs this sampling plan was computed from, so a
    # later run with different GENOMES/num/group_csv doesn't silently reuse
    # a stale plan computed for a different genome directory.
    f.write("# CONFIG_FINGERPRINT\n")
    f.write(f"{args.genomes}|{args.num}|{args.group_csv or ''}\n")
    f.write("# SAMPLES\n")
    for s in SAMPLES:
        f.write(f"{s}\n")
    f.write("# SELECTED_SAMPLES\n")
    for s in SELECTED_SAMPLES:
        f.write(f"{s}\n")
    f.write("# od\n")
    for k, v in od.items():
        f.write(f"{k},{v}\n")
    f.write("# od_e\n")
    for k, v in od_e.items():
        f.write(f"{k},{v}\n")

print("Precomputing done and values saved in sampling_output.txt")