import glob
from collections import OrderedDict
import random,os
from pathlib import Path
import subprocess
import csv
from collections import defaultdict

group_map = {}
grouped_species = defaultdict(list)
ungrouped_species = []

with open(config["GROUP_CSV"], newline="") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        species = row["species"]
        group = row["group"].strip()
        group = group if group else None
        group_map[species] = group
        if group:
            grouped_species[group].append(species)
        else:
            ungrouped_species.append(species)

SAMPLES = list(group_map.keys())
n = len(SAMPLES)
c = len(grouped_species)
m = sum(len(splist) for splist in grouped_species.values())

prob_dict = {}
for species in SAMPLES:
    group = group_map[species]
    if group:
        m_i = len(grouped_species[group])
        prob = (m / n) * (1 / c) * (1 / m_i)
    else:
        prob = 1 / n
    prob_dict[species] = prob

total = sum(prob_dict.values())
for s in prob_dict:
    prob_dict[s] /= total

od = OrderedDict([(s, 0) for s in SAMPLES])
sampled_species = random.choices(SAMPLES, weights=[prob_dict[s] for s in SAMPLES], k=num)
for s in sampled_species:
    od[s] += 1

temInt=1

od_e = OrderedDict([(key,0) for key in SAMPLES])

for i in range(n):
	od_e[SAMPLES[i]]=od[SAMPLES[i]]+temInt
	od[SAMPLES[i]]=temInt
	temInt=od_e[SAMPLES[i]]

def get_index_s(wildcards):
	return od[wildcards]

def get_index_e(wildcards):
	if od_e[wildcards] == num:
		return od_e[wildcards]+1
	return od_e[wildcards]

rule sequence_select:
	input:
		config["GENOMES"] + "/{sample}." + ("fa.gz" if EXTENSION[0]=="gz" else "fa")
	params:
		LENGTH=config["LENGTH"],
		KFAC=lambda wildcards: get_index_s(wildcards.sample),
		KFAC_e=lambda wildcards: get_index_e(wildcards.sample),
		THRES=config["UPPER_CASE"]
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{sample}.sample.txt"
	threads: lambda wildcards: int(config['num_threads'])
	output:
        	config["OUT_DIR"]+"/samples/{sample}_temp.fa"
	shell:
			'''
			echo "We are starting to sample {input}"
			echo "./workflow/scripts/sampling/build/sampling -i {input} -o {output} -l {params.LENGTH} -s {params.KFAC} -e {params.KFAC_e} -t {params.THRES}"
			time ./workflow/scripts/sampling/build/sampling -i {input} -o {output} -l {params.LENGTH} -s {params.KFAC} -e {params.KFAC_e} -t {params.THRES}
			'''

rule sequence_merge:
	input:
		expand(config["OUT_DIR"]+"/samples/{sample}_temp.fa", sample=SAMPLES),
	params:
		gene_dir = config["OUT_DIR"]+"/samples",
		plotdir = config["OUT_DIR"]+"/plots",
		statdir = config["OUT_DIR"]+"/statistics"
	output:
        	config["OUT_DIR"]+"/samples/out.fa",
			report(config["OUT_DIR"]+"/plots/sampling.png",caption="../report/sampling.rst",category='Sampling Report')
	shell:
		"python3 workflow/scripts/sequence_merge.py {params.gene_dir} {output} {params.plotdir} {params.statdir}"
