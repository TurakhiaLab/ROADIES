import glob
from collections import OrderedDict
import random,os
from pathlib import Path
import subprocess

num_species = len(os.listdir(config["GENOMES"]))
od = OrderedDict([(key,0) for key in SAMPLES])
num_genomes = len(SAMPLES)
for i in range(num):
	index = random.randint(0,num_genomes-1)
	od[SAMPLES[index]] = od[SAMPLES[index]]+1
temInt=1

od_e = OrderedDict([(key,0) for key in SAMPLES])

for i in range(num_genomes):
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
		expand(config["OUT_DIR"]+"/samples/{sample}_temp.fa", sample=SAMPLES)
	params:
		gene_dir = config["OUT_DIR"]+"/samples",
		plotdir = config["OUT_DIR"]+"/plots",
		statdir = config["OUT_DIR"]+"/statistics"
	output:
        	config["OUT_DIR"]+"/samples/out.fa"
	shell:
		"python3 workflow/scripts/sequence_merge.py {params.gene_dir} {output} {params.plotdir} {params.statdir}"

