import glob
from collections import OrderedDict
import random,os
from pathlib import Path
import subprocess
num_species = len(os.listdir(config["GENOMES"]))
if num_species != config["TO_ALIGN"]:
	print("Doing sampling of {0} species of {1} total with weight of {2} from {3}".format(config["TO_ALIGN"],num_species,config["WEIGHTED"],config["GENOMES"]))
	cmd = "python3 workflow/scripts/weighted_s.py -k {0} --weighted {1} --genomes {2} --to_align {3} --num_species {4}".format(num,config["WEIGHTED"],config["GENOMES"],config["TO_ALIGN"],num_species)
	print(cmd)
	subprocess.call(cmd.split(),shell=False)
	od = {}
	od_e = {}
	print("opening gene id")
	with open("gene_ids.csv",'r') as g:
		lines = g.readlines()
		for i in range(len(lines)):
			s = lines[i].strip().split(',')
			species = s[0].replace('.fa','')
			start = int(s[1])
			end = int(s[2])
			od[species] = start
			od_e[species] = end
	print(od)
	print(od_e)
else:
	print("Not sampling")
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
subset_map = {}
subset_file = config["SUBSET"]
subset_dir = config["SUBSET_DIR"]
#if subset file does not exist before run
print("subset",config["SUBSET"])
if config["SUBSET"] == None or config["SUBSET"] == "0" or config["SUBSET"] == 0:
	subset_file = "subsets/subsets.txt"
	subset_dir = "subsets"
	print("No subset file specified in config creating one at subsets/subsets.txt")
	os.system('mkdir -p subsets')
	os.system('touch subsets/subsets.txt')
	map_exist = 0

else:
	map_exist = 1
	print("Reading in mapping",config["SUBSET"])
	with open(config["SUBSET"],'r') as f:
		lines = f.readlines()
		for line in lines:
			s = line.strip().split()
			sub = s[0]
			s2 = s[1].strip().split('/')
			g = s2[len(s2)-1].replace('.fa','')
			#print(g)
			if g in subset_map:
				subset_map[g].append(sub)
			else:
				subset_map[g] = [sub]
	for m in subset_map:
		print(m,len(subset_map[m]))
print("subset file is ",subset_file)
print("subset dir is ", subset_dir)
def does_map_exist():
	return map_exist
def is_mapped(wildcards):
	print("w",wildcards)
	#file = config["GENOMES"]
	if wildcards in subset_map:	
		print(wildcards,len(subset_map[wildcards]))
		return len(subset_map[wildcards])
	else:
		return 0
rule sequence_select:
	input:
		config["GENOMES"]+"/{sample}.fa"
	params:
		LENGTH=config["LENGTH"],
		KFAC=lambda wildcards: get_index_s(wildcards.sample),
		KFAC_e=lambda wildcards: get_index_e(wildcards.sample),
		THRES=config["UPPER_CASE"]
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{sample}.sample.txt"
	output:
		config["OUT_DIR"]+"/samples/{sample}_temp.fa"
	threads:1
	shell:
		'''
		echo "We are starting to sample {input}"
		echo "./sampling/build/sampling -i {input} -o {output} -l {params.LENGTH} -s {params.KFAC} -e {params.KFAC_e} -t {params.THRES}"
		time ./sampling/build/sampling -i {input} -o {output} -l {params.LENGTH} -s {params.KFAC} -e {params.KFAC_e} -t {params.THRES}
		'''

rule sequence_merge:
	input:
		expand(config["OUT_DIR"]+"/samples/{sample}_temp.fa", sample=SAMPLES),
	params:
		gene_dir = config["OUT_DIR"]+"/samples",
		plotdir = config["OUT_DIR"]+"/plots"
	conda: 
		"../envs/plots.yaml"
	output:
		config["OUT_DIR"]+"/samples/out.fa",
		report(config["OUT_DIR"]+"/plots/sampling.png",caption="../report/sampling.rst",category='Sampling Report')
	shell:
		"python3 workflow/scripts/sequence_merge.py {params.gene_dir} {output} {params.plotdir}"
rule assign_genes:
	input:
		genes = config["OUT_DIR"]+"/samples/out.fa",
	output:
		config["OUT_DIR"]+"/samples/{sample}_genes.fa"
	params:
		species = "{sample}"
	shell:
		"python3 workflow/scripts/weighted_out.py {input.genes} {params.species} {output}"
