import glob
from collections import OrderedDict
import random,os
from pathlib import Path
import subprocess
if config["WEIGHTED"] != 0:
	print("Attemping to do weighted sampling")
	if Path(config["QUARTETS"]).is_file():
		cmd = "python3 workflow/scripts/weighted_s.py -k {0} --quartets {1} --weighted {2} --genomes {3} --id_out {4} --align_out {5} --to_align {6} --species_out {7}".format(config["KREG"],config["QUARTETS"],config["WEIGHTED"],config["GENOMES"],config["GENE_IDS"],config["SPECIES_LISTS"],config["TO_ALIGN"],config["SPECIES_IDS"])
		print(cmd)
		subprocess.call(cmd.split(),shell=False)
		od = {}
		od_e = {}
		with open(config["GENE_IDS"],'r') as g:
			lines = g.readlines()
			for i in range(len(lines)):
				s = lines[i].strip().split(',')
				#$print(s)
				species = s[0].replace('.fa','')
				start = int(s[1])
				end = int(s[2])
				#print(species,start,end)
				od[species] = start
				od_e[species] = end
		print(od)
		print(od_e)
	else:
		print("Can't find freqQuad.csv for weighted sampling exiting (Either incorrect file location in config or need to run ROADIES to get needed file)")
		exit()

	#print(od)
else:
	print("Doing unweighted sampling")
	od = OrderedDict([(key,0) for key in SAMPLES])
	num=config["KREG"]
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

#print(od_e)
def get_index_s(wildcards):
	#print(wildcards,type(wildcards))
	return od[wildcards]

def get_index_e(wildcards):
	if od_e[wildcards] == int(config["KREG"]):
		return od_e[wildcards]+1
	return od_e[wildcards]
rule sequence_select:
	input:
		config["GENOMES"]+"/{sample}.fa"
	params:
		LENGTH=config["LENGTH"],
		KFAC=lambda wildcards: get_index_s(wildcards.sample),
		KFAC_e=lambda wildcards: get_index_e(wildcards.sample),
		THRES=config["UPPER_CASE"]
	#conda:
        	#"../envs/bio.yaml"
	threads:16
	output:
        	config["OUT_DIR"]+"/samples/{sample}_temp.fa"
	shell:
			'''
			echo "We are starting to sample {input}"
			echo "./sampling/build/sampling -i {input} -o {output} -l {params.LENGTH} -s {params.KFAC} -e {params.KFAC_e} -t {params.THRES}"
			time ./sampling/build/sampling -i {input} -o {output} -l {params.LENGTH} -s {params.KFAC} -e {params.KFAC_e} -t {params.THRES}
			'''
        	#"python3 workflow/scripts/sampling.py --input {input[]} -s {params.KFAC} -t {params.THRES} -e {params.KFAC_e} -l {params.LENGTH} --output {output}"

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
		ids = config["SPECIES_IDS"]
	output:
		config["OUT_DIR"]+"/samples/{sample}_genes.fa"
	params:
		species = "{sample}"
	shell:
		"python3 workflow/scripts/weighted_out.py {input.genes} {input.ids} {params.species} {output}"