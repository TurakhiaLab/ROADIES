import glob
from collections import OrderedDict
import random

#SAMPLES = glob_wildcards(config["PATH"]+"/{samples}.fa").samples
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

def get_index_s(wildcards):
	return od[wildcards]

def get_index_e(wildcards):
	return od_e[wildcards]
rule sequence_select:
	input:
		config["PATH"]+"/{sample}.fa"
	params:
		LENGTH=config["LENGTH"],
		KFAC=lambda wildcards: get_index_s(wildcards.sample),
		KFAC_e=lambda wildcards: get_index_e(wildcards.sample),
		THRES=config["THRESHOLD"]
	conda:
        	"../envs/bio.yaml"
	output:
        	temp(config["OUT"]+"/samples/{sample}_temp.fa")
	shell:
        	"python3 workflow/scripts/sampling.py --input {input[0]} -s {params.KFAC} -t {params.THRES} -e {params.KFAC_e} -l {params.LENGTH} --output {output}"

rule sequence_merge:
	input:
		expand(config["OUT"]+"/samples/{sample}_temp.fa", sample=SAMPLES),
	params:
		dir = config["OUT"]+"/samples",
		plotdir = config["OUT"]+"/plots"
	conda: 
		"../envs/plots.yaml"
	output:
        	config["OUT"]+"/samples/out.fa",
			report(config["OUT"]+"/plots/sampling.png",caption="../report/sampling.rst",category='Sampling Report')
	shell:
		"python3 workflow/scripts/sequence_merge.py {params.dir} {output} {params.plotdir}"
