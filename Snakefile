import glob
from collections import OrderedDict
import random
SAMPLES = glob_wildcards(config["PATH"]+"/{samples}.fa").samples
od = OrderedDict([(key,0) for key in SAMPLES])
num=config["KREG"]
num_genomes = len(SAMPLES)
for i in range(num):
    index = random.randint(0,num_genomes-1)
    od[SAMPLES[index]] = od[SAMPLES[index]]+1

def get_index(wildcards):
    return od[wildcards]

rule all:
    input:expand(config["OUT"]+"/alignments/{samples}.maf",samples=SAMPLES)	

rule lastz:
	input:
		config["OUT"]+"/samples/out.fasta",
		config["PATH"]+"/{sample}.fa"
	output:
		config["OUT"]+"/alignments/{sample}.maf"
	shell:
		"lastz_32 {input[1]}[multiple] {input[0]}[multiple] --filter=coverage:70 --filter=identity:70 --step=20 --format=maf --output={output}"

rule sequence_merge:
    input:
        expand(config["OUT"]+"/samples/{sample}_temp.fa", sample=SAMPLES),
	dir = config["OUT"]
    output:
        config["OUT"]+"/samples/out.fasta"
    shell:
       "cat {input.dir}/samples/*.fa >> {output}" 

rule sequence_select:
    input:
        config["PATH"]+"/{sample}.fa"
    params:
        LENGTH=config["LENGTH"],
        KFAC=lambda wildcards: get_index(wildcards.sample)
    output:
        temp(config["OUT"]+"/samples/{sample}_temp.fa")
    shell:
        '''
        cd sequence_select
        python my_select.py --input ../{input[0]} -k {params.KFAC} -l {params.LENGTH} --output ../{output}
        cd ..
        '''
