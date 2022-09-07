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
temp=0
od_e = OrderedDict([(key,0) for key in SAMPLES])
for i in range(num_genomes):
    od_e[SAMPLES[i]]=od[SAMPLES[i]]+temp
    od[SAMPLES[i]]=temp
    temp=od_e[SAMPLES[i]]
def get_index_s(wildcards):
    return od[wildcards]

def get_index_e(wildcards):
    return od_e[wildcards]

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
        KFAC=lambda wildcards: get_index_s(wildcards.sample),
        KFAC_e=lambda wildcards: get_index_e(wildcards.sample)
    output:
        temp(config["OUT"]+"/samples/{sample}_temp.fa")
    shell:
        '''
        cd sequence_select
        python new_select.py --input ../{input[0]} -s {params.KFAC} -e {params.KFAC_e} -l {params.LENGTH} --output ../{output} 
        cd ..
        '''

rule clean:
    input:
        dir=config["OUT"]
    shell:
        '''
        cd sequence_select
        rm out.fasta
        rm *.fa
        cd ..
        cd input.dir
        rm *.maf
        cd ..
        '''
