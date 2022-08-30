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
    input:expand(config["OUT"]+"/{samples}.maf",samples=SAMPLES)	

rule lastz:
	input:
		"sequence_select/out.fasta",
		config["PATH"]+"/{sample}.fa"
	output:
		config["OUT"]+"/{sample}.maf"
	shell:
		"lastz_32 {input[1]}[multiple] {input[0]}[multiple] --filter=coverage:70 --filter=identity:70 --step=20 --format=maf --output={output}"

rule sequence_merge:
    input:
        expand("sequence_select/{sample}_temp.fa", sample=SAMPLES)
    output:
        "sequence_select/out.fasta"
    shell:
        '''
        cd sequence_select
        cat *.fa >> ../{output} 
        cd ..
        '''

rule sequence_select:
    input:
        config["PATH"]+"/{sample}.fa"
    threads:
        64
    params:
        LENGTH=config["LENGTH"],
        KFAC=lambda wildcards: get_index(wildcards.sample)
    output:
        temp("sequence_select/{sample}_temp.fa")
    shell:
        '''
        cd sequence_select
        python my_select.py --input ../{input[0]} -k {params.KFAC} -l {params.LENGTH} --output ../{output}
        cd ..
        '''
