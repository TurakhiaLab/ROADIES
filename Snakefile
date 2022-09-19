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
IDS = list(range(1,num+1))
def get_count(gene):
	f = open(gene,'r')
	lines = f.read_lines()
	count = 0	
	for l in lines:
		if '>' in l:
			count = count+1
	print(count)
	if count > 1:
		return 1
	else:
		return 0
rule all:
	input:
		config["OUT"]+"/speciesTree.newick"
		# config["OUT"]+"/geneTree/gene_tree_merged.newick"
		# expand(config["OUT"]+"/geneTree/gene_tree_{id}.newick",id=IDS)
		# expand(config["OUT"]+"/msa/gene_aln_{id}.fa",id=IDS)
		#expand(config["OUT"]+"/msa/gene_aln_{id}.treefile",id=IDS)
rule astral:
	input:
		config["OUT"]+"/geneTree/gene_tree_merged.newick"
	output:
		config["OUT"]+"/speciesTree.newick"
	shell:
		# "astral -i {input} -o {output}"
		"java -Djava.library.path=A-pro/ASTRAL-MP/lib -jar A-pro/ASTRAL-MP/astral.1.1.6.jar -i {input} -o {output} -a mapping.txt --polylimit 3"
rule mergeTrees:
	input:
		expand(config["OUT"]+"/geneTree/gene_tree_{id}.newick",id=IDS)
	output:
		config["OUT"]+"/geneTree/gene_tree_merged.newick"
	shell:
		"./treeBuilder/astralWrapper.sh {output} {input}"
rule iqtree:
	input:
		config["OUT"]+"/msa/gene_aln_{id}.fa"
	output:
		config["OUT"]+"/geneTree/gene_tree_{id}.newick"
	params:
		logDir = config["OUT"]+"/geneTree/"
	shell:
		"./treeBuilder/iqtWrapper.sh {input} {output} {params.logDir}"
# rule fastTree:
# 	input:
# 		config["OUT"]+"/msa/gene_aln_{id}.fa"
# 	output:
# 		config["OUT"]+"/geneTree/gene_tree_{id}.newick"
# 	shell:
# 		"./treeBuilder/ftWrapper.sh {input} {output}"
rule mafft:
	input:
		config["OUT"]+"/genes/gene_{id}.fa"
	output:
		config["OUT"]+"/msa/gene_aln_{id}.fa"
	conda: 
		"envr.yaml"
	shell:
		"./treeBuilder/mafftWrapper.sh {input} {output}"
rule lastz2fasta:
	input:
		expand(config["OUT"]+"/alignments/{sample}.maf",sample=SAMPLES)   
	output:
		expand(config["OUT"]+"/genes/gene_{id}.fa",id=IDS)
	params:
		k = num,
		out = config["OUT"]+"/genes",
		p = config["OUT"]+"/alignments",
		m = config["MIN_ALIGN"]
	conda:
		"envr.yaml"
	shell:
		"python lastz_align/lastz2fasta.py -k {params.k} --path {params.p} --outdir {params.out} -m {params.m}"
		
		
rule lastz:
	input:
		config["OUT"]+"/samples/out.fa",
		config["PATH"]+"/{sample}.fa"
	output:
		config["OUT"]+"/alignments/{sample}.maf"
	conda:
		"envr.yaml"
	shell:
		"lastz_32 {input[1]}[multiple] {input[0]}[multiple] --filter=coverage:90 --filter=identity:90 --step=20 --format=maf --output={output}"

rule sequence_merge:
	input:
		expand(config["OUT"]+"/samples/{sample}_temp.fa", sample=SAMPLES),
		dir = config["OUT"]
	output:
        	config["OUT"]+"/samples/out.fa"
	shell:
		"cat {input.dir}/samples/*.fa >> {output}" 

rule sequence_select:
	input:
		config["PATH"]+"/{sample}.fa"
	params:
		LENGTH=config["LENGTH"],
		KFAC=lambda wildcards: get_index_s(wildcards.sample),
		KFAC_e=lambda wildcards: get_index_e(wildcards.sample),
		THRES=config["THRESHOLD"]
	conda:
        	"envr.yaml"
	output:
        	temp(config["OUT"]+"/samples/{sample}_temp.fa")
	shell:
		"python sequence_select/new_select.py --input {input[0]} -s {params.KFAC} -t {params.THRES} -e {params.KFAC_e} -l {params.LENGTH} --output {output}" 

rule clean:
	input:
		dir=config["OUT"]
	shell:
		'''
		p = `pwd`
		cd {input.dir}
		cd samples
		rm *.fa
		cd ..
		cd alignments
		rm *.maf
		cd $p
		'''
