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
	with open("gene_ids.csv",'r') as g:
		lines = g.readlines()
		for i in range(len(lines)):
			s = lines[i].strip().split(',')
			species = s[0].replace('.fa','')
			start = int(s[1])
			end = int(s[2])
			od[species] = start
			od_e[species] = end
	#print(od)
	#print(od_e)
else:
	print("Doing unweighted sampling")
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
		config["GENOMES"]+"/{sample}.fa.gz"
	params:
		LENGTH=config["LENGTH"],
		KFAC=lambda wildcards: get_index_s(wildcards.sample),
		KFAC_e=lambda wildcards: get_index_e(wildcards.sample),
		THRES=config["UPPER_CASE"]
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{sample}.sample.txt"
	output:
        	config["OUT_DIR"]+"/samples/{sample}_temp.fa"
	threads:config["CORES"]
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


num_species = len(os.listdir(config["GENOMES"]))
if config["TO_ALIGN"] != num_species:
	g = config["OUT_DIR"]+"/samples/{sample}_genes.fa"
else:
	g = config["OUT_DIR"]+"/samples/out.fa"

rule lastz2fasta:
	input:
		expand(config["OUT_DIR"]+"/alignments/{sample}.maf",sample=SAMPLES)   
	output:
		expand(config["OUT_DIR"]+"/genes/gene_{id}.fa",id=IDS),
		report(config["OUT_DIR"]+"/plots/num_genes.png",caption="../report/num_genes_p.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/statistics/homologues.csv",caption="../report/homologues.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/statistics/num_genes.csv",caption="../report/num_genes_t.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/statistics/num_gt.txt",caption="../report/num_gt.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/plots/gene_dup.png",caption="../report/gene_dup.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/plots/homologues.png",caption="../report/homologues_p.rst",category="Genes Report")
	params:
		k = num,
		out = config["OUT_DIR"]+"/genes",
		p = config["OUT_DIR"]+"/alignments",
		m = config["MIN_ALIGN"],
		plotdir = config["OUT_DIR"]+"/plots",
		statdir = config["OUT_DIR"]+"/statistics",
		tool = config["MSA"],
		d = config["MAX_DUP"]
	shell:
		"python workflow/scripts/lastz2fasta.py -k {params.k} --path {params.p} --outdir {params.out} -m {params.m} --plotdir {params.plotdir} --statdir {params.statdir} -d {params.d} --tool {params.tool}" 


rule lastz:
	input:
		genes = g,
		genome = config["GENOMES"]+"/{sample}.fa.gz"
	output:
		config["OUT_DIR"]+"/alignments/{sample}.maf"
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{sample}.lastz.txt"
	threads: 2
	params:
		species = "{sample}",
		identity = config['IDENTITY'],
		coverage = config['COVERAGE'],
		continuity = config['CONTINUITY'],
		align_dir = config['OUT_DIR']+ "/alignments",
		max_dup = 2*int(config['MAX_DUP']),
		steps = config["STEPS"]
	shell:
		'''
		if [[ "{input.genome}" == *.gz ]]; then
			lastz_32 <(gunzip -dc {input.genome})[multiple] {input.genes} --coverage={params.coverage} --continuity={params.continuity} --filter=identity:{params.identity} --format=maf --output={output} --ambiguous=iupac --step={params.steps} --notransition --queryhspbest={params.max_dup} 
		else
			lastz_32 {input.genome}[multiple] {input.genes} --coverage={params.coverage} --continuity={params.continuity} --filter=identity:{params.identity} --format=maf --output={output} --ambiguous=iupac --step={params.steps} --notransition --queryhspbest={params.max_dup} 
		fi
		'''

rule mashtree:
	input:
		config["OUT_DIR"]+"/genes/gene_{id}.fa"
	output:
		config["OUT_DIR"]+"/genes/gene_{id}.dnd"
	params:
		m=config["MIN_ALIGN"],
		n=config["OUT_DIR"],
		max_len=int(1.5*config["LENGTH"]),
		prefix = "gene_{id}",
		suffix = "fa.aln",
		gene_dir = config["OUT_DIR"] + "/genes/gene_{id}"
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{id}.pasta.txt"
	threads: 8
	shell:
		'''
		if [[ `grep -n '>' {input} | wc -l` -gt {params.m} ]] || [[ `awk 'BEGIN{{l=0;n=0;st=0}}{{if (substr($0,1,1) == ">") {{st=1}} else {{st=2}}; if(st==1) {{n+=1}} else if(st==2) {{l+=length($0)}}}} END{{if (n>0) {{print int((l+n-1)/n)}} else {{print 0}} }}' {input}` -gt {params.max_len} ]]
		then
			mashtree --mindepth 0 --numcpus 32 --kmerlength 10 {params.gene_dir}/*.fa > {output}
		fi
		touch {output}
		'''

rule mergeTrees:
	input:
		expand(config["OUT_DIR"]+"/genes/gene_{id}.dnd",id=IDS)
	output:
		config["OUT_DIR"]+"/genetrees/gene_tree_merged.nwk"
	params:
		msa_dir = config["OUT_DIR"]+"/genes",
	shell:
		'''
		cat {params.msa_dir}/*.dnd > {output}
		'''

rule astral:
	input:
		config["OUT_DIR"]+"/genetrees/gene_tree_merged.nwk"
	output:
		config["OUT_DIR"]+"/roadies.nwk"
	params:
		genes = config["OUT_DIR"]+"/genes",
		stats = config["OUT_DIR"]+"/roadies_stats.nwk"
	benchmark:
		config["OUT_DIR"]+"/benchmarks/astral.txt"
	shell:
		'''
		ASTER-Linux/bin/astral-pro -i {input} -o {output} -a {params.genes}/mapping.txt
		ASTER-Linux/bin/astral-pro -u 3 -i {input} -o {params.stats} -a {params.genes}/mapping.txt
		'''