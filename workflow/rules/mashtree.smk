import glob
from collections import OrderedDict
import random,os
from pathlib import Path
import subprocess

num_species = len(os.listdir(config["GENOMES"]))

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


g = config["OUT_DIR"]+"/samples/out.fa"

rule lastz:
	input:
		gene = g,
		genome = config["GENOMES"] + "/{sample}." + ("fa.gz" if EXTENSION[0]=="gz" else "fa")
	output:
		config["OUT_DIR"]+"/alignments/{sample}.maf"
	threads: lambda wildcards: int(config['num_threads'])
	params:
		identity_low = config['IDENTITY_LOW'],
		identity_high = config['IDENTITY_HIGH'],
		coverage_low = config['COVERAGE_LOW'],
		coverage_high = config['COVERAGE_HIGH'],
		continuity_low = config['CONTINUITY_LOW'],
		continuity_high = config['CONTINUITY_HIGH'],
		align_dir = config['OUT_DIR']+ "/alignments",
		id_dir = config['OUT_DIR']+ "/statistics/sampling_ids.csv",
		max_dup = 2*int(config['MAX_DUP']),
		steps = config["STEPS"],
		mash_dir = config["OUT_DIR"] + "/mash_distances.txt",
		instances = config['NUM_INSTANCES'],
		score_file = config["SCORES"]
	shell:
		"python workflow/scripts/custom_lastz_fast.py {params.id_dir} {params.mash_dir} {input.gene} {input.genome} {output} --steps {params.steps} --max_dup {params.max_dup} --identity_low {params.identity_low} --coverage_low {params.coverage_low} --continuity_low {params.continuity_low} --identity_high {params.identity_high} --coverage_high {params.coverage_high} --continuity_high {params.continuity_high} --num_cpus {threads} --num_instances {params.instances} --score_file {params.score_file} --align_dir {params.align_dir}"


rule lastz2fasta:
	input:
		expand(config["OUT_DIR"]+"/alignments/{sample}.maf",sample=SAMPLES)   
	output:
		expand(config["OUT_DIR"]+"/genes/gene_{id}.fa",id=IDS),
		report(config["OUT_DIR"]+"/plots/num_genes.png",caption="../report/num_genes_p.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/statistics/homologs.csv",caption="../report/homologs.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/statistics/num_genes.csv",caption="../report/num_genes_t.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/statistics/num_gt.txt",caption="../report/num_gt.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/plots/gene_dup.png",caption="../report/gene_dup.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/plots/homologs.png",caption="../report/homologs_p.rst",category="Genes Report")
	params:
		k = num,
		out = config["OUT_DIR"]+"/genes",
		p = config["OUT_DIR"]+"/alignments",
		m = config["MIN_ALIGN"],
		plotdir = config["OUT_DIR"]+"/plots",
		statdir = config["OUT_DIR"]+"/statistics",
		mode = mode,
		d = config["MAX_DUP"]
	shell:
		"python workflow/scripts/lastz2fasta.py -k {params.k} --path {params.p} --outdir {params.out} -m {params.m} --plotdir {params.plotdir} --statdir {params.statdir} -d {params.d} --tool {params.mode}" 

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
		config["OUT_DIR"]+"/benchmarks/{id}.mashtree.txt"
	threads: lambda wildcards: int(config['num_threads'])
	shell:
		'''
		if [[ `grep -n '>' {input} | wc -l` -gt {params.m} ]] || [[ `awk 'BEGIN{{l=0;n=0;st=0}}{{if (substr($0,1,1) == ">") {{st=1}} else {{st=2}}; if(st==1) {{n+=1}} else if(st==2) {{l+=length($0)}}}} END{{if (n>0) {{print int((l+n-1)/n)}} else {{print 0}} }}' {input}` -gt {params.max_len} ]]
		then
			mashtree --mindepth 0 --numcpus {threads} --kmerlength 10 {params.gene_dir}/*.fa > {output}
		fi
		touch {output}
		'''

rule mergeTrees:
	input:
		expand(config["OUT_DIR"]+"/genes/gene_{id}.dnd",id=IDS)
	output:
		original_list=config["OUT_DIR"]+"/genetrees/original_list.txt",
		merged_list=config["OUT_DIR"]+"/genetrees/gene_tree_merged.nwk"
	params:
		msa_dir = config["OUT_DIR"]+"/genes",
		plotdir = config["OUT_DIR"]+"/plots",
		statdir = config["OUT_DIR"]+"/statistics"
	shell:
		'''
		for file in {params.msa_dir}/*.dnd; do
            id=$(echo $file | sed 's/.*gene_//;s/.dnd//')
            cat $file >> {output.merged_list}
            echo "$id, $(cat $file)" >> {output.original_list}
        done
		'''
