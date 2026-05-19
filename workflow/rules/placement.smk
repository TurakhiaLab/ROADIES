import glob
from collections import OrderedDict
import random,os
from pathlib import Path
import subprocess

num_species = len(os.listdir(config["GENOMES"]))
num_genomes = len(SAMPLES)

rule pasta:
	input:
		input_sequence = config["OUT_DIR"]+"/genes/gene_{id}.fa",
	output:
		gene_tree = config["OUT_DIR"]+"/genes/gene_{id}.fa.aln.raxml.bestTree"
	params:
		m=MIN_ALIGN,
		n=config["OUT_DIR"],
		max_len=int(1.5*config["LENGTH"]),
		prefix = "gene_{id}",
		suffix = "fa.aln",
		outdir = config["OUT_DIR"]+"/genes",
		workdir = config["OUT_DIR"]+"/genes/gene_{id}",
		msa = config["OUT_DIR"]+"/genes/gene_{id}.fa.aln",
		ref_msa = config["REF_DIR"]+"/genes/gene_{id}_filtered.fa.aln",
        ref_gene_tree = config["REF_DIR"]+"/genes/gene_{id}_filtered.fa.aln.raxml.bestTree",
        ref_model = config["REF_DIR"]+"/genes/gene_{id}_filtered.fa.aln.raxml.bestModel",
		ref_sequences = config["REF_DIR"]+"/genes/gene_{id}.fa",
		roadies_root = lambda wildcards: os.path.abspath(os.path.join(workflow.basedir, ".."))
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{id}.pasta.txt"
	threads: lambda wildcards: int(16)
	shell:
		'''
		if [[ `grep -n '>' {input.input_sequence} | wc -l` -gt 0 ]]
		then
			if [[ -s {params.ref_gene_tree} ]]
			then
				./workflow/scripts/placement.sh {input.input_sequence} {threads} {params.workdir} {params.ref_msa} {params.ref_gene_tree} {params.ref_model} {params.msa} {output.gene_tree} {params.ref_sequences} {params.roadies_root}

			else
				cp {params.ref_gene_tree} {output.gene_tree}
			fi
		else
			cp {params.ref_gene_tree} {output.gene_tree}
		fi
		'''

rule mergeTrees:
	input:
		expand(config["OUT_DIR"]+"/genes/gene_{id}.fa.aln.raxml.bestTree",id=IDS)
	output:
		original_list=config["OUT_DIR"]+"/genetrees/original_list.txt",
		merged_list=config["OUT_DIR"]+"/genetrees/gene_tree_merged.nwk"
	params:
		msa_dir = config["OUT_DIR"]+"/genes",
		plotdir = config["OUT_DIR"]+"/plots",
		statdir = config["OUT_DIR"]+"/statistics"
	shell:
		'''
		for file in {params.msa_dir}/*.fa.aln.raxml.bestTree; do
            id=$(echo $file | sed 's/.*gene_//;s/.fa.aln.raxml.bestTree//')
            cat $file >> {output.merged_list}
            echo "$id, $(cat $file)" >> {output.original_list}
        done
		'''
