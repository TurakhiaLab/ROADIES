import glob
from collections import OrderedDict
import random,os
from pathlib import Path
import subprocess

num_species = len(os.listdir(config["GENOMES"]))
num_genomes = len(SAMPLES)

def fasta_input_for_id(wildcards):
	batch = (int(wildcards.id) - 1) // config["BATCH_SIZE"] + 1
	_ = checkpoints.lastz2fasta_batch.get(batch=batch)
	return config["OUT_DIR"] + f"/genes/gene_{wildcards.id}.fa"

# Whether a locus gets a real placement.sh run (twilight/epa-ng/raxml-ng, all
# genuinely multi-threaded) or the instant touch/cp fallback is decided by
# whether REF_DIR has a usable reference gene tree for it - REF_DIR is a
# static, already-built backbone that doesn't change during this run, so
# that's knowable at DAG-planning time, before the job runs. Measured via
# full SLURM accounting across ~18700 completed pasta jobs (flat 16cpu/64G
# before this change): real placements average 30.6% CPU (of 16) / 168MB
# peak 512MB; fallback jobs average 5.0% CPU / 1.5MB peak 144MB. Both were
# drastically over-provisioned; fallback especially so (0.002% of the mem
# it was requesting).
def pasta_resources(wildcards):
	ref_gene_tree = f"{config['REF_DIR']}/genes/gene_{wildcards.id}_filtered.fa.aln.raxml.bestTree"
	try:
		has_ref_tree = os.path.getsize(ref_gene_tree) > 0
	except OSError:
		has_ref_tree = False
	if has_ref_tree:
		return {"threads": 8, "mem_mb": 2000}
	return {"threads": 1, "mem_mb": 500}

rule pasta:
	input:
		input_sequence = fasta_input_for_id,
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
		roadies_root = lambda wildcards: os.path.abspath(os.path.join(workflow.basedir, "..")),
		gpu = gpu
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{id}.pasta.txt"
	threads: lambda wildcards: pasta_resources(wildcards)["threads"]
	resources:
		mem_mb=lambda wildcards: pasta_resources(wildcards)["mem_mb"]
	shell:
		'''
		if [[ -s {params.ref_gene_tree} ]]
		then
			if [[ `grep -n '>' {input.input_sequence} | wc -l` -gt 0 ]]
			then
				./workflow/scripts/placement.sh {input.input_sequence} {threads} {params.workdir} {params.ref_msa} {params.ref_gene_tree} {params.ref_model} {params.msa} {output.gene_tree} {params.ref_sequences} {params.roadies_root} {params.gpu}

			else
				cp {params.ref_gene_tree} {output.gene_tree}
			fi
		else
			# Backbone locus has no completed reference gene tree (missing
			# entirely, e.g. from a premature backbone stop, or skipped
			# during backbone sampling) - skip placement for this locus
			# instead of failing the run.
			touch {output.gene_tree}
		fi
		'''

rule mergeTrees:
	input:
		expand(config["OUT_DIR"]+"/genes/gene_{id}.fa.aln.raxml.bestTree",id=IDS)
	output:
		original_list=config["OUT_DIR"]+"/genetrees/original_list.txt",
		merged_list=config["OUT_DIR"]+"/genetrees/gene_tree_merged.nwk",
		mapping=config["OUT_DIR"]+"/genes/mapping.txt"
	resources:
		mem_mb=4000
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
		cat {params.msa_dir}/mapping_batch_*.txt > {output.mapping}
		'''
