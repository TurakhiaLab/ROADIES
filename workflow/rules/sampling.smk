rule sequence_select:
	wildcard_constraints:
		sample="|".join(SELECTED_SAMPLES)
	input:
		genome_file = config["GENOMES"] + "/{sample}." + ("fa.gz" if EXTENSION[0]=="gz" else "fa")
	params:
		LENGTH=config["LENGTH"],
		KFAC=lambda wildcards: od[wildcards.sample],
		KFAC_e=lambda wildcards: od_e[wildcards.sample],
		THRES=config["UPPER_CASE"]
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{sample}.sample.txt"
	output: sample_file = config["OUT_DIR"]+"/samples/{sample}_temp.fa"
	shell:
		'''
		echo "We are starting to sample {input}"
		echo "./workflow/scripts/sampling/build/sampling -i {input.genome_file} -o {output.sample_file} -l {params.LENGTH} -s {params.KFAC} -e {params.KFAC_e} -t {params.THRES}"
		time ./workflow/scripts/sampling/build/sampling -i {input.genome_file} -o {output.sample_file} -l {params.LENGTH} -s {params.KFAC} -e {params.KFAC_e} -t {params.THRES}
		'''

rule sequence_merge:
	input:
		expand(config["OUT_DIR"]+"/samples/{sample}_temp.fa", sample=SELECTED_SAMPLES),
	params:
		gene_dir = config["OUT_DIR"]+"/samples",
		plotdir = config["OUT_DIR"]+"/plots",
		statdir = config["OUT_DIR"]+"/statistics"
	output:
		config["OUT_DIR"] + "/samples/out.fa",
		report(config["OUT_DIR"]+"/plots/sampling.png",caption="../report/sampling.rst",category='Sampling Report')
	shell:
		"python3 workflow/scripts/sequence_merge.py {params.gene_dir} {output} {params.plotdir} {params.statdir}"
