
rule mafft:
	input:
		config["OUT_DIR"]+"/genes/gene_{id}.fa"
	output:
		config["OUT_DIR"]+"/msa/gene_aln_{id}.fa"
	params:
		m=config["MIN_ALIGN"]
	conda: 
		"../envs/mafft.yaml"
	shell:
		'''
		if [[ `grep -n '>' {input} | wc -l` -gt {params.m} ]]
		then
			#mafft --anysymbol --localpair --retree 2 --maxiterate 10 {input} > {output}
			mafft --anysymbol --auto {input} > {output}
		else
			touch {output}
		fi

		'''
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
		d = config["MAX_DUP"]
	conda:
		"../envs/plots.yaml"
	shell:
		"python workflow/scripts/lastz2fasta.py -k {params.k} --path {params.p} --outdir {params.out} -m {params.m} --plotdir {params.plotdir} --statdir {params.statdir} -d {params.d}" 
		
		
rule lastz:
	input:
		genes = config["OUT_DIR"]+"/samples/out.fa",
		genome = config["GENOMES"]+"/{sample}.fa"
	output:
		config["OUT_DIR"]+"/alignments/{sample}.maf"
	conda:
		"../envs/lastz.yaml"
	params:
		species = "{sample}",
		identity = config['IDENTITY'],
		coverage = config['COVERAGE'],
		align_dir = config['OUT_DIR']+ "/alignments"
	shell:
		'''
		if [[ `stat --printf="%s" {input.genome}` -gt 3000000000 ]]
		then
			echo "File size of {input.genome} is `stat --printf="%s" {input.genome}` which is greater than the lastz limit so subsetting fasta"
			python3 workflow/scripts/get_names.py {input.genome} {params.align_dir} 2
			echo "Aligning {params.align_dir}/{params.species}.0.subset"
			lastz_32 {input.genome}[subset={params.align_dir}/{params.species}.0.subset,multiple] {input.genes}[multiple] --filter=coverage:{params.coverage} --filter=identity:{params.identity} --format=maf --output={output} --ambiguous=iupac
			echo "Aligning {params.align_dir}/{params.species}.1.subset"
			lastz_32 {input.genome}[subset={params.align_dir}/{params.species}.1.subset,multiple] {input.genes}[multiple] --filter=coverage:{params.coverage} --filter=identity:{params.identity} --format=maf --output={output} --ambiguous=iupac
			cat {params.align_dir}/{params.species}.0.maf >> {output}
			cat {params.align_dir}/{params.species}.1.maf >> {output}
		else
			echo "File size of {input.genome} is `stat --printf="%s" {input.genome}` so aligning normally"
			lastz_32 {input.genome}[multiple] {input.genes}[multiple] --filter=coverage:{params.coverage} --filter=identity:{params.identity} --format=maf --output={output} --ambiguous=iupac --step=100
		fi
		'''
