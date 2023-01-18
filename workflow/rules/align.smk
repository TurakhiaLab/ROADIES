
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
		"workflow/scripts/mafftWrapper.sh {input} {output} {params.m}"
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
		config["OUT_DIR"]+"/samples/out.fa",
		config["GENOMES"]+"/{sample}.fa"
	output:
		config["OUT_DIR"]+"/alignments/{sample}.maf"
	conda:
		"../envs/lastz.yaml"
	params:
		s = "{sample}",
		i = config['IDENTITY'],
		c = config['COVERAGE'],
		d = config['OUT_DIR']+ "/alignments"
	shell:
		"./workflow/scripts/lastzwrapper.sh {input[1]} {input[0]} {params.c} {params.i} {output} {params.d} {params.s}" 