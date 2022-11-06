
rule mafft:
	input:
		config["OUT"]+"/genes/gene_{id}.fa"
	output:
		config["OUT"]+"/msa/gene_aln_{id}.fa"
	params:
		m=config["MIN_ALIGN"]
	conda: 
		"../envs/mafft.yaml"
	shell:
		"workflow/scripts/mafftWrapper.sh {input} {output} {params.m}"
rule lastz2fasta:
	input:
		expand(config["OUT"]+"/alignments/{sample}.maf",sample=SAMPLES)   
	output:
		expand(config["OUT"]+"/genes/gene_{id}.fa",id=IDS),
		report(config["OUT"]+"/plots/num_genes.png",caption="../report/num_genes_p.rst",category="Genes Report"),
		report(config["OUT"]+"/statistics/homologues.csv",caption="../report/homologues.rst",category="Genes Report"),
		report(config["OUT"]+"/statistics/num_genes.csv",caption="../report/num_genes_t.rst",category="Genes Report"),
		report(config["OUT"]+"/statistics/num_gt.txt",caption="../report/num_gt.rst",category="Genes Report"),
		report(config["OUT"]+"/plots/gene_dup.png",caption="../report/gene_dup.rst",category="Genes Report"),
		report(config["OUT"]+"/plots/homologues.png",caption="../report/homologues_p.rst",category="Genes Report")


	params:
		k = num,
		out = config["OUT"]+"/genes",
		p = config["OUT"]+"/alignments",
		m = config["MIN_ALIGN"],
		plotdir = config["OUT"]+"/plots",
		statdir = config["OUT"]+"/statistics"
	conda:
		"../envs/plots.yaml"
	shell:
		"python workflow/scripts/lastz2fasta.py -k {params.k} --path {params.p} --outdir {params.out} -m {params.m} --plotdir {params.plotdir} --statdir {params.statdir}"
		
		
rule lastz:
	input:
		config["OUT"]+"/samples/out.fa",
		config["PATH"]+"/{sample}.fa"
	output:
		config["OUT"]+"/alignments/{sample}.maf"
	conda:
		"../envs/lastz.yaml"
	params:
		i = config['IDENTITY'],
		c = config['COVERAGE']
	shell:
		"lastz_32 {input[1]}[multiple] {input[0]}[multiple] --filter=coverage:{params.c} --filter=identity:{params.i} --format=maf --output={output}"
