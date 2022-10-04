
rule mafft:
	input:
		config["OUT"]+"/genes/gene_{id}.fa"
	output:
		config["OUT"]+"/msa/gene_aln_{id}.fa"
	conda: 
		"../envs/mafft.yaml"
	shell:
		"workflow/scripts/mafftWrapper.sh {input} {output}"
rule lastz2fasta:
	input:
		expand(config["OUT"]+"/alignments/{sample}.maf",sample=SAMPLES)   
	output:
		expand(config["OUT"]+"/genes/gene_{id}.fa",id=IDS),
		report(config["OUT"]+"/plots/num_genes.png",caption="num_genes_p.rst",category="Genes Report"),
		report(config["OUT"]+"/statistics/homologues.csv",caption="homologues.rst",category="Genes Report"),
		report(config["OUT"]+"/statistics/num_genes.csv",caption="num_genes_t.rst",category="Genes Report"),
		report(config["OUT"]+"/statistics/num_gt.txt",caption="num_gt.rst",category="Genes Report")
	params:
		k = num,
		out = config["OUT"]+"/genes",
		p = config["OUT"]+"/alignments",
		m = config["MIN_ALIGN"]
	conda:
		"../envs/plots.yaml"
	shell:
		"python workflow/scripts/lastz2fasta.py -k {params.k} --path {params.p} --outdir {params.out} -m {params.m}"
		
		
rule lastz:
	input:
		config["OUT"]+"/samples/out.fa",
		config["PATH"]+"/{sample}.fa"
	output:
		config["OUT"]+"/alignments/{sample}.maf"
	conda:
		"../envs/lastz.yaml"
	shell:
		"lastz_32 {input[1]}[multiple] {input[0]}[multiple] --filter=coverage:75 --filter=identity:75 --format=maf --output={output}"
