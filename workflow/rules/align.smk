
rule mafft:
	input:
		config["OUT"]+"/genes/gene_{id}.fa"
	output:
		config["OUT"]+"/msa/gene_aln_{id}.fa"
	conda: 
		"../envs/align.yaml"
	shell:
		"workflow/scripts/mafftWrapper.sh {input} {output}"
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
		"../envs/align.yaml"
	shell:
		"python workflow/scripts/lastz2fasta.py -k {params.k} --path {params.p} --outdir {params.out} -m {params.m}"
		
		
rule lastz:
	input:
		config["OUT"]+"/samples/out.fa",
		config["PATH"]+"/{sample}.fa"
	output:
		config["OUT"]+"/alignments/{sample}.maf"
	conda:
		"../envs/align.yaml"
	shell:
		"lastz_32 {input[1]}[multiple] {input[0]}[multiple] --filter=coverage:90 --filter=identity:90 --format=maf --output={output}"
