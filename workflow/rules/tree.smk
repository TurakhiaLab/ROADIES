
rule mergeTrees:
	input:
		expand(config["OUT_DIR"]+"/geneTree/gene_tree_{id}.newick",id=IDS)
	output:
		config["OUT_DIR"]+"/geneTree/gene_tree_merged.newick"
	conda: 
		"../envs/tree.yaml"
	shell:
		"workflow/scripts/astralWrapper.sh {output} {input}"
rule iqtree:
	input:
		config["OUT_DIR"]+"/msa/gene_aln_{id}.fa"
	output:
		config["OUT_DIR"]+"/geneTree/gene_tree_{id}.newick"
	params:
		logDir = config["OUT_DIR"]+"/geneTree/",
		m = config["MIN_ALIGN"]
	conda: 
		"../envs/tree.yaml"	
	shell:
		"workflow/scripts/iqtWrapper.sh {input} {output} {params.logDir} {params.m}"
