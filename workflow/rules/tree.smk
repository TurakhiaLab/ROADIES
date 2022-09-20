
rule mergeTrees:
	input:
		expand(config["OUT"]+"/geneTree/gene_tree_{id}.newick",id=IDS)
	output:
		config["OUT"]+"/geneTree/gene_tree_merged.newick"
	conda: 
		"../envs/tree.yaml"
	shell:
		"workflow/scripts/astralWrapper.sh {output} {input}"
rule iqtree:
	input:
		config["OUT"]+"/msa/gene_aln_{id}.fa"
	output:
		config["OUT"]+"/geneTree/gene_tree_{id}.newick"
	params:
		logDir = config["OUT"]+"/geneTree/"
	conda: 
		"../envs/tree.yaml"	
	shell:
		"workflow/scripts/iqtWrapper.sh {input} {output} {params.logDir}"
