rule astral:
	input:
		config["OUT_DIR"]+"/geneTree/gene_tree_merged.newick"
	output:
		config["OUT_DIR"]+"/speciesTree.newick"
	params:
		genes = config["OUT_DIR"]+"/genes"
	shell:
		# "astral -i {input} -o {output}"
		"""
		pwd
		ASTER-Linux/bin/astral-pro -i {input} -o {output} -a {params.genes}/mapping.txt
		"""
