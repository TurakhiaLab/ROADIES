rule astral:
	input:
		config["OUT"]+"/geneTree/gene_tree_merged.newick"
	output:
		config["OUT"]+"/speciesTree.newick"
	params:
		p = config["polylimit"],
		genes = config["OUT"]+"/genes"
	shell:
		# "astral -i {input} -o {output}"
		"""
		pwd
		ASTER-LINUX/bin/astral-pro -i {input} -o {output} -a {params.genes}/mapping.txt --polylimit {params.p}
		"""
