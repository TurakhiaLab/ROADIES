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
		java -D"java.library.path=A-pro/ASTRAL-MP/lib" -jar A-pro/ASTRAL-MP/astral.1.1.6.jar -i {input} -o {output} -a {params.genes}/mapping.txt --polylimit {params.p}
		"""
