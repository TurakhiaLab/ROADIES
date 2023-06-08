rule astral:
	input:
		config["OUT_DIR"]+"/genetrees/gene_tree_merged.nwk"
	output:
		config["OUT_DIR"]+"/roadies.nwk"
	params:
		genes = config["OUT_DIR"]+"/genes"
	shell:
		'''
		ASTER-Linux/bin/astral-pro -i {input} -o {output} -a {params.genes}/mapping.txt
		
		'''
