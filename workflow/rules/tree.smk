
rule mergeTrees:
	input:
		expand(config["OUT_DIR"]+"/msa/gene_aln_{id}.fa.treefile",id=IDS)
	output:
		config["OUT_DIR"]+"/geneTree/gene_tree_merged.nwk"
	params:
		msa_dir = config["OUT_DIR"]+"/msa"
	conda: 
		"../envs/tree.yaml"
	shell:
		'''
		cat {params.msa_dir}/*.treefile > {output}
		'''

rule iqtree:
	input:
		msa = config["OUT_DIR"]+"/msa/gene_aln_{id}.fa"
	output:
		gene_tree = config["OUT_DIR"]+"/msa/gene_aln_{id}.fa.treefile"
	params:
		m = config["MIN_ALIGN"],
		l = config["LENGTH"],
		g= config["GAPS"]
	conda: 
		"../envs/tree.yaml"	
	shell:
		'''
		python workflow/scripts/filter_msa.py {input.msa} {params.l} {params.g}
		if [[ `grep -n '>' {input.msa} | wc -l` -gt {params.m} ]]
		then
			iqtree -s {input.msa} -m GTR+I+G -redo
		else
			touch {output.gene_tree}
		fi
		'''
