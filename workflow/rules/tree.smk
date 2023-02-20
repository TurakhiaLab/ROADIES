
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
		max_len = int(100*config["LENGTH"]/config["IDENTITY"])
	conda: 
		"../envs/tree.yaml"	
	shell:
		'''
		if [[ `grep -n '>' {input.msa} | wc -l` -gt {params.m} ]] && [[ `awk 'BEGIN{{l=0;n=0;st=0}}{{if (substr($0,1,1) == ">") {{st=1}} else {{st=2}}; if(st==1) {{n+=1}} else if(st==2) {{l+=length($0)}}}} END{{if (n>0) {{print int((l+n-1)/n)}} else {{print 0}} }}' {input.msa}` -lt {params.max_len} ]]
		then
			iqtree -s {input.msa} -m GTR+I+G -redo
		else
			touch {output.gene_tree}
		fi
		'''
