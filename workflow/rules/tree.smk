
rule raxmlng:
	input:
		msa = config["OUT_DIR"]+"/genes/gene_{id}_filtered.fa.aln"
	output:
		gene_tree = config["OUT_DIR"]+"/genes/gene_{id}_filtered.fa.aln.raxml.bestTree"
	params:
		m = config["MIN_ALIGN"],
		max_len = int(100*config["LENGTH"]/config["IDENTITY"])
	threads: 4
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{id}.iqtree.txt"
	shell:
		'''
		if [[ `grep -n '>' {input.msa} | wc -l` -gt {params.m} ]] && [[ `awk 'BEGIN{{l=0;n=0;st=0}}{{if (substr($0,1,1) == ">") {{st=1}} else {{st=2}}; if(st==1) {{n+=1}} else if(st==2) {{l+=length($0)}}}} END{{if (n>0) {{print int((l+n-1)/n)}} else {{print 0}} }}' {input.msa}` -lt {params.max_len} ]]
		then
			raxml-ng --msa {input.msa} --model GTR+G+F --redo --threads 4
		else
			touch {output.gene_tree}
		fi
		'''

rule mergeTrees:
	input:
		expand(config["OUT_DIR"]+"/genes/gene_{id}_filtered.fa.aln.raxml.bestTree",id=IDS)
	output:
		original_list=config["OUT_DIR"]+"/genetrees/original_list.txt",
		merged_list=config["OUT_DIR"]+"/genetrees/gene_tree_merged.nwk"
	params:
		msa_dir = config["OUT_DIR"]+"/genes",
	shell:
		'''
		for file in {params.msa_dir}/*.fa.aln.raxml.bestTree; do
            id=$(echo $file | sed 's/.*gene_//;s/_filtered.fa.aln.raxml.bestTree//')
            cat $file >> {output.merged_list}
            echo "$id, $(cat $file)" >> {output.original_list}
        done
		'''

