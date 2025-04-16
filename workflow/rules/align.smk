num_species = len(os.listdir(config["GENOMES"]))

g = config["OUT_DIR"]+"/samples/out.fa"
		
rule lastz:
	input:
		genes = g,
		genome = config["GENOMES"] + "/{sample}." + ("fa.gz" if EXTENSION[0]=="gz" else "fa")
	output:
		config["OUT_DIR"]+"/alignments/{sample}.maf"
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{sample}.lastz.txt"
	threads: lambda wildcards: int(config['num_threads'])
	params:
		species = "{sample}",
		identity = config['IDENTITY'],
		identity_deep = config['IDENTITY_DEEP'],
		coverage = config['COVERAGE'],
		continuity = config['CONTINUITY'],
		align_dir = config['OUT_DIR']+ "/alignments",
		max_dup = 2*int(config['MAX_DUP']),
		steps = config["STEPS"],
		deep_mode = str(deep_mode),
		scores = config['SCORES']
	shell:
		'''
		if [[ "{input.genome}" == *.gz ]]; then
			if [[ "{params.deep_mode}" == "True" ]]; then
				lastz_40 <(gunzip -dc {input.genome})[multiple] {input.genes} --coverage={params.coverage} --continuity={params.continuity} --filter=identity:{params.identity_deep} --format=maf --output={output} --ambiguous=iupac --step={params.steps} --queryhspbest={params.max_dup} --scores={params.scores}
			else
				lastz_40 <(gunzip -dc {input.genome})[multiple] {input.genes} --coverage={params.coverage} --continuity={params.continuity} --filter=identity:{params.identity} --format=maf --output={output} --ambiguous=iupac --step={params.steps} --queryhspbest={params.max_dup}
			fi
		else
			if [[ "{params.deep_mode}" == "True" ]]; then
				lastz_40 {input.genome}[multiple] {input.genes}  --coverage={params.coverage} --continuity={params.continuity} --filter=identity:{params.identity_deep} --format=maf --output={output} --ambiguous=iupac --step={params.steps} --queryhspbest={params.max_dup} --scores={params.scores}
			else
				lastz_40 {input.genome}[multiple] {input.genes}  --coverage={params.coverage} --continuity={params.continuity} --filter=identity:{params.identity} --format=maf --output={output} --ambiguous=iupac --step={params.steps} --queryhspbest={params.max_dup} 
			fi
		fi
		'''

rule lastz2fasta:
	input:
		expand(config["OUT_DIR"]+"/alignments/{sample}.maf",sample=SAMPLES)   
	output:
		expand(config["OUT_DIR"]+"/genes/gene_{id}.fa",id=IDS),
		report(config["OUT_DIR"]+"/plots/num_genes.png",caption="../report/num_genes_p.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/statistics/homologs.csv",caption="../report/homologs.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/statistics/num_genes.csv",caption="../report/num_genes_t.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/statistics/num_gt.txt",caption="../report/num_gt.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/plots/gene_dup.png",caption="../report/gene_dup.rst",category="Genes Report"),
		report(config["OUT_DIR"]+"/plots/homologs.png",caption="../report/homologs_p.rst",category="Genes Report")
	params:
		k = num,
		out = config["OUT_DIR"]+"/genes",
		p = config["OUT_DIR"]+"/alignments",
		m = MIN_ALIGN,
		plotdir = config["OUT_DIR"]+"/plots",
		statdir = config["OUT_DIR"]+"/statistics",
		d = config["MAX_DUP"]
	shell:
		"python workflow/scripts/lastz2fasta.py -k {params.k} --path {params.p} --outdir {params.out} -m {params.m} --plotdir {params.plotdir} --statdir {params.statdir} -d {params.d}" 


rule pasta:
	input:
		config["OUT_DIR"]+"/genes/gene_{id}.fa"
	output:
		config["OUT_DIR"]+"/genes/gene_{id}.fa.aln"
	params:
		m=MIN_ALIGN,
		n=config["OUT_DIR"],
		max_len=int(1.5*config["LENGTH"]),
		prefix = "gene_{id}",
		suffix = "fa.aln"
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{id}.pasta.txt"
	threads: lambda wildcards: int(config['num_threads'])
	conda: 
		"../envs/msa.yaml"
	shell:
		'''
		if [[ `grep -n '>' {input} | wc -l` -gt {params.m} ]] || [[ `awk 'BEGIN{{l=0;n=0;st=0}}{{if (substr($0,1,1) == ">") {{st=1}} else {{st=2}}; if(st==1) {{n+=1}} else if(st==2) {{l+=length($0)}}}} END{{if (n>0) {{print int((l+n-1)/n)}} else {{print 0}} }}' {input}` -gt {params.max_len} ]]
		then
			input_file={input}
			output_file={output}
			reference=""
			all_matched=true

			while IFS= read -r line; do
				line=$(echo "$line" | tr '[:lower:]' '[:upper:]')
  				if [[ "$line" != ">"* ]]; then
    				if [ -z "$reference" ]; then
      					reference="$line"
    				elif [ "$line" != "$reference" ]; then
       					all_matched=false
      					break
    				fi
 	 			fi
			done < "$input_file"

			if [ "$all_matched" = true ]; then
				cp "$input_file" "$output_file"
			else
				python pasta/run_pasta.py -i {input} -j {params.prefix} --alignment-suffix={params.suffix} --num-cpus {threads}
			fi
		fi
		touch {output}

		'''
		
rule filtermsa:
	input:
		config["OUT_DIR"]+"/genes/gene_{id}.fa.aln"
	output:
		config["OUT_DIR"]+"/genes/gene_{id}_filtered.fa.aln"
	params:
		m = config["FILTERFRAGMENTS"],
		n = config["MASKSITES"],
	conda:
		"../envs/filtermsa.yaml"
	shell:
		'''
		python pasta/run_seqtools.py -masksitesp {params.n} -filterfragmentsp {params.m} -infile {input} -outfile {output}
			
		'''

