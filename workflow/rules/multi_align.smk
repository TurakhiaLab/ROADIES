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
				run_pasta.py -i {input} -j {params.prefix} --alignment-suffix={params.suffix} --num-cpus {threads}
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
	shell:
		'''
		run_seqtools.py -masksitesp {params.n} -filterfragmentsp {params.m} -infile {input} -outfile {output}
			
		'''

