num_species = len(os.listdir(config["GENOMES"]))

if (mode == "placement"):
	g = config["REF_DIR"]+"/samples/out.fa"
else:
	g = config["OUT_DIR"]+"/samples/out.fa"
		
rule lastz:
	input:
		genes = g,
		genome = config["GENOMES"] + "/{sample}." + ("fa.gz" if EXTENSION[0]=="gz" else "fa")
	output:
		config["OUT_DIR"]+"/alignments/{sample}.maf"
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{sample}.lastz.txt"
	threads: 2
	resources:
		mem_mb=20000
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
		scores = lambda wildcards: os.path.join(workflow.basedir, "..", config.get("SCORES", "HOXD55.q"))
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
		d = config["MAX_DUP"],
		mode = mode,
		gpu = gpu
	shell:
		"python workflow/scripts/lastz2fasta.py -k {params.k} --path {params.p} --outdir {params.out} -m {params.m} --plotdir {params.plotdir} --statdir {params.statdir} -d {params.d} --tool {params.mode} --gpu {params.gpu}" 
