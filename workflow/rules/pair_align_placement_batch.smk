def compute_lastz_batch_resources(wildcards, input):
	import os
	genome_path = input.genome
	try:
		size = os.path.getsize(genome_path)
	except OSError:
		size = 0
	if size > 5_000_000_000:
		return {"threads": 64, "mem_mb": 256000}  # extra large genome
	elif size > 999_999_999:
		return {"threads": 16, "mem_mb": 64000}   # large genome
	else:
		# Observed OOM (SLURM State=OUT_OF_MEMORY) at 16000: our largest
		# current genome (630MB compressed) peaked ~14.7GB RSS, too close to
		# a 16GB cap. 32000 gives real margin; still cheap relative to node
		# capacity (~2TB), so no reason to cut it close for genomes <1GB.
		return {"threads": 4, "mem_mb": 32000}    # default (covers all 22 current query genomes, all <1GB compressed)

rule lastz_batch:
	input:
		genes = config["REF_DIR"] + "/samples/out_batch_{batch}.fa",
		genome = config["GENOMES"] + "/{sample}." + ("fa.gz" if EXTENSION[0]=="gz" else "fa")
	output:
		config["OUT_DIR"]+"/alignments/{sample}_batch_{batch}.maf"
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{sample}_batch_{batch}.lastz.txt"
	threads: lambda wildcards, input: compute_lastz_batch_resources(wildcards, input)["threads"]
	resources:
		mem_mb=lambda wildcards, input: compute_lastz_batch_resources(wildcards, input)["mem_mb"]
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

# One checkpoint per batch: waits for every query genome's alignment of that
# batch's loci, then parses all of them into per-locus gene_{id}.fa files.
# Checkpointed (not a plain rule) because pasta's input for a given locus id
# can only be resolved once we know which batch that id falls in and that
# batch's parsing has actually run - see fasta_input_for_id in placement.smk.
checkpoint lastz2fasta_batch:
	input:
		lambda wildcards: expand(config["OUT_DIR"] + "/alignments/{sample}_batch_" + str(wildcards.batch) + ".maf", sample=SAMPLES)
	output:
		touch(config["OUT_DIR"] + "/genes/batch_{batch}.done")
	resources:
		mem_mb=8000
	params:
		batch = "{batch}",
		batch_size = config["BATCH_SIZE"],
		k = num,
		path = config["OUT_DIR"]+"/alignments",
		outdir = config["OUT_DIR"]+"/genes",
		m = MIN_ALIGN,
		plotdir = config["OUT_DIR"]+"/plots",
		statdir = config["OUT_DIR"]+"/statistics",
		d = config["MAX_DUP"],
		mode = mode
	shell:
		'''
		mkdir -p {params.plotdir} {params.statdir} {params.outdir}
		python workflow/scripts/lastz2fasta_batch.py --batch {params.batch} --batch_size {params.batch_size} -k {params.k} --path {params.path} --outdir {params.outdir} -m {params.m} --plotdir {params.plotdir} --statdir {params.statdir} -d {params.d} --tool {params.mode}
		'''
