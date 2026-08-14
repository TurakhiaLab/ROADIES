num_species = len(os.listdir(config["GENOMES"]))

if (mode == "placement"):
    g = config["REF_DIR"] + "/samples/out.fa"
else:
    g = config["OUT_DIR"] + "/samples/out.fa"

rule kegalign:
    input:
        genes = g,
        genome = config["GENOMES"] + "/{sample}." + ("fa.gz" if EXTENSION[0]=="gz" else "fa")
    output:
        maf = config["OUT_DIR"] + "/alignments/{sample}.maf"
    benchmark:
        config["OUT_DIR"] + "/benchmarks/{sample}.lastz.txt"
    threads: lambda wildcards: int(48)
    resources:
        gpu = gpu
    params:
        align_dir = config["OUT_DIR"] + "/alignments",
        identity = config['IDENTITY'],
        identity_deep = config['IDENTITY_DEEP'],
        coverage = config['COVERAGE'],
        continuity = config['CONTINUITY'],
        max_dup = 2*int(config['MAX_DUP']),
        steps = config["STEPS"],
        deep_mode = str(deep_mode),
        scores_path = lambda wildcards: os.path.join(workflow.basedir, "..", config.get("SCORES", "HOXD55.q")),
        num_gpu = gpu,
        # The shell block below cd's into a per-sample work directory before
        # referencing these, so they must be absolute - Snakemake's
        # {input.genome}/{input.genes}/{output.maf} wildcards are otherwise
        # substituted as the relative paths they're declared with, which
        # would resolve against the wrong directory once inside $sample_workdir.
        genome_abs = lambda wildcards, input: os.path.abspath(input.genome),
        genes_abs = lambda wildcards, input: os.path.abspath(input.genes),
        maf_abs = lambda wildcards, output: os.path.abspath(output.maf)
    conda:
        "../envs/kegalign.yaml"
    shell:
        """
		exec > >(tee {wildcards.sample}_timing.log) 2>&1
        sample_workdir={params.align_dir}/{wildcards.sample}
        mkdir -p $sample_workdir
        cd $sample_workdir
        mkdir -p work
        cd work

        /usr/bin/time faToTwoBit <(gzip -cdfq {params.genome_abs}) ref.2bit
        /usr/bin/time faToTwoBit <(gzip -cdfq {params.genes_abs}) query.2bit

        cd ..

        /usr/bin/time -v kegalign {params.genome_abs} {params.genes_abs} work/ \
            --num_gpu {params.num_gpu} \
            --num_threads {threads} > {wildcards.sample}_lastz-commands.txt

		if [[ "{params.deep_mode}" == "True" ]]; then
			awk '{{
			sub(/ 2> /,
				" --coverage={params.coverage} --continuity={params.continuity} --filter=identity:{params.identity_deep} --ambiguous=iupac --step={params.steps} --queryhspbest={params.max_dup} --scores={params.scores_path} 2> ");
				print
			}}' {wildcards.sample}_lastz-commands.txt \
			> {wildcards.sample}_lastz-commands.final.sh
		else
			awk '{{
			sub(/ 2> /,
				" --coverage={params.coverage} --continuity={params.continuity} --filter=identity:{params.identity} --ambiguous=iupac --step={params.steps} --queryhspbest={params.max_dup} 2> ");
				print
			}}' {wildcards.sample}_lastz-commands.txt \
			> {wildcards.sample}_lastz-commands.final.sh
		fi

        chmod +x {wildcards.sample}_lastz-commands.final.sh

        /usr/bin/time -v parallel --max-procs {threads} \
            < {wildcards.sample}_lastz-commands.final.sh

        (echo "##maf version=1"; cat *.maf-) > {params.maf_abs}

        rm -rf $sample_workdir
        """


rule lastz2fasta:
	input:
		expand(config["OUT_DIR"]+"/alignments/{sample}.maf",sample=SAMPLES)
	output:
		expand(config["OUT_DIR"]+"/genes/gene_{id}.fa",id=IDS),
		config["OUT_DIR"]+"/genes/mapping.txt",
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
