import glob
SAMPLES = glob_wildcards(config["PATH"]+"/{samples}.fa").samples
print(SAMPLES)
rule all:
	input:expand(config["OUT"]+"/{samples}.maf",samples=SAMPLES)	
rule lastz:
	input:
		"sequence_select/out.fasta",
		config["PATH"]+"/{sample}.fa"
	output:
		config["OUT"]+"/{sample}.maf"
	shell:
		"lastz_32 {input[1]}[multiple] {input[0]}[multiple] --filter=coverage:70 --filter=identity:70 --step=20 --format=maf --output={output}"

rule sequence_select:
    input:
        "sequence_select/index.csv",
        PATH=config["PATH"]
    threads:
        16
    params:
        LENGTH=config["LENGTH"] 
    output:
        "sequence_select/out.fasta"
    shell:
        '''
	echo "{SAMPLES}"
        cd sequence_select
        ./get_seq.sh ../{input[0]} ../{input.PATH} ../{output} {params.LENGTH}
        cd ..
        '''

rule get_seeding:
    input:
        config["PATH"]
    params:
        config["KREG"]
    output:
        "sequence_select/index.csv"
    shell: 
        '''
        cd sequence_select
        python3 seeding.py -k {params[0]} --output ../{output[0]} --path ../{input[0]}
        cd ..
        '''
