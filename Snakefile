rule all:
    input:
        "./sequence_select/out.fasta"
print(config["KREG"])
rule sequence_select:
    input:
        "./sequence_select/index.csv",
        PATH=config["PATH"]
    params:
        LENGTH=config["LENGTH"] 
    output:
        "./sequence_select/out.fasta"
    shell:
        '''
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
        "./sequence_select/index.csv"
    shell: 
        '''
        cd sequence_select
        python3 seeding.py -k {params[0]} --output ../{output[0]} --path ../{input[0]}
        cd ..
        '''