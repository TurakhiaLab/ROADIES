rule all:
    input:
        "./sequence_select/out.fasta"

rule sequence_select:
    input:
        config["PATH"],
        index="./sequence_select/index.csv"
    output:
        "./sequence_select/out.fasta"
    shell:
        '''
        cd sequence_select
        ./get_seq.sh ../{input.index} ../{input[0]} ../{output}
        cd ..
        '''

rule get_seeding:
    input:
        config["PATH"]
    output:
        "./sequence_select/index.csv"
    shell: 
        '''
        cd sequence_select
        python3 seeding.py --output ../{output[0]} --path ../{input[0]}
        cd ..
        '''