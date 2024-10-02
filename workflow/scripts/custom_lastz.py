import argparse
import os
import subprocess
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process species and mash distance data.')
    parser.add_argument('species_id_csv', help='CSV file with species IDs.')
    parser.add_argument('mash_distance_csv', help='CSV file with mash distances.')
    parser.add_argument('out_fa', help='FASTA file with gene sequences.')
    parser.add_argument('input_genome_file', help='Directory containing genome sequences.')
    parser.add_argument('output_maf_file', help='Final output maf file.')
    parser.add_argument("--steps", type=int, default=4)
    parser.add_argument("--max_dup", type=int, default=10)
    parser.add_argument("--identity_low", type=int, default=40)
    parser.add_argument("--coverage_low", type=int, default=85)
    parser.add_argument("--continuity_low", type=int, default=85)
    parser.add_argument("--identity_high", type=int, default=65)
    parser.add_argument("--coverage_high", type=int, default=85)
    parser.add_argument("--continuity_high", type=int, default=85)
    parser.add_argument("--num_cpus", type=int, default=1)
    parser.add_argument("--num_instances", type=int, default=4)
    parser.add_argument("--score_file", default="HOXD55.q")
    parser.add_argument("--align_dir", default="output_files/alignments")
    
    return parser.parse_args()

def read_species_ids(species_id_csv):
    species_ids = {}
    with open(species_id_csv, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            species = parts[0]
            ids = parts[1:]  # Get all IDs
            species_ids[species] = ids
    return species_ids

def read_mash_distances(mash_distance_csv):
    mash_distances = pd.read_csv(mash_distance_csv, index_col=0)
    return mash_distances

def segregate_species(mash_distances, species, threshold):
    distances = mash_distances.loc[species]
    group_A = distances[distances > threshold].index.tolist()
    group_B = distances[distances <= threshold].index.tolist()
    return group_A, group_B

def get_gene_ids(species_ids, species_list):
    gene_ids = []
    for species in species_list:
        gene_ids.extend(species_ids.get(species, []))
    return gene_ids

def read_sequences(out_fa):
    sequences = {}
    with open(out_fa, 'r') as f:
        seq_id = ''
        seq = ''
        for line in f:
            if line.startswith('>'):
                if seq_id:
                    sequences[seq_id] = seq
                # seq_id = line.strip()[1:]  # Remove '>'
                seq_id = line.strip()[1:].split('_')[-1]
                seq = ''
            else:
                seq += line.strip()
        if seq_id:
            sequences[seq_id] = seq  # Add the last sequence
    return sequences

def write_fasta(sequences, output_file):
    with open(output_file, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f'>gene_{seq_id}\n{seq}\n')

def run_lastz(input_genome_file, gene_file, output_file, params, group):

    species_name = os.path.basename(input_genome_file)
    if species_name.endswith('.fa.gz'):
        command = f"workflow/scripts/lastz_32 <(gunzip -dc {os.path.join(input_genome_file)})[multiple] " \
            f"{gene_file} --coverage={params.coverage_low if group == 'A' else params.coverage_high} " \
            f"--continuity={params.continuity_low if group == 'A' else params.continuity_high} " \
            f"--filter=identity:{params.identity_low if group == 'A' else params.identity_high} " \
            f"--format=maf --output={output_file} " \
            f"--ambiguous=iupac --step={params.steps} --queryhspbest={params.max_dup} "
    elif species_name.endswith('.fa'):
        command = f"workflow/scripts/lastz_32 {os.path.join(input_genome_file)}[multiple] " \
            f"{gene_file} --coverage={params.coverage_low if group == 'A' else params.coverage_high} " \
            f"--continuity={params.continuity_low if group == 'A' else params.continuity_high} " \
            f"--filter=identity:{params.identity_low if group == 'A' else params.identity_high} " \
            f"--format=maf --output={output_file} " \
            f"--ambiguous=iupac --step={params.steps} --queryhspbest={params.max_dup} "
            
    if group == 'A':
        command += f"--scores={params.score_file} "

    subprocess.run(command, shell=True, executable='/bin/bash')

def merge_maf_files(file1, file2, output_file):
    with open(output_file, 'w') as outfile:
        # Function to write only non-header lines to the output file
        def write_non_header_lines(file):
            with open(file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        outfile.write(line)

        def write_all_lines(file):
            with open(file, 'r') as f:
                for line in f:
                    outfile.write(line)

        # Write contents of the first file excluding headers
        write_all_lines(file1)

        # Write contents of the second file excluding headers
        write_non_header_lines(file2)

def main():
    args = parse_arguments()
    
    species_ids = read_species_ids(args.species_id_csv)
    mash_distances = read_mash_distances(args.mash_distance_csv)
    sequences = read_sequences(args.out_fa)

    # Extract the species name from the input genome file
    # species_name = os.path.basename(args.input_genome_file).replace('.fa', '')
    species_name = os.path.basename(args.input_genome_file)
    if species_name.endswith('.fa.gz'):
        species_name = species_name.replace('.fa.gz', '')
    elif species_name.endswith('.fa'):
        species_name = species_name.replace('.fa', '')

    # Process only the specified species
    group_A, group_B = segregate_species(mash_distances, species_name, threshold=0.1)  # Set your threshold
    
    gene_ids_A = get_gene_ids(species_ids, group_A)
    gene_ids_B = get_gene_ids(species_ids, group_B)

    sequences_A = {id_: sequences[id_] for id_ in gene_ids_A if id_ in sequences}
    sequences_B = {id_: sequences[id_] for id_ in gene_ids_B if id_ in sequences}

    write_fasta(sequences_A, f'group_A_{species_name}.fa')
    write_fasta(sequences_B, f'group_B_{species_name}.fa')

    run_lastz(args.input_genome_file, f'group_A_{species_name}.fa', f'output_A_{species_name}.maf', args, 'A')
    run_lastz(args.input_genome_file, f'group_B_{species_name}.fa', f'output_B_{species_name}.maf', args, 'B')

    if not os.path.exists(args.align_dir):
        os.makedirs(args.align_dir)
    merge_maf_files(f'output_A_{species_name}.maf', f'output_B_{species_name}.maf', args.output_maf_file) 

if __name__ == '__main__':
    main()