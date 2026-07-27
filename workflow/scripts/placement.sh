#!/bin/bash

seqFile=$1
threads=$2
workDir=$3
ref_msa=$4
ref_gene_tree=$5
ref_model=$6
output_msa=$7
output_gene_trees=$8
refseqFile=$9
roadies_root=${10}

export OMP_NUM_THREADS=$threads

mkdir -p $workDir
mkdir $workDir/iter0_msa_input

mkdir $workDir/iter0_tree_output

touch $workDir/iter0_msa.aln

${roadies_root}/TWILIGHT/bin/twilight -a $ref_msa -i $seqFile -o $workDir/iter0_msa.aln -C $threads --match 40 --mismatch -7 --transition 17 --gap-open -140 --gap-extend -10 --overwrite

TIP_QUERY=$(mktemp)
TIP_REF=$(mktemp)

grep '^>' "$seqFile" | sed 's/^>//' > "$TIP_QUERY"

grep '^>' "$ref_msa" | sed 's/^>//' > "$TIP_REF"

OUT_QUERY="$workDir/iter0_output_msa_from_query.fa"
OUT_REF="$workDir/iter0_output_msa_from_ref.fa"

# Extract sequences from $output_msa
extract_sequences() {
    TIPLIST="$1"
    OUTPUT="$3"
    INPUT="$2"
    awk -v tips="$TIPLIST" 'BEGIN {
        while ((getline < tips) > 0) {
            wanted[$1] = 1
        }
        close(tips)
    }
    /^>/ {
        keep = 0
        header = substr($0, 2)
        if (header in wanted) {
            keep = 1
        }
    }
    {
        if (keep) print
    }' "$INPUT" > "$OUTPUT"
}

# Run extraction
extract_sequences "$TIP_QUERY" "$workDir/iter0_msa.aln" "$OUT_QUERY"
extract_sequences "$TIP_REF" "$workDir/iter0_msa.aln" "$OUT_REF"

epa-ng --ref-msa $workDir/iter0_output_msa_from_ref.fa --tree $ref_gene_tree --query $workDir/iter0_output_msa_from_query.fa --model $ref_model --threads $2 --outdir $workDir/iter0_tree_output --redo #--no-heur

gappa examine graft --jplace-path $workDir/iter0_tree_output/epa_result.jplace --out-dir $workDir/iter0_tree_output --fully-resolve --allow-file-overwriting

cat $seqFile $refseqFile > $workDir/iter1_input.fa

${roadies_root}/TWILIGHT/bin/twilight -t $workDir/iter0_tree_output/epa_result.newick -i $workDir/iter1_input.fa -o $output_msa -C $threads --match 40 --mismatch -7 --transition 17 --gap-open -140 --gap-extend -10 --overwrite

OUT_QUERY_ITR1="$workDir/iter1_output_msa_from_query.fa"
OUT_REF_ITR1="$workDir/iter1_output_msa_from_ref.fa"

# Run extraction
extract_sequences "$TIP_QUERY" "$output_msa" "$OUT_QUERY_ITR1"
extract_sequences "$TIP_REF" "$output_msa" "$OUT_REF_ITR1"

# mkdir $workDir/iter1_tree_output

${roadies_root}/MLIPPER/MLIPPER --tree-alignment $workDir/iter1_output_msa_from_ref.fa --query-alignment $workDir/iter1_output_msa_from_query.fa --tree $ref_gene_tree --best-model $ref_model --commit-to-tree $output_gene_trees --local-spr --batch-insert-size 5 --local-spr-radius 4 --local-spr-rounds 1 --gpu-auto

# cp $workDir/iter1_tree_output/gene_tree.raxml.bestTree $output_gene_trees


# #----------------------------------------------------------------

# #!/bin/bash
# set -euo pipefail
# shopt -s nullglob

# # --------------------------
# # Input arguments
# # --------------------------
# seqFile="$1"
# threads="$2"
# workDir="$3"
# ref_msa="$4"
# ref_gene_tree="$5"
# ref_model="$6"
# output_msa="$7"
# output_gene_trees="$8"
# refseqFile="$9"

# export OMP_NUM_THREADS="$threads"

# # --------------------------
# # Create necessary directories
# # --------------------------
# mkdir -p "$workDir" "$workDir/iter0_msa_input" "$workDir/iter0_tree_output" "$workDir/query_parts"

# # --------------------------
# # Function to extract sequences by tip list
# # --------------------------
# extract_sequences() {
#     local tiplist="$1"
#     local msa="$2"
#     local output="$3"

#     awk -v tips="$tiplist" '
#     BEGIN {
#         while ((getline < tips) > 0) wanted[$1]=1
#         close(tips)
#     }
#     /^>/ {
#         keep=0
#         header=substr($0,2)
#         if(header in wanted) keep=1
#     }
#     { if(keep) print }
#     ' "$msa" > "$output"
# }

# # --------------------------
# # Temporary tip files
# # --------------------------
# TIP_QUERY=$(mktemp)
# TIP_REF=$(mktemp)

# grep '^>' "$seqFile" | sed 's/^>//' > "$TIP_QUERY"
# grep '^>' "$ref_msa" | sed 's/^>//' > "$TIP_REF"

# # --------------------------
# # Initial Twilight alignment
# # --------------------------
# initial_msa="$workDir/iter0_msa.aln"

# twilight -a "$ref_msa" -i "$seqFile" -o "$initial_msa" -C "$threads" --match 40 --mismatch -7 --transition 17 --gap-open -140 --gap-extend -10 --overwrite

# # --------------------------
# # Extract sequences from initial MSA
# # --------------------------
# OUT_QUERY="$workDir/iter0_output_msa_from_query.fa"
# OUT_REF="$workDir/iter0_output_msa_from_ref.fa"

# extract_sequences "$TIP_QUERY" "$initial_msa" "$OUT_QUERY"
# extract_sequences "$TIP_REF" "$initial_msa" "$OUT_REF"

# # --------------------------
# # EPA-ng placement
# # --------------------------
# epa_outdir="$workDir/iter0_tree_output"

# epa-ng --ref-msa "$OUT_REF" --tree "$ref_gene_tree" --query "$OUT_QUERY" --model "$ref_model" --threads "$threads" --outdir "$epa_outdir" --redo

# gappa examine graft --jplace-path "$epa_outdir/epa_result.jplace" --out-dir "$epa_outdir" --fully-resolve

# # --------------------------
# # Prepare next input for Twilight
# # --------------------------
# cat "$seqFile" "$refseqFile" > "$workDir/iter1_input.fa"

# twilight -t "$epa_outdir/epa_result.newick" -i "$workDir/iter1_input.fa" -o "$output_msa" -C "$threads" --match 40 --mismatch -7 --transition 17 --gap-open -140 --gap-extend -10 --overwrite

# # --------------------------
# # Extract sequences from iter1 MSA
# # --------------------------
# OUT_QUERY="$workDir/iter1_output_msa_from_query.fa"
# OUT_REF="$workDir/iter1_output_msa_from_ref.fa"

# extract_sequences "$TIP_QUERY" "$output_msa" "$OUT_QUERY"
# extract_sequences "$TIP_REF" "$output_msa" "$OUT_REF"

# # Split query sequences into separate files
# awk -v workDir="$workDir" '/^>/{
#     f=sprintf("%s/query_parts/%s.fa", workDir, substr($0,2))
# }
# { print > f }' "$OUT_QUERY"

# # --------------------------
# # Initialize backbone
# # --------------------------
# cp "$OUT_REF" "$workDir/backbone_msa.fa"
# cp "$ref_gene_tree" "$workDir/backbone_tree.nwk"

# # --------------------------
# # Iterative placement of queries
# # --------------------------
# for query in "$workDir"/query_parts/*.fa; do
#     base=$(basename "$query" .fa)
#     iterDir="$workDir/${base}_placement"
#     mkdir -p "$iterDir"

#     epa-ng --ref-msa "$workDir/backbone_msa.fa" --tree "$workDir/backbone_tree.nwk" --query "$query" --model "$ref_model" --threads "$threads" --outdir "$iterDir"

#     gappa examine graft --jplace-path "$iterDir/epa_result.jplace" --out-dir "$iterDir" --fully-resolve

#     # Update backbone tree
#     [ -f "$iterDir/epa_result.newick" ] || { echo "EPA result missing for $base"; exit 1; }
#     cp "$iterDir/epa_result.newick" "$workDir/backbone_tree.nwk"

#     # Update backbone MSA
#     cat "$workDir/backbone_msa.fa" "$query" > "$workDir/tmp.fa"
#     mv "$workDir/tmp.fa" "$workDir/backbone_msa.fa"
# done

# # --------------------------
# # Final outputs
# # --------------------------
# cp "$workDir/backbone_tree.nwk" "$output_gene_trees"

# # --------------------------
# # Cleanup temporary files
# # --------------------------
# rm -f "$TIP_QUERY" "$TIP_REF"
