#!/bin/bash
set -euo pipefail

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
gpu=${11}

export OMP_NUM_THREADS=$threads

# Fail fast with a clear message instead of a silent "command not found" that
# only surfaces several steps later as an unrelated-looking error (e.g. an
# empty alignment file or a raxml-ng "file not found").
if [[ ! -x "${roadies_root}/TWILIGHT/bin/twilight" ]]; then
	echo "ERROR: ${roadies_root}/TWILIGHT/bin/twilight not found or not executable." >&2
	echo "Placement mode needs TWILIGHT built - run roadies_env.sh again (it builds this automatically) or see the ROADIES_XP guide." >&2
	exit 1
fi
if [[ "$gpu" -gt 0 && ! -x "${roadies_root}/MLIPPER/MLIPPER" ]]; then
	echo "ERROR: ${roadies_root}/MLIPPER/MLIPPER not found or not executable." >&2
	echo "GPU placement needs MLIPPER built - run roadies_env.sh again with CUDA/libpll available, or: bash MLIPPER/install/setup_host.sh" >&2
	exit 1
fi

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

mkdir -p $workDir/iter1_tree_output

if [[ "$gpu" -gt 0 ]]; then
	# GPU placement: MLIPPER commits queries directly onto the reference gene tree
	${roadies_root}/MLIPPER/MLIPPER --tree-alignment $OUT_REF_ITR1 --query-alignment $OUT_QUERY_ITR1 --tree $ref_gene_tree --best-model $ref_model --commit-to-tree $output_gene_trees --local-spr --batch-insert-size 5 --local-spr-radius 4 --local-spr-rounds 1 --gpu-auto
else
	# CPU placement: re-optimize the gene tree under a topological constraint from the reference tree
	raxml-ng --msa $output_msa --model GTR+G+F --threads auto{$threads} --worker 1 --tree-constraint $ref_gene_tree --prefix $workDir/iter1_tree_output/gene_tree --redo --blopt nr_safe
	cp $workDir/iter1_tree_output/gene_tree.raxml.bestTree $output_gene_trees
fi
