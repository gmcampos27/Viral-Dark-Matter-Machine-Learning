#!/bin/bash
set -euo pipefail

echo "--- HMMER (VFam & Pfam) ---"

# ========= CONFIG =========
vfam_db="database/vFam-A_2014.hmm"
pfam_db="database/Pfam-A.hmm"
threads=16
e_value="1e-5"

input_dir="final_contigs_filtered"
output_vfam="results_annotation/vfam"
output_pfam="results_annotation/pfam"
pep_dir="results_annotation/proteins"

mkdir -p "$output_vfam" "$output_pfam" "$pep_dir"
# ==========================

for fasta in "$input_dir"/*_scaffolds_min500.fasta; do
    [ -e "$fasta" ] || continue

    name=$(basename "$fasta" _scaffolds_min500.fasta)
    echo "Sample: $name"

    # 1. Translation (6 frames)
    pep="$pep_dir/${name}_pep.fasta"
    transeq -sequence "$fasta" -outseq "$pep" -clean -frame 6

    pep_clean="$pep_dir/${name}_pep_clean.fasta"
    sed 's/_[1-6]$//; s/ /_/g' "$pep" > "$pep_clean"

    # 2. HMMSEARCH VFam
    echo "   > VFam..."
    hmmsearch --cpu "$threads" \
      --domtblout "$output_vfam/${name}_vfam.domtbl" \
      --noali \
      -E "$e_value" \
      "$vfam_db" "$pep_clean" > /dev/null

    # 3. HMMSEARCH Pfam
    echo "   > Pfam..."
    hmmsearch --cpu "$threads" \
      --domtblout "$output_pfam/${name}_pfam.domtbl" \
      --noali \
      -E "$e_value" \
      "$pfam_db" "$pep_clean" > /dev/null

    # If you want to save space
    #rm "$pep" "$pep_clean"

    echo "✅ $name"
done

echo "--- Pipeline finished! ---"