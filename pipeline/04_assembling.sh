#!/bin/bash
set -euxo pipefail

# ========= CONFIG =========
input_dir="dark_matter"
output_base="assembly_results"
final_dir="final_contigs_filtered"
threads=16
memory=64
min_length=500

mkdir -p "$output_base" "$final_dir"
# ==========================

for file1 in "$input_dir"/*_dark_1.fastq.gz; do
    [ -e "$file1" ] || continue

    name=$(basename "$file1" _dark_1.fastq.gz)
    file2="$input_dir/${name}_dark_2.fastq.gz"
    sample_out="$output_base/${name}_spades"
    
    echo "--- Contigs, Assemble!: $name ---"

    # SPAdes
    spades.py --meta -1 "$file1" -2 "$file2" -k 21,33,55,77 -t "$threads" -m "$memory" -o "$sample_out"

    # Contigs > 500 nt
    if [ -f "$sample_out/scaffolds.fasta" ]; then
        
        grep '>' "$sample_out/scaffolds.fasta" | awk -v len="$min_length" '{
            split($0, a, "_"); 
            if (a[4] >= len) print $1
        }' | sed 's/>//g' > "$sample_out/ids_above_500.txt"

        seqtk subseq "$sample_out/scaffolds.fasta" "$sample_out/ids_above_500.txt" | \
        sed "s/>/>${name}_/g" > "$final_dir/${name}_scaffolds_min500.fasta"

        echo "✅ $name"
        
        # If you would like to save space
        # rm -rf "$sample_out"
    else
        echo "❌ $name"
    fi
done