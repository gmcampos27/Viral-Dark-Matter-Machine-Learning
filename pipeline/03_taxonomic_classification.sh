#!/bin/bash
set -euxo pipefail

echo "--- Metavirome Analyses ---"
echo "Started at: $(date)"

# ========= CONFIG =========
dirdb="database/kraken2/"  #Please read the Kraken2 manual to install the best database
THREADS=24
dir="noHost/"
dirout="metadata"
dir_dark="dark_matter"

mkdir -p "$dirout" "$dir_dark"
# ==========================

for file1 in "$dir"/*_host_removed_R1.fastq.gz; do
    [ -e "$file1" ] || continue
    
    name=$(basename "$file1" _host_removed_R1.fastq.gz)
    file2="$dir/${name}_host_removed_R2.fastq.gz"

    echo "--- Sample: $name ---"

    # Kraken2
    # --unclassified-out generates dark matter files
    kraken2 --db "$dirdb" \
        --threads "$THREADS" \
        --paired "$file1" "$file2" \
        --use-names \
        --confidence 0.1 \
        --report "$dirout/${name}_kraken2_report.tsv" \
        --output "$dirout/${name}_kraken2_output.txt" \
        --unclassified-out "$dir_dark/${name}_dark_#.fastq"

    echo "Zipping files..."
    gzip -f "$dir_dark/${name}_dark_1.fastq" "$dir_dark/${name}_dark_2.fastq"

    echo "✅ $name"
done

echo "--- Pipeline finished at: $(date) ---"
echo "You can now see which viruses were at you samples"