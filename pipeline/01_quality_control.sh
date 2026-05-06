#!/bin/bash
set -e

# Configs
input_dir="rawData" #where your files are
output_qc="QC"
trimmed_dir="trimmedData"
threads=12

mkdir -p "$output_qc/before" "$output_qc/after" "$trimmed_dir"

echo "---- FastQC of raw reads ----"
fastqc -t "$threads" "$input_dir"/*.fastq.gz -o "$output_qc/before"
multiqc "$input_dir" -o "$output_qc/before"

echo "---- fastp: filtering and trimming ----"
for file1 in "$input_dir"/*R1_001.fastq.gz; do #change to your file names
    [ -e "$file1" ] || continue

    filename_base=$(basename "$file1" _R1_001.fastq.gz)
    file2="$input_dir/${filename_base}_R2_001.fastq.gz"

    if [ -e "$file2" ]; then
        echo "Sample: $filename_base"
        
        fastp -w "$threads" \
            -i "$file1" -I "$file2" \
            -o "$trimmed_dir/${filename_base}_R1_trimmed.fastq.gz" \
            -O "$trimmed_dir/${filename_base}_R2_trimmed.fastq.gz" \
            -q 30 -g -x -c -D --dup_calc_accuracy 1 \
            -h "$output_qc/after/${filename_base}_fastp.html" \
            -j "$output_qc/after/${filename_base}_fastp.json"
    else
        echo "ERROR $filename_base not found..."
    fi
done

echo "---- FastQC Post-Trimming ----"
fastqc -t "$threads" "$trimmed_dir"/*.fastq.gz -o "$output_qc/after"
multiqc "$output_qc" -o "$output_qc/after"

echo "Quality Control complete"