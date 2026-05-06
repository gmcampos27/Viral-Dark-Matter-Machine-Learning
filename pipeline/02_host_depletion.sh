#!/bin/bash
set -e

echo "--- Host Depletion using BWA ---"

dir="trimmedData" 
dirout="noHost/"
ref="genomes/Homo_sapiens/GRCh38_latest_genomic.fna"
threads=16

mkdir -p "$dirout"

# It is necessary to have the Human Genome indexed, it was used GRCh38_latest_genomic.fna, avaliable at https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/
if [ ! -f "${ref}.bwt" ]; then
    echo "indexing reference..."
    bwa index "$ref"
fi

for file1 in "$dir"/*_R1_trimmed.fastq.gz; do 
    if [ -f "$file1" ]; then
        name=$(basename "$file1" _R1_trimmed.fastq.gz)
        file2="${dir}/${name}_R2_trimmed.fastq.gz"
        
        echo "Sample: $name"

        # BWA MEM -> Samtools ( -f 4: unmapped) -> BAM
        bwa mem -t "$threads" "$ref" "$file1" "$file2" | \
        samtools view -@ "$threads" -b -f 4 > "$dirout/${name}_unmapped.bam"

        echo "BAM to FASTQ..."
        samtools fastq -@ "$threads" \
            -1 "$dirout/${name}_host_removed_R1.fastq.gz" \
            -2 "$dirout/${name}_host_removed_R2.fastq.gz" \
            -0 /dev/null -s /dev/null -n \
            "$dirout/${name}_unmapped.bam"

        echo "$name finished."
    fi
done

echo "--- Host Depletion complete ---"