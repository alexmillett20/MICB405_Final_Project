#!/bin/bash
# Author: Brendan Ng
# MICB405 Final Project Star Alignment Script
# Date created: October 30, 2025
# Last updated: October 30, 2025

# This scripts performs STAR alignment on paired-end RNA-seq data.
# It assumes the input FASTQ files are named in the format:
# SAMPLE_1.fastq.gz and SAMPLE_2.fastq.gz for read pairs.
# The fastq files and STAR index directory are stored on the MICB 405 server.

# ChatGPT was used to help generate the structure of this script. Modifications were made to the STAR command parameters to suit project needs.

# ============ CONFIGURATION ============
RAW_DIR="/work/data/raw-data"
OUT_DIR="/work/data/aligned"
STAR_INDEX="/work/data/STARIndex"
THREADS=8   # Adjust based on available cores

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

echo "Starting STAR alignments..."
echo "Input dir: $RAW_DIR"
echo "Output dir: $OUT_DIR"
echo "Using STAR index: $STAR_INDEX"
echo "----------------------------------------"

# Loop over all paired-end read sets
for R1 in "$RAW_DIR"/*_1.fastq.gz; do
    # Get the sample name by removing path and _1.fastq.gz
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"

    # Check that both pairs exist
    if [[ ! -f "$R2" ]]; then
        echo "⚠️  Missing pair for $R1, skipping..."
        continue
    fi

    echo "🔹 Processing sample: $SAMPLE"

    STAR \
        --runThreadN $THREADS \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$OUT_DIR/${SAMPLE}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --quantMode GeneCounts

    echo "✅ Finished $SAMPLE"
    echo "----------------------------------------"
done

echo "🎉 All alignments completed!"

