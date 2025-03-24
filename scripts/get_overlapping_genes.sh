#!/bin/bash

# Check if input file is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.gff output_overlapping_genes.txt"
    exit 1
fi

# Input and output files
INPUT_GFF="$1"
OUTPUT_OVERLAPS="$2"
GENE_BED="genes.bed"

# Step 1: Extract only protein-coding gene coordinates and convert to BED format
awk -F '\t' '$3 == "gene" && $9 ~ /gene_biotype=protein_coding/ {
    match($9, /ID=([^;]+)/, id);
    gene_id = (id[1] != "" ? id[1] : "unknown");
    print $1, $4-1, $5, gene_id, ".", $7;
}' OFS="\t" "$INPUT_GFF" > "$GENE_BED"

# Step 2: Identify overlapping genes using bedtools
~/miniconda3/bin/bedtools intersect -a "$GENE_BED" -b "$GENE_BED" -wo | awk '$4 != $10 {print $4}' | sort -u > "$OUTPUT_OVERLAPS"

# Cleanup
rm -f "$GENE_BED"

echo "Overlapping protein-coding gene list saved to: $OUTPUT_OVERLAPS"
