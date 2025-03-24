#!/bin/bash

# Check for correct number of arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <gff_file> <sample_size>"
    exit 1
fi

# Input file and sample size
gff_file="$1"
sample_size="$2"

# Check if the file exists
if [[ ! -f "$gff_file" ]]; then
    echo "Error: File '$gff_file' does not exist."
    exit 1
fi

# Count occurrences of each gene ID (multi-exonic genes have more than one occurrence)
multi_exonic_genes=$(awk -F'\t' '{print $9}' "$gff_file" | sort | uniq -c | awk '$1 > 1 {print $2}')

# Sample as many lines as the value specified
sampled_ids=$(echo "$multi_exonic_genes" | shuf -n "$sample_size" 2>/dev/null)

# If there are fewer multi-exonic genes than requested, adjust
if [[ $(echo "$sampled_ids" | wc -l) -lt "$sample_size" ]]; then
    echo "Warning: Not enough multi-exonic genes to sample $sample_size. Sampling all available."
    sampled_ids="$multi_exonic_genes"
fi

# Create a subset of the input with only the sampled multi-exonic gene IDs
subset=$(echo "$sampled_ids" | while read -r id; do grep -w "$id" "$gff_file"; done)

# Output the result to stdout
echo "$subset"
