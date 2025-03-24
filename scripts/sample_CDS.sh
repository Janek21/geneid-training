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

# Extract unique gene IDs from the last column
unique_ids=$(awk -F'\t' '{print $9}' "$gff_file" | sort | uniq)

# Sample as many lines as the value specified
sampled_ids=$(echo "$unique_ids" | shuf -n "$sample_size")

# Create a subset of the input with only the sampled IDs
subset=$(echo "$sampled_ids" | while read -r id; do grep -w "$id" "$gff_file"; done)

# Output the result to stdout
echo "$subset"
