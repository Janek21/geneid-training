#!/bin/bash

process_file() {
    local input_file="$1"

    # Check if the input file exists
    if [[ ! -f "$input_file" ]]; then
        echo "Error: File '$input_file' does not exist."
        exit 1
    fi

    # Extract lines where the 3rd column is "CDS" and save to tmp.gff
    awk '$3 == "CDS"' "$input_file" > tmp.gff

    # Extract fields and parse the 9th column to get the desired values, saving to output file
    awk -F'\t|;' '{split($9, a, "="); print $1, $2, $3, $4, $5, $6, $7, $8, a[2]}' OFS="\t" tmp.gff 

    # Remove temporary file
    rm tmp.gff
} 

# Check if input argument is provided
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Input file
input_file="$1"

# Determine if the file is gzipped
if [[ "$input_file" == *.gz ]]; then
    ann_name=$(basename "$input_file" .gz)
    gunzip -c "$input_file" > "$ann_name"
    process_file "$ann_name" 
    rm "$ann_name"  # Clean up uncompressed temporary file
else
    process_file "$input_file" 
fi
