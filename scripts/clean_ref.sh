#!/bin/bash

# Check for correct number of arguments
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <input_fasta>"
    exit 1
fi

# Input and output files
input_fasta="$1"

# Check if the input file exists
if [[ ! -f "$input_fasta" ]]; then
    echo "Error: File '$input_fasta' does not exist."
    exit 1
fi

# Determine if the input file is gzipped
if [[ "$input_fasta" == *.gz ]]; then
    # Process gzipped file, result file is uncompressed
    gunzip -c "$input_fasta" | awk '/^>/ {print $1; next} {print}' 
else
    # Process plain file
    awk '/^>/ {print $1; next} {print}' "$input_fasta" 
fi

