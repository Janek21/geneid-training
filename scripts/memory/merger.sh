#!/bin/bash
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --job-name=merger
#SBATCH --qos=test
#SBATCH --time=90
#SBATCH --mem=12G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Annotation merging (AGAT), one array task per species/source. Command file from notebook 3.
source ../scripts/array_lib.sh
run_array mergeAnn.txt '^agat config --expose'
