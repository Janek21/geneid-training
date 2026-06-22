#!/bin/bash
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --job-name=busco
#SBATCH --qos=test
#SBATCH --time=90
#SBATCH --mem=12G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# BUSCO evaluation, one array task per species/source. Same block layout for both
# flavours, so one wrapper handles either command file:
#   sbatch busco.sh buscoScoring.txt     # regular predictions (notebook 2)
#   sbatch busco.sh M_buscoScoring.txt   # merged annotations  (notebook 3)
source ../scripts/array_lib.sh
run_array "${1:?usage: sbatch busco.sh <buscoScoring.txt|M_buscoScoring.txt>}" '^agat config --expose'
