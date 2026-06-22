#!/bin/bash
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --job-name=predictor
#SBATCH --qos=test
#SBATCH --time=90
#SBATCH --mem=12G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# geneid prediction, one array task per species/source. Command file from notebook 2.
source ../scripts/array_lib.sh
run_array total_preds.txt '^time geneid'
