#!/bin/bash
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --job-name=trainer
#SBATCH --qos=test
#SBATCH --time=180
#SBATCH --mem=12G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

trainFile="${1:-total_training.txt}"
# geneid training, one array task per species. Command file from notebook 1.
source ../scripts/array_lib.sh
run_array  '^echo "Specie:'
