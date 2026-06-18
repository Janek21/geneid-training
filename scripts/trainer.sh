#!/bin/bash

#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --job-name=trainer

## qos determining
#SBATCH --qos=test

#SBATCH --time=180

## Mem+cpu (per array task = per species)
#SBATCH --mem=12G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

set -euo pipefail

# Command file produced by notebook 1. Layout:
#   line(s) before the first species block  -> shared preamble (e.g. mkdir -p ../results/trainedParams)
#   one block per species, each starting with:  echo "Specie: <name>"
# Each block holds that species' singularity train + mv + ln lines.
CMD_FILE="total_training.txt"

# Max array tasks to run concurrently (cluster-friendly throttle). Empty = no cap.
THROTTLE="6"

# ── dispatch: submitted without an array index -> become a job array ──
# A plain `sbatch trainer.sh` lands here (SLURM_ARRAY_TASK_ID unset); we count the
# species blocks and resubmit this same script as a 0..N-1 array, then exit. The
# resubmitted tasks each carry an index and fall through to the work section.
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    nblocks=$(grep -c '^echo "Specie:' "$CMD_FILE")
    if [ "$nblocks" -eq 0 ]; then
        echo "No species blocks found in $CMD_FILE; nothing to train."
        exit 0
    fi
    range="0-$((nblocks - 1))"
    [ -n "$THROTTLE" ] && range="${range}%${THROTTLE}"
    echo "Dispatching training array for $nblocks species (--array=$range)"
    exec sbatch --array="$range" "$0"
fi

#record start
echo ">STARTING at $(date)"

# ── select this task's species block ──
# Start line of every species block, in file order.
mapfile -t starts < <(grep -n '^echo "Specie:' "$CMD_FILE" | cut -d: -f1)
nblocks=${#starts[@]}
idx="$SLURM_ARRAY_TASK_ID"

if [ "$idx" -ge "$nblocks" ]; then
    echo "Array index $idx >= species count $nblocks; nothing to do."
    exit 0
fi

# Preamble = everything before the first species block (folder creation etc.).
preamble_end=$(( ${starts[0]} - 1 ))

# Block spans from its start line to just before the next block (or EOF for the last).
start=${starts[$idx]}
if [ $((idx + 1)) -lt "$nblocks" ]; then
    end=$(( ${starts[$((idx + 1))]} - 1 ))
else
    end=$(wc -l < "$CMD_FILE")
fi

echo ">Task $idx running lines ${start}-${end} of $CMD_FILE"

# Run the shared preamble (idempotent), then this species' block.
{
    [ "$preamble_end" -ge 1 ] && sed -n "1,${preamble_end}p" "$CMD_FILE"
    sed -n "${start},${end}p" "$CMD_FILE"
} | bash

#record memory usage
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}") #transfer to mb
echo ">Peak memory was $peak_mem_mb MegaBytes"

#record end
echo ">ENDING at $(date)"
