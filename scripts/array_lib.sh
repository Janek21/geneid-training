#!/bin/bash
# Shared SLURM job-array runner for the command files produced by the notebooks.
#
# Every command file has the same shape: a short shared preamble (folder creation,
# variable setup) followed by one self-contained block per species/source. A dedicated
# wrapper in job_commands/ (trainer.sh, predictor.sh, merger.sh, busco.sh) sets its
# SBATCH header, then calls
#   run_array <command_file> <block_start_regex>
# where <block_start_regex> matches the first line of each block. Wrappers run with
# cwd=job_commands/, so they source this file as ../scripts/array_lib.sh.
#
# Behaviour:
#   * `sbatch wrapper.sh` with no array index  -> count blocks, resubmit as a 0..N-1 array
#   * each array task runs the preamble (idempotent) + only its own block, in one shell
#     (so the preamble's mkdir/vars are always in place before the block runs)
#
# THROTTLE caps concurrent array tasks (default 6; set empty for no cap).

run_array() {
    set -euo pipefail
    local cmd_file="$1" block_re="$2"
    local throttle="${THROTTLE:-}"

    [ -f "$cmd_file" ] || { echo "Command file not found: $cmd_file" >&2; exit 1; }

    # ── dispatch: no array index yet -> become an array sized to the block count ──
    if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
        local n; n=$(grep -c "$block_re" "$cmd_file")
        [ "$n" -gt 0 ] || { echo "No blocks in $cmd_file; nothing to do."; exit 0; }
        local range="0-$((n - 1))"
        [ -n "$throttle" ] && range="${range}%${throttle}"
        echo "Dispatching array for $n blocks (--array=$range)"
        exec sbatch --array="$range" "$0" "$cmd_file"
    fi

    echo ">STARTING at $(date)"
    local start_time; start_time=$(date +%s)

    # Start line of every block, in file order.
    local starts; mapfile -t starts < <(grep -n "$block_re" "$cmd_file" | cut -d: -f1)
    local n=${#starts[@]} idx="$SLURM_ARRAY_TASK_ID"
    [ "$idx" -lt "$n" ] || { echo "Index $idx >= block count $n; nothing to do."; exit 0; }

    local preamble_end=$(( ${starts[0]} - 1 ))        # lines before the first block
    local start=${starts[$idx]} end
    if [ $((idx + 1)) -lt "$n" ]; then
        end=$(( ${starts[$((idx + 1))]} - 1 ))
    else
        end=$(wc -l < "$cmd_file")
    fi
    echo ">Task $idx running lines ${start}-${end} of $cmd_file"

    # Preamble (idempotent) + this block, piped to one shell so ordering holds.
    {
        [ "$preamble_end" -ge 1 ] && sed -n "1,${preamble_end}p" "$cmd_file"
        sed -n "${start},${end}p" "$cmd_file"
    } | bash

    # Peak memory + wall time for this task.
    local cgroup_dir peak_mem peak_mem_mb
    cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
    peak_mem=$(cat "/sys/fs/cgroup${cgroup_dir}/memory.peak")
    peak_mem_mb=$(awk "BEGIN {printf \"%.2f\", $peak_mem / 1048576}")
    echo ">Peak memory was $peak_mem_mb MegaBytes"
    echo ">Took $(( ($(date +%s) - start_time) / 60 )) minutes"
    echo ">ENDING at $(date)"
}
