#!/bin/bash

# stop_snakemake.sh
# ------------------
# Terminate an active Snakemake workflow launched via run_snakemake.sh
#  1) Sends SIGTERM to running snakemake processes started in the given workdir
#  2) Cancels any Slurm jobs that were submitted by that workflow (detected via log files)
#
# Usage:
#   scripts/stop_snakemake.sh [WORKDIR]
#
#   WORKDIR  Path to the workflow directory (defaults to current directory)
#
# This script requires the following executables to be available in $PATH:
#   - pgrep / pkill (from procps)
#   - squeue, scancel (Slurm commands)
#
# ------------------------------------------------------------

set -euo pipefail

WORKDIR="${1:-$(pwd)}"
LOGDIR="$WORKDIR/logs_cluster"

if [ ! -d "$WORKDIR" ]; then
    echo "[ERROR] Provided WORKDIR does not exist: $WORKDIR" >&2
    exit 1
fi

# ------------------------------------------------------------
# 1. Stop local Snakemake processes
# ------------------------------------------------------------
# Detect snakemake processes that were launched from WORKDIR
# First try to find processes with the exact workdir pattern
SNK_PIDS=$(pgrep -f "snakemake.*$WORKDIR" || true)

# If no exact match, look for run_snakemake.sh processes in this directory
if [ -z "$SNK_PIDS" ]; then
    # Get run_snakemake.sh processes
    RUN_SNK_PIDS=$(pgrep -f "run_snakemake.sh" || true)
    
    # Check if any run_snakemake.sh process is running in our workdir
    for pid in $RUN_SNK_PIDS; do
        if [ -n "$pid" ]; then
            # Get the working directory of the process
            PROC_PWD=$(pwdx "$pid" 2>/dev/null | cut -d: -f2 | xargs || true)
            if [ "$PROC_PWD" = "$WORKDIR" ]; then
                SNK_PIDS="$SNK_PIDS $pid"
            fi
        fi
    done
    
    # Also get any child snakemake python processes
    for pid in $RUN_SNK_PIDS; do
        if [ -n "$pid" ]; then
            PROC_PWD=$(pwdx "$pid" 2>/dev/null | cut -d: -f2 | xargs || true)
            if [ "$PROC_PWD" = "$WORKDIR" ]; then
                # Find child processes
                CHILD_PIDS=$(pgrep -P "$pid" || true)
                SNK_PIDS="$SNK_PIDS $CHILD_PIDS"
            fi
        fi
    done
fi

# Clean up the PID list (remove duplicates and empty values)
SNK_PIDS=$(echo $SNK_PIDS | tr ' ' '\n' | sort -u | grep -v '^$' | tr '\n' ' ')

if [ -n "$SNK_PIDS" ]; then
    echo "[INFO] Terminating Snakemake processes: $SNK_PIDS"
    # First try graceful SIGTERM
    kill -15 $SNK_PIDS 2>/dev/null || true
    
    # Continuously check if processes have stopped, with 30s timeout
    echo "[INFO] Waiting for processes to stop gracefully (timeout: 30s)..."
    TIMEOUT=30
    ELAPSED=0
    ALL_STOPPED=false
    
    while [ $ELAPSED -lt $TIMEOUT ]; do
        STILL_RUNNING=()
        for pid in $SNK_PIDS; do
            if kill -0 "$pid" 2>/dev/null; then
                STILL_RUNNING+=("$pid")
            fi
        done
        
        if [ ${#STILL_RUNNING[@]} -eq 0 ]; then
            echo "[INFO] All processes have stopped gracefully."
            ALL_STOPPED=true
            break
        fi
        
        sleep 1
        ELAPSED=$((ELAPSED + 1))
        
        # Print progress every 10 seconds
        if [ $((ELAPSED % 10)) -eq 0 ]; then
            echo "[INFO] Still waiting... (${ELAPSED}s elapsed, processes still running: ${STILL_RUNNING[*]})"
        fi
    done
    
    # Force kill if still alive after timeout
    if [ "$ALL_STOPPED" = false ]; then
        echo "[WARN] Timeout reached. Force killing remaining processes..."
        for pid in $SNK_PIDS; do
            if kill -0 "$pid" 2>/dev/null; then
                echo "[WARN] Force killing process $pid"
                kill -9 "$pid" 2>/dev/null || true
            fi
        done
    fi
else
    echo "[INFO] No running Snakemake process detected for $WORKDIR"
fi

# ------------------------------------------------------------
# 2. Cancel Slurm jobs belonging to this workflow
# ------------------------------------------------------------

if ! command -v squeue &>/dev/null; then
    echo "'squeue' command not found. Cannot cancel Slurm jobs." >&2
    exit 1
fi

# Gather all active Slurm job IDs for current user
ACTIVE_IDS=$(squeue -u "$USER" -h -o "%A")

# Build list of job IDs to cancel (those with matching log files)
IDS_TO_CANCEL=()

if [ -d "$LOGDIR" ]; then
    for id in $ACTIVE_IDS; do
        if [ -e "$LOGDIR/slurm-${id}.out" ] || [ -e "$LOGDIR/slurm-${id}-error.out" ]; then
            IDS_TO_CANCEL+=("$id")
        fi
    done
fi

if [ ${#IDS_TO_CANCEL[@]} -gt 0 ]; then
    echo "[INFO] Cancelling Slurm jobs: ${IDS_TO_CANCEL[*]}"
    scancel "${IDS_TO_CANCEL[@]}" || true
else
    echo "[INFO] No related Slurm jobs found to cancel."
fi

echo "[INFO] Done." 