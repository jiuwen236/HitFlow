#!/bin/bash

# Run the Snakemake workflow locally or on Slurm
# Usage: scripts/run_snakemake.sh [--cores 32] [--local] [--slurm [PARTITION]] [--report] [options]
#   --local: Run locally
#   --slurm: Run on Slurm; PARTITION is optional (default: the cluster's default partition)
#   --report: Generate a workflow report
#   options: Additional options to pass to snakemake

# Remember the current Snakemake binary if already on PATH
SNAKEMAKE_CMD="$(command -v snakemake 2>/dev/null || true)"

# Ensure Snakemake is available; if not, load the designated conda environment
if [[ -z "$SNAKEMAKE_CMD" ]]; then
    source ~/.bashrc
    conda activate snakemake
    SNAKEMAKE_CMD="$(command -v snakemake 2>/dev/null || true)"
fi

if [[ -z "$SNAKEMAKE_CMD" ]]; then
    echo "Error: Snakemake executable not found. Please install or activate the appropriate environment." >&2
    exit 1
fi

# Ensure we're using the correct conda even if snakemake is available
CONDA_BASE=$(conda info --base 2>/dev/null || echo "")
if [[ -n "$CONDA_BASE" && -d "$CONDA_BASE/bin" ]]; then
    export PATH="$PATH:$CONDA_BASE/bin"
fi

# Ensure the Snakemake environment has the cluster executor plugin available
PLUGIN_NAME="snakemake-executor-plugin-cluster-generic"
PYTHON_BIN="$(dirname "$SNAKEMAKE_CMD")/python"
if [[ ! -x "$PYTHON_BIN" ]]; then
    PYTHON_BIN="$(dirname "$SNAKEMAKE_CMD")/python3"
fi
if [[ ! -x "$PYTHON_BIN" ]]; then
    PYTHON_BIN="$(command -v python3 2>/dev/null || command -v python 2>/dev/null || true)"
fi

if [[ -n "$PYTHON_BIN" ]]; then
    if ! "$PYTHON_BIN" -m pip show "$PLUGIN_NAME" >/dev/null 2>&1; then
        echo "Installing missing dependency: $PLUGIN_NAME"
        if ! "$PYTHON_BIN" -m pip install --upgrade "$PLUGIN_NAME"; then
            echo "Error: Failed to install $PLUGIN_NAME. Please install it manually." >&2
            exit 1
        fi
    fi
else
    echo "Warning: Could not determine Python interpreter for installing $PLUGIN_NAME; please ensure it is available." >&2
fi

# Guarantee a sane PATH inside Apptainer/Singularity containers even with --cleanenv.
# Users can still override this by exporting SINGULARITYENV_PATH manually.
# if [[ -z "${SINGULARITYENV_PATH:-}" ]]; then
#     export SINGULARITYENV_PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
# fi

# ---------------------------
# Default argument values
# ---------------------------
J=100
CORES=32
RUN_MODE="local"        # "local" or "slurm"
SLURM_PARTITION=""
GENERATE_REPORT=false
EXTRA_ARGS=()

# ---------------------------
# Parse command-line options
# ---------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cores)
            CORES="$2"
            shift 2
            ;;
        --j)
            J="$2"
            shift 2
            ;;
        --local)
            RUN_MODE="local"
            shift
            ;;
        --slurm)
            RUN_MODE="slurm"
            if [[ $# -ge 2 && "$2" != -* ]]; then
                SLURM_PARTITION="$2"
                shift
            fi
            shift
            ;;
        --report)
            GENERATE_REPORT=true
            shift
            ;;
        *)
            EXTRA_ARGS+=("$1")
            shift
            ;;
    esac
done

# ---------------------------
# Build common Snakemake flags
# ---------------------------
SMK_FLAGS=("--rerun-incomplete" "--keep-going")

# ---------------------------
# Launch Snakemake
# ---------------------------
echo "Launching Snakemake (mode: $RUN_MODE, cores: $CORES) ..."

SNAKEMAKE_EXIT_CODE=0
if [[ "$RUN_MODE" == "local" ]]; then
    "$SNAKEMAKE_CMD" \
        --cores "$CORES" \
        "${SMK_FLAGS[@]}" \
        "${EXTRA_ARGS[@]}" || SNAKEMAKE_EXIT_CODE=$?
else
    PARTITION_OPT=""
    [[ -n "$SLURM_PARTITION" ]] && PARTITION_OPT="-p ${SLURM_PARTITION} "
    "$SNAKEMAKE_CMD" \
        --executor cluster-generic \
        --cluster-generic-submit-cmd "sbatch ${PARTITION_OPT}--cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.time} -o logs_cluster/slurm-%j.out --error logs_cluster/slurm-%j.out" \
        -j "$J" --cores "$CORES" \
        --latency-wait 60 --retries 1 \
        "${SMK_FLAGS[@]}" \
        "${EXTRA_ARGS[@]}" || SNAKEMAKE_EXIT_CODE=$?
fi

# ---------------------------
# Generate workflow report if requested
# ---------------------------
if [[ "$GENERATE_REPORT" == true && "$SNAKEMAKE_EXIT_CODE" -eq 0 ]]; then
    mkdir -p logs_cluster
    if "$SNAKEMAKE_CMD" --report logs_cluster/workflow_report.html; then
        echo "Finished, created logs_cluster/workflow_report.html"
    else
        SNAKEMAKE_EXIT_CODE=$?
    fi
fi

exit "$SNAKEMAKE_EXIT_CODE"
