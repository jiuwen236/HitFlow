#!/usr/bin/env bash
# import_shapeit5.sh
# This script ensures phase_common from SHAPEIT5 is available and its directory is on PATH.
# It can be sourced within Snakemake shell commands.

# Check if phase_common is already in PATH
if command -v phase_common >/dev/null 2>&1; then
    echo "[import_shapeit5] phase_common found in PATH. Using system installation."
    return 0 2>/dev/null || exit 0
fi

# Resolve the directory where this script resides
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Same location as the other downloaded binaries: <workflow>/bin/
SHAPEIT_DIR="${SCRIPT_DIR}/../../bin/shapeit5"
PHASE_COMMON="${SHAPEIT_DIR}/phase_common"

# Create directory if it does not exist
mkdir -p "${SHAPEIT_DIR}"

if [ -x "${PHASE_COMMON}" ]; then
    echo "[import_shapeit5] Found existing phase_common at ${PHASE_COMMON}."
else
    echo "[import_shapeit5] phase_common not found, downloading static binary..."
    DOWNLOAD_URL="https://github.com/odelaneau/shapeit5/releases/download/v5.1.1/phase_common_static"
    # Use curl if available, otherwise fallback to wget
    if command -v curl >/dev/null 2>&1; then
        curl -L -o "${PHASE_COMMON}" "${DOWNLOAD_URL}"
    elif command -v wget >/dev/null 2>&1; then
        wget -O "${PHASE_COMMON}" "${DOWNLOAD_URL}"
    else
        echo "Error: neither curl nor wget available to download phase_common." >&2
        return 1 2>/dev/null || exit 1
    fi
    chmod +x "${PHASE_COMMON}"
    echo "[import_shapeit5] Download completed and executable permission set."
fi

# Prepend directory to PATH if not already present
case ":$PATH:" in
    *":${SHAPEIT_DIR}:"*) ;; # already in PATH
    *) export PATH="${SHAPEIT_DIR}:$PATH" ;;
esac

echo "[import_shapeit5] PATH updated to include ${SHAPEIT_DIR}."
