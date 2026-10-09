#!/usr/bin/env bash
#
# Usage: ./run_benchmarks.sh <num_parallel_jobs> <joblist_file>
#
# Runs all Julia examples listed in joblist_file (one .jl path per line)
# in parallel, one job per core slot, with line-buffered stdout/stderr
# logging to ./logs/.
#
# Example:
#   ./run_benchmarks.sh 8 joblist.txt

set -euo pipefail

# ---- Parse arguments -------------------------------------------------
if [ $# -lt 2 ]; then
    echo "Usage: $0 <num_parallel_jobs> <joblist_file>" >&2
    exit 1
fi

NJOBS="$1"
JOBLIST="$2"

if ! [[ "$NJOBS" =~ ^[0-9]+$ ]] || [ "$NJOBS" -lt 1 ]; then
    echo "Error: num_parallel_jobs must be a positive integer, got '$NJOBS'" >&2
    exit 1
fi

# ---- Check dependencies ----------------------------------------------
for cmd in parallel stdbuf julia; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: required command '$cmd' not found in PATH" >&2
        exit 1
    fi
done

# ---- Check job list ----------------------------------------------------
if [ ! -f "$JOBLIST" ]; then
    echo "Error: job list '$JOBLIST' not found" >&2
    exit 1
fi

if [ ! -s "$JOBLIST" ]; then
    echo "Error: job list '$JOBLIST' is empty" >&2
    exit 1
fi

# ---- Prepare logs directory (clear if it already exists) --------------
LOGS_DIR="logs"

if [ -d "$LOGS_DIR" ]; then
    echo "Logs directory '$LOGS_DIR' already exists — clearing its contents..."
    rm -rf "${LOGS_DIR:?}"/*
else
    mkdir -p "$LOGS_DIR"
fi

# ---- Run all jobs in parallel ------------------------------------------
echo "Running $(wc -l < "$JOBLIST") job(s) with up to $NJOBS in parallel..."

parallel -j "$NJOBS" \
    'stdbuf -oL -eL julia --threads 1 {} > '"$LOGS_DIR"'/{/.}.out 2> '"$LOGS_DIR"'/{/.}.err' \
    :::: "$JOBLIST"

echo "Done. Logs in '$LOGS_DIR/'."
