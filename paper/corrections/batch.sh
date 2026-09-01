#!/usr/bin/env bash
# Run a list of <accession>:<mode> cells in sequence. Cells are chained with ';' rather
# than '&&' so one bad cell does not take the rest of the batch down with it.
set -uo pipefail
BENCH=/scratch/user/uqmhal11/lrge-issue34-benchmark
for cell in "$@"; do
    acc=${cell%%:*}; mode=${cell##*:}
    bash "$BENCH/run_one.sh" "$acc" "$mode" 8
done
echo "[BATCH-DONE] $*"
