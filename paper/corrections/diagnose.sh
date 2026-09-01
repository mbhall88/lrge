#!/usr/bin/env bash
# Diagnostics, not part of the user-facing matrix.
#  * -v surfaces the debug line carrying the skew score on runs where auto does NOT engage.
#  * --normalize always separates "the detector missed it" from "normalizing would not have helped".
set -uo pipefail
BENCH=/scratch/user/uqmhal11/lrge-issue34-benchmark
BIN="$BENCH/bin/lrge-254e146"
for acc in "$@"; do
    case "$acc" in SRR8618952) preset=pb ;; *) preset=ont ;; esac
    IN="/scratch/user/uqmhal11/lrge-issue29/${acc}/${acc}.fastq.gz"
    for mode in auto always; do
        d="$BENCH/diagnostics/${acc}/${mode}"; mkdir -p "$d"
        t=$(mktemp -d -p "$d")
        "$BIN" -v -P "$preset" -s 4556 -T 10000 -Q 5000 -t 8 -D "$t" --normalize "$mode" "$IN" \
            > "$d/estimate.out" 2> "$d/lrge.err"
        rm -rf "$t"
        echo "[DIAG] $acc $mode estimate=$(tr -d '\n' < "$d/estimate.out")"
        grep -E "Depth skew|Depth normalization|less than the sum|Using [0-9]+ target" "$d/lrge.err" | sed 's/^/        /'
    done
done
echo "[DIAG-DONE]"
