#!/usr/bin/env bash
# Q1 evidence: how stable is the skew score when it is computed from only 69 reads?
# Re-seeding redraws the 1% detection sample, so the spread across seeds is the
# sampling noise the threshold has to survive.
set -uo pipefail
BENCH=/scratch/user/uqmhal11/lrge-issue34-benchmark
BIN="$BENCH/bin/lrge-254e146"
acc=${1:-SRR26715166}
IN="/scratch/user/uqmhal11/lrge-issue29/${acc}/${acc}.fastq.gz"
d="$BENCH/seedvar/$acc"; mkdir -p "$d"
for seed in 4556 1 2 3 7 11 42 99 1234 31337 271828 8675309; do
    t=$(mktemp -d -p "$d")
    "$BIN" -v -P ont -s "$seed" -T 10000 -Q 5000 -t 8 -D "$t" --normalize auto "$IN" \
        > "$d/seed${seed}.out" 2> "$d/seed${seed}.err"
    rm -rf "$t"
    line=$(grep -oE 'Depth skew (not )?detected \(99\.9th percentile minimizer count is [0-9.]+x the median; sampled reads: [0-9]+\)' "$d/seed${seed}.err" | head -1)
    printf '[SEED] %s seed=%-8s est=%-9s %s\n' "$acc" "$seed" "$(tr -d '\n' < "$d/seed${seed}.out")" "${line:-NO-DETECTOR-LINE}"
done
echo "[SEEDVAR-DONE] $acc"
