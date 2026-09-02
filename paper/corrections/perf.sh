#!/usr/bin/env bash
# Confirm (or refute) that the auto-mode overhead is CountMinSketch cache traffic
# rather than the minimizer computation itself.
set -uo pipefail
BENCH=/scratch/user/uqmhal11/lrge-issue34-benchmark
BIN="$BENCH/bin/lrge-254e146"
IN=/scratch/user/uqmhal11/lrge-issue29/SRR8618952/SRR8618952.fastq.gz
d="$BENCH/perf"; mkdir -p "$d"; cd "$d"

for mode in never auto; do
    t=$(mktemp -d -p "$d")
    perf stat -e task-clock,cycles,instructions,cache-references,cache-misses,LLC-load-misses,LLC-store-misses \
        -o "$d/stat.${mode}.txt" \
        "$BIN" -P pb -s 4556 -T 10000 -Q 5000 -t 8 -D "$t" --normalize "$mode" "$IN" \
        > "$d/${mode}.out" 2> "$d/${mode}.err"
    echo "[PERF-STAT] $mode exit=$? est=$(tr -d '\n' < "$d/${mode}.out")"
    rm -rf "$t"
done

# Where does auto actually spend its time?
t=$(mktemp -d -p "$d")
perf record -F 199 -g --call-graph dwarf -o "$d/auto.data" -- \
    "$BIN" -P pb -s 4556 -T 10000 -Q 5000 -t 8 -D "$t" --normalize auto "$IN" \
    > /dev/null 2> "$d/record.err"
echo "[PERF-RECORD] exit=$?"
rm -rf "$t"
perf report -i "$d/auto.data" --stdio --sort symbol --percent-limit 1 2> /dev/null | head -40 > "$d/report.txt"
echo "[PERF-DONE]"
