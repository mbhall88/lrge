#!/usr/bin/env bash
# Isolate the cost of the detection pass on an input where auto does not engage.
# Alternates modes across three repeats so page-cache warmth cannot favour one of them.
set -uo pipefail
BENCH=/scratch/user/uqmhal11/lrge-issue34-benchmark
BIN="$BENCH/bin/lrge-254e146"
acc=SRR8618952
IN="/scratch/user/uqmhal11/lrge-issue29/${acc}/${acc}.fastq.gz"
d="$BENCH/timing"; mkdir -p "$d"
for rep in 1 2 3; do
    for mode in never auto; do
        t=$(mktemp -d -p "$d")
        /usr/bin/time -v -o "$d/${acc}.${mode}.${rep}.time" \
            "$BIN" -P pb -s 4556 -T 10000 -Q 5000 -t 8 -D "$t" --normalize "$mode" "$IN" \
            > "$d/${acc}.${mode}.${rep}.out" 2> "$d/${acc}.${mode}.${rep}.err"
        rm -rf "$t"
        el=$(grep 'Elapsed (wall clock)' "$d/${acc}.${mode}.${rep}.time" | awk '{print $NF}')
        rss=$(grep 'Maximum resident set size' "$d/${acc}.${mode}.${rep}.time" | awk '{print $NF}')
        printf '[TIME] %s rep=%s mode=%-5s wall=%s rss_kb=%s est=%s\n' \
            "$acc" "$rep" "$mode" "$el" "$rss" "$(tr -d '\n' < "$d/${acc}.${mode}.${rep}.out")"
    done
done
echo "[TIME-DONE]"
