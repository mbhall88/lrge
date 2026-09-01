#!/usr/bin/env bash
# Matrix item 4: on a control where `auto` does not engage, the reads it selects must be
# byte-identical to `never`. Keeps the temp dirs so the selected FASTA files can be hashed.
set -uo pipefail
BENCH=/scratch/user/uqmhal11/lrge-issue34-benchmark
BIN="$BENCH/bin/lrge-254e146"
SEED=4556; TARGETS=10000; QUERIES=5000

for acc in "$@"; do
    case "$acc" in
        SRR8618952) preset=pb ;;
        *)          preset=ont ;;
    esac
    INPUT="/scratch/user/uqmhal11/lrge-issue29/${acc}/${acc}.fastq.gz"
    for mode in auto never; do
        d="$BENCH/identity/${acc}/${mode}"; rm -rf "$d"; mkdir -p "$d"
        "$BIN" -P "$preset" -s "$SEED" -T "$TARGETS" -Q "$QUERIES" -t 8 \
               -C -D "$d" --normalize "$mode" "$INPUT" \
               > "$d/estimate.out" 2> "$d/lrge.err"
        echo "[IDENT] $acc $mode exit=$? estimate=$(tr -d '\n' < "$d/estimate.out")"
    done
    for f in $(cd "$BENCH/identity/${acc}/auto" && find . -type f \( -name '*.fq' -o -name '*.fa' -o -name '*.fasta' -o -name '*.fastq' \) | sort); do
        a=$(md5sum "$BENCH/identity/${acc}/auto/$f" 2>/dev/null | cut -d' ' -f1)
        b=$(md5sum "$BENCH/identity/${acc}/never/$f" 2>/dev/null | cut -d' ' -f1)
        if [ "$a" = "$b" ] && [ -n "$a" ]; then echo "[MD5-SAME] $acc $f $a"; else echo "[MD5-DIFF] $acc $f auto=$a never=$b"; fi
    done
    ls -l "$BENCH/identity/${acc}/auto" "$BENCH/identity/${acc}/never"
done
echo "[IDENT-DONE]"
