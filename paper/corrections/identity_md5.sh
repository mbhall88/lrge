#!/usr/bin/env bash
# Compare selected-read files between the auto and never runs of identity_check.sh.
# lrge nests its intermediates in a randomly-named subdirectory, so match on basename.
set -uo pipefail
BENCH=/scratch/user/uqmhal11/lrge-issue34-benchmark
for acc in "$@"; do
    for base in target.fa query.fa; do
        a=$(find "$BENCH/identity/$acc/auto"  -name "$base" -exec md5sum {} \; 2>/dev/null | cut -d' ' -f1)
        b=$(find "$BENCH/identity/$acc/never" -name "$base" -exec md5sum {} \; 2>/dev/null | cut -d' ' -f1)
        sa=$(find "$BENCH/identity/$acc/auto"  -name "$base" -printf '%s' 2>/dev/null)
        sb=$(find "$BENCH/identity/$acc/never" -name "$base" -printf '%s' 2>/dev/null)
        if [ -z "$a" ] || [ -z "$b" ]; then verdict=MISSING
        elif [ "$a" = "$b" ]; then verdict=IDENTICAL
        else verdict=DIFFERENT; fi
        printf '%-13s %-9s %-10s auto=%s (%s B)  never=%s (%s B)\n' "$acc" "$base" "$verdict" "${a:-none}" "${sa:-0}" "${b:-none}" "${sb:-0}"
    done
done
