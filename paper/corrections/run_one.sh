#!/usr/bin/env bash
# One benchmark cell: <accession> <mode>, where mode is auto|auto_F|never|always.
# Paper settings throughout: seed 4556, 10000 targets, 5000 queries, 8 threads.
set -uo pipefail

acc=$1; mode=$2; threads=${3:-8}

BENCH=/scratch/user/uqmhal11/lrge-issue34-benchmark
BIN="$BENCH/bin/lrge-254e146"
COMMIT=254e1462cfa910f507313fe8b718fcab320e275d
INPUT="/scratch/user/uqmhal11/lrge-issue29/${acc}/${acc}.fastq.gz"
SEED=4556; TARGETS=10000; QUERIES=5000

# SRR8618952 is the PacBio control; every other input in this matrix is ONT.
case "$acc" in
    SRR8618952) preset=pb ;;
    *)          preset=ont ;;
esac

case "$mode" in
    auto)   extra=(--normalize auto) ;;
    auto_F) extra=(--normalize auto -F) ;;
    never)  extra=(--normalize never) ;;
    always) extra=(--normalize always) ;;
    *)      echo "[ERROR] unknown mode $mode" >&2; exit 2 ;;
esac

outdir="$BENCH/runs/${acc}/${mode}"; mkdir -p "$outdir"
test -s "$INPUT" || { echo "[ERROR] missing input $INPUT" >&2; exit 3; }

tmpdir=$(mktemp -d -p "$outdir"); trap 'rm -rf "$tmpdir"' EXIT

cmd=("$BIN" -P "$preset" -s "$SEED" -T "$TARGETS" -Q "$QUERIES" -t "$threads" -D "$tmpdir" "${extra[@]}" "$INPUT")
printf '%s\n' "${cmd[*]}" > "$outdir/cmd.txt"
echo "$COMMIT" > "$outdir/commit.txt"
echo "${SLURM_JOB_ID:-none}" > "$outdir/jobid.txt"

/usr/bin/time -v -o "$outdir/time.txt" "${cmd[@]}" > "$outdir/estimate.out" 2> "$outdir/lrge.err"
status=$?
echo "$status" > "$outdir/exit_status.txt"

echo "[CELL] $acc $mode exit=$status estimate=$(tr -d '\n' < "$outdir/estimate.out")"
