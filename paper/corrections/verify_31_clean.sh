#!/usr/bin/env bash
# Issue #31 verification, built from a tree containing ONLY the #31 change.
set -uo pipefail
cd "/scratch/temp/27749624/claude-1100000562/-scratch-user-uqmhal11-lrge/01a9cc44-050e-460c-9be2-8cb707278dae/scratchpad/lrge-31" && cargo build --release 2>&1 | tail -2
LRGE="/scratch/temp/27749624/claude-1100000562/-scratch-user-uqmhal11-lrge/01a9cc44-050e-460c-9be2-8cb707278dae/scratchpad/lrge-31/target/release/lrge"
run () {
  local acc=$1 tag=$2 plat=$3; shift 3
  cd "/scratch/user/uqmhal11/lrge-issue29/$acc"
  "$LRGE" -t 8 -vv -P "$plat" -s 4556 -T 10000 -Q 5000 "$@" \
      -o "${acc}.${tag}.size" "${acc}.fastq.gz" 2> "${acc}.${tag}.log"
  grep -F 'Estimate for' "${acc}.${tag}.log" | awk '{print $NF}' > "${acc}.${tag}.perread.txt"
  echo "[RESULT] $acc $tag estimate=$(cat "${acc}.${tag}.size") n_inf=$(grep -c inf "${acc}.${tag}.perread.txt")"
}
run SRR16767125 clean-default pb
run SRR16767125 clean-F pb -F
run SRR26465560 clean-default ont
run SRR8618952  clean-default pb
run SRR8618952  clean-F pb -F
run SRR26465526 clean-default ont
run SRR12247681 clean-default ont
echo "[DONE]"
