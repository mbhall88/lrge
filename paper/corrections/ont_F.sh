#!/usr/bin/env bash
# Does -F (post-#31) do anything for the ONT sub-0.5x failures, or are they purely mechanism 1?
set -euxo pipefail
run=$1; threads=${2:-8}
WORKDIR=/scratch/user/uqmhal11/lrge-issue29
LRGE=/scratch/user/uqmhal11/lrge/target/release/lrge
cd "$WORKDIR/$run"
variant () {
  local tag=$1; shift
  "$LRGE" -t "$threads" -vv -P ont -s 4556 -T 10000 -Q 5000 "$@" \
      -o "${run}.${tag}.size" "${run}.fastq.gz" 2> "${run}.${tag}.log"
  grep -F 'Estimate for' "${run}.${tag}.log" | awk '{print $NF}' > "${run}.${tag}.perread.txt"
  echo "[RESULT] $run $tag estimate=$(cat "${run}.${tag}.size") n_inf=$(grep -c inf "${run}.${tag}.perread.txt")"
}
variant fixed-F -F
variant fixed-F-strict -F --max-overhang-ratio 0.05
echo "[DONE] $run"
