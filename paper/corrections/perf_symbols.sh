#!/usr/bin/env bash
# Re-record with debug symbols. CARGO_PROFILE_RELEASE_DEBUG keeps release optimisation
# and only adds the debug info perf needs to name symbols. Separate target dir so the
# pinned benchmark binary is untouched.
set -euxo pipefail
WORKTREE=/scratch/temp/27749624/lrge-issue34
BENCH=/scratch/user/uqmhal11/lrge-issue34-benchmark
IN=/scratch/user/uqmhal11/lrge-issue29/SRR8618952/SRR8618952.fastq.gz
d="$BENCH/perf"; mkdir -p "$d"

cd "$WORKTREE"
CARGO_TARGET_DIR=/scratch/temp/27749624/lrge-issue34-symtarget \
CARGO_PROFILE_RELEASE_DEBUG=true \
CARGO_PROFILE_RELEASE_STRIP=false \
  cargo +1.89.0 build --release --locked -p lrge
BIN=/scratch/temp/27749624/lrge-issue34-symtarget/release/lrge

cd "$d"
t=$(mktemp -d -p "$d")
perf record -F 499 -o "$d/sym.data" -- \
    "$BIN" -P pb -s 4556 -T 10000 -Q 5000 -t 8 -D "$t" --normalize auto "$IN" > /dev/null 2> "$d/sym.err"
rm -rf "$t"
perf report -i "$d/sym.data" --stdio --no-children --sort symbol --percent-limit 0.5 2>/dev/null \
    | grep -v '^#' | grep -v '^$' | head -35 > "$d/sym_report.txt"
echo "[SYM-DONE]"
cat "$d/sym_report.txt"
