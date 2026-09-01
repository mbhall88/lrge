#!/usr/bin/env bash
# Build the pinned PR #45 binary once; every benchmark job invokes this absolute path.
set -euxo pipefail
WORKTREE=/scratch/temp/27749624/lrge-issue34
BENCH=/scratch/user/uqmhal11/lrge-issue34-benchmark
COMMIT=254e1462cfa910f507313fe8b718fcab320e275d

cd "$WORKTREE"
test "$(git rev-parse HEAD)" = "$COMMIT"
test -z "$(git status --porcelain)"

cargo +1.89.0 build --release --locked -p lrge
install -m 0755 target/release/lrge "$BENCH/bin/lrge-254e146"
"$BENCH/bin/lrge-254e146" --version
echo "[BUILD-OK] $COMMIT"
