#!/usr/bin/env bash
# profile_stage.sh <name> <cpu> <hmmsearch args...>
#   perf-records one hmmsearch run (build-prof), prints wall/CPU/RSS, the aggregated
#   pipeline pass rates and the per-stage CPU attribution (main=reader/output thread,
#   worker=pipeline threads). Artifacts land in $OUT (default bench-results/profiling).
# Note: perf adds a ~1.2 s wall-time floor; take wall times from scale.sh instead.
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
H=${HMMSEARCH:-$ROOT/build-prof/src/hmmsearch}
OUT=${OUT:-$ROOT/bench-results/profiling}; mkdir -p "$OUT"
name=$1; cpu=$2; shift 2
cd "$OUT"
/usr/bin/time -f "TIME wall=%e user=%U sys=%S maxrss_kb=%M" -o "$name.time" \
  perf record -q -F 997 -g -o "$name.data" "$H" --cpu "$cpu" -o "$name.out" "$@" 2> "$name.err"
{ echo "### $name  (--cpu $cpu $*)"; cat "$name.time"; python3 "$HERE/pstats.py" "$name.out"
  perf script -i "$name.data" -F comm,pid,tid,ip,sym 2>/dev/null | python3 "$HERE/stages.py" --threads | head -32; } > "$name.stages"
cat "$name.stages"
