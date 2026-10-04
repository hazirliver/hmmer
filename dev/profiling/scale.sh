#!/usr/bin/env bash
# scale.sh <label> <reps> <hmmsearch args...>
#   wall/CPU/RSS of hmmsearch (build-dev, no perf) at --cpu 0,1,2,4,8,12,16; best of <reps>.
#   cpu_util = (user+sys)/wall; a plateau near 1 + (worker share) means the reader is the limit.
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
H=${HMMSEARCH:-$ROOT/build-dev/src/hmmsearch}
label=$1; reps=$2; shift 2
for c in ${CPUS:-0 1 2 4 8 12 16}; do
  best=999999; line=""
  for r in $(seq "$reps"); do
    out=$( { /usr/bin/time -f "%e %U %S %M" "$H" --cpu "$c" -o /dev/null "$@" > /dev/null; } 2>&1 | tail -1)
    w=${out%% *}
    if (( $(echo "$w < $best" | bc) )); then best=$w; line=$out; fi
  done
  read -r w u s m <<< "$line"
  printf "%-14s cpu=%-2s wall=%7.2f user=%7.2f sys=%5.2f cpu_util=%5.2fx rss=%6dMB\n" \
    "$label" "$c" "$w" "$u" "$s" "$(echo "($u+$s)/$w" | bc -l)" $((m/1024))
done
