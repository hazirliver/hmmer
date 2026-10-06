#!/usr/bin/env bash
# shard_queries.sh <hmmsearch> <P> <T> <outdir> [hmmsearch options...] <hmmfile> <seqdb>
#   Split <hmmfile> into P contiguous shards of about equal total model length
#   (sum of M), run P hmmsearch processes concurrently with --cpu T, and print
#   the wall time of the whole set. Each process writes <outdir>/s<i>.{out,tbl,dom}.
#   Every process sees all targets, so each query's results don't depend on the
#   split; concatenating the shards' bodies in order gives the single-run output.
set -euo pipefail
H=$1 P=$2 T=$3 OUT=$4; shift 4
args=("$@"); hmm=${args[-2]}; db=${args[-1]}; opts=("${args[@]:0:${#args[@]}-2}")
mkdir -p "$OUT"
awk -v P="$P" -v out="$OUT" '
  /^HMMER3/ { n++ } { rec[n] = rec[n] $0 "\n" }
  /^LENG/   { M[n] = $2; tot += $2 }
  END { s = 0; acc = 0
        for (i = 1; i <= n; i++) {
          if (acc >= (s + 1) * tot / P && s < P - 1) s++
          printf "%s", rec[i] > (out "/q" s ".hmm"); acc += M[i] } }' "$hmm"
t0=$(date +%s.%N)
for ((i = 0; i < P; i++)); do
  [[ -f $OUT/q$i.hmm ]] || continue
  "$H" --cpu "$T" "${opts[@]}" -o "$OUT/s$i.out" --tblout "$OUT/s$i.tbl" --domtblout "$OUT/s$i.dom" "$OUT/q$i.hmm" "$db" > /dev/null &
done
wait
t1=$(date +%s.%N)
printf "P=%-2d T=%-2d %s wall=%.2f s\n" "$P" "$T" "${opts[*]}" "$(echo "$t1 - $t0" | bc)"
