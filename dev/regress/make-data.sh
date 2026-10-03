#!/usr/bin/env bash
# Provenance of dev/regress/data/protdb.fa (the synthetic protein target db).
# It is committed and must NOT be regenerated casually: goldens depend on it.
# Run from the repo root with a 3.4 build in build-dev:  dev/regress/make-data.sh
set -euo pipefail
B=${HMMER_BUILD_DIR:-build-dev}; E=$B/easel/miniapps; D=dev/regress/data
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
for m in globins4 fn3 Pkinase RRM_1 LuxC; do
  $B/src/hmmbuild "$tmp/$m.hmm" "$D/$m.sto" > /dev/null
  $B/src/hmmemit --seed 11 -N 40 "$tmp/$m.hmm" | sed "s/^>/>emit_${m}_/"
done > "$tmp/emit.fa"
$E/esl-reformat fasta "$D/7LESS_DROME" > "$tmp/7less.fa"
cat "$D/globins45.fa" "$tmp/7less.fa" "$D/HBB_HUMAN" "$tmp/emit.fa" > "$tmp/real.fa"
$E/esl-shuffle --seed 3 -N 4 "$tmp/real.fa"                | sed 's/^>/>shuf_/' > "$tmp/shuf.fa"
$E/esl-shuffle --seed 5 -G --amino -N 1500 -L 350          | sed 's/^>/>rnd_/'  > "$tmp/rnd.fa"
cat "$tmp/real.fa" "$tmp/shuf.fa" "$tmp/rnd.fa" > "$tmp/all.fa"
$E/esl-reformat fasta "$tmp/all.fa" > "${1:-$D/protdb.fa}"
