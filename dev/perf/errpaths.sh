#!/usr/bin/env bash
# dev/perf/errpaths.sh A B [--cpus "0 4 16"] [--db FILE]
#
# Compare two hmmsearch binaries on multi-query runs that end in an error, or
# that exercise query batching, at several thread counts. A and B are paths
# (or snapshot names under bench-results/bins/). Each run executes in its own
# directory with the same file names, so the outputs can be compared byte for
# byte; dates and timings are masked. Prints one line per case, exits 1 on a
# difference.
#
# Cases, built from pfam_1k.hmm (first 40 models) in a fresh temporary dir:
#   --cut_ga with GA lines missing on model 1 / 23     (p7_pli_NewModel() fails)
#   HMM file truncated mid-model / mixed alphabets     (p7_hmmfile_Read() fails)
#   plain runs with and without -Z, --cut_ga, --dbcache, -A
set -uo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
D=$(ls -d ~/.cache/hmmer-bench/derived-* | head -1)
CPUS="0 4 16"; DB=$D/sprot_10k.fa
resolve() { if [[ -x $1 ]]; then realpath "$1"; elif [[ -x $1/hmmsearch ]]; then realpath "$1/hmmsearch"; else realpath "$ROOT/bench-results/bins/$1/hmmsearch"; fi; }
A=$(resolve "$1"); B=$(resolve "$2"); shift 2
while [[ $# -gt 0 ]]; do case $1 in --cpus) CPUS=$2; shift ;; --db) DB=$2; shift ;; *) echo "unknown option $1" >&2; exit 2 ;; esac; shift; done

W=$(mktemp -d); trap 'rm -r "$W"' EXIT
cd "$W"
awk 'BEGIN{n=0} /^HMMER3/{n++} n<=40{print}' "$D/pfam_1k.hmm" > p40.hmm
awk 'BEGIN{n=0} /^HMMER3/{n++} {if (n==1  && /^GA /) next; print}' p40.hmm > noga1.hmm
awk 'BEGIN{n=0} /^HMMER3/{n++} {if (n==23 && /^GA /) next; print}' p40.hmm > noga23.hmm
head -c $(( $(stat -c %s p40.hmm) * 3 / 5 )) p40.hmm > trunc.hmm
awk 'BEGIN{n=0} /^HMMER3/{n++} n<=20{print}' p40.hmm > p20.hmm
awk 'BEGIN{n=0} /^HMMER3/{n++} n<=1{print}' "$D/dfam_100.hmm" > dna1.hmm
cat p20.hmm dna1.hmm p20.hmm > mixed.hmm
mkdir dA dB; ln -s "$A" dA/hmmsearch; ln -s "$B" dB/hmmsearch
for f in *.hmm; do ln -s ../$f dA/$f; ln -s ../$f dB/$f; done

norm() { sed -E 's/^# (CPU time|Mc\/sec|Date|Current dir).*//' "$1"; }
bad=0
for t in "--cut_ga noga1.hmm" "--cut_ga noga23.hmm" "--dbcache --cut_ga noga23.hmm" "trunc.hmm" "mixed.hmm" \
         "--cut_ga p40.hmm" "p40.hmm" "-Z 10000 p40.hmm" "--dbcache p40.hmm" "--dbcache -Z 10000 --cut_ga p40.hmm"; do
  for c in $CPUS; do
    for x in dA dB; do
      (cd $x && ./hmmsearch --cpu "$c" --tblout r.tbl --domtblout r.dom -A r.sto -o r.out $t "$DB" > r.stdout 2> r.err; echo $? > r.rc)
    done
    diffs=""
    for f in tbl dom sto out stdout err rc; do cmp -s <(norm dA/r.$f) <(norm dB/r.$f) || diffs="$diffs $f"; done
    [[ -n $diffs ]] && bad=1
    verdict="identical"; [[ -n $diffs ]] && verdict="DIFF:$diffs"
    printf "%-38s cpu=%-2s rc=%s queries=%-3s %-10s %s\n" "[$t]" "$c" "$(cat dB/r.rc)" "$(grep -c '^Query:' dB/r.out)" \
      "$verdict" "$(head -c 60 dB/r.err | tr '\n' ' ')"
  done
done
exit $bad
