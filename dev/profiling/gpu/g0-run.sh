#!/bin/bash
# G0 ceilings: all workloads of the task's list. Run detached; appends its exit status.
cd "$(dirname "$0")/../../.." || exit 1
REF=${REF:-5935bd9d-perf}; PAT=${PAT:-/home/asokolov/Projects/hmmer-g0/build/src/hmmsearch}; OUT=${OUT:-bench-results/g0}
mkdir -p $OUT
G="python3 dev/profiling/gpu/g0.py --ref $REF --patched $PAT --reps ${REPS:-3}"
$G --workloads '^(X5_human_pfam1k_cN|X1_p100_50k_cN|X2_m200_sprot_cN|X4_pfamA_h3m_cN)$' --out $OUT/main.json
$G --workloads '^(W1_hitters_10k|W5_med_1k_max)$' --cpu-override 16 --out $OUT/w.json
$G --file annot.tsv --workloads '^(A_human_L1|A_human_L4|A_bacteria_L2)$' --out $OUT/annot.json
echo "exit $?"
