#!/bin/bash
cd "$(dirname "$0")/../../.." || exit 1
mkdir -p bench-results/g0
python3 dev/profiling/gpu/g0.py --ref ${REF:-5935bd9d-perf} --patched ${PAT:-/home/asokolov/Projects/hmmer-g0/build/src/hmmsearch} --reps ${REPS:-2} --file dev/profiling/gpu/g0-big.tsv --workloads "${WL:-.}" --out bench-results/g0/big-${WL_TAG:-all}.json
echo "exit $?"
