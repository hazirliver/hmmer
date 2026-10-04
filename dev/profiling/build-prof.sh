#!/usr/bin/env bash
# Build a profiling copy of HMMER in build-prof/: same -O3 as build-dev, plus
# symbols and frame pointers so perf call chains can attribute time to stages.
# Outputs are byte-identical to build-dev (checked when this was written).
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
mkdir -p "$ROOT/build-prof" && cd "$ROOT/build-prof"
[[ -f Makefile ]] || ../configure CFLAGS="-O3 -g -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer" > configure.log
make -j"$(nproc)" > make.log
echo "built $ROOT/build-prof/src/hmmsearch"
