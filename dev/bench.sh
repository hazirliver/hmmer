#!/usr/bin/env bash
# Speed benchmark + output-regression on large real datasets.
#   dev/bench.sh prepare           download pinned data (~2 GB) and derive datasets (once)
#   dev/bench.sh run [S|M|L]       benchmark the build in build-dev (default tier S)
#   dev/bench.sh compare A.json B.json
# See dev/bench/README.md.
exec python3 "$(dirname "${BASH_SOURCE[0]}")/bench/bench.py" "$@"
