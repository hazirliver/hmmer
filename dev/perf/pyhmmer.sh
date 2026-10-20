#!/usr/bin/env bash
# Install pyhmmer in its own venv and make a snapshot that series.py and ab.sh run like an hmmsearch binary.
#
#   dev/perf/pyhmmer.sh [--ftz] [--cpus N] [VERSION]       default 0.12.3
#
# venv:     ~/.cache/hmmer-bench/pyhmmer-VERSION (the PyPI wheel, pinned; pip resolves psutil)
# snapshot: bench-results/bins/pyhmmer-VERSION/hmmsearch, a wrapper that runs dev/perf/pyhmmsearch.py
#           with the venv's python -I (no user site, no current dir on sys.path), plus info.txt.
# --ftz:    snapshot pyhmmer-VERSION-ftz instead, with flush-to-zero and denormals-are-zero on, as in
#           hmmsearch (PYHMMSEARCH_FTZ=1, see pyhmmsearch.py); pyhmmer as installed leaves them off.
# --cpus N: snapshot pyhmmer-VERSION-cN, with N threads wherever a workload asks for --cpu above 0
#           (PYHMMSEARCH_CPUS=N), e.g. 8 = this VM's physical cores, pyhmmer's own default.
# Prints the snapshot dir. Compare with: dev/perf/series.py series-p0 pyhmmer-VERSION <fork> --mask-tail \
#   --file vs_pyhmmer.tsv (see dev/profiling/PYHMMER.md).
set -euo pipefail
ftz= cpus=
while [ $# -gt 0 ]; do
  case $1 in
    --ftz)  ftz=-ftz; shift ;;
    --cpus) cpus=${2:?--cpus N}; shift 2 ;;
    *)      break ;;
  esac
done
ver=${1:-0.12.3}
root=$(cd "$(dirname "$0")/../.." && pwd)
data=${HMMER_BENCH_DATA:-$HOME/.cache/hmmer-bench}
results=${HMMER_BENCH_RESULTS:-$root/bench-results}
venv=$data/pyhmmer-$ver
out=$results/bins/pyhmmer-$ver$ftz${cpus:+-c$cpus}

[ -x "$venv/bin/python" ] || python3 -m venv "$venv" >&2
"$venv/bin/pip" install -q --disable-pip-version-check "pyhmmer==$ver" >&2
got=$("$venv/bin/python" -I -c 'import pyhmmer; print(pyhmmer.__version__)')
[ "$got" = "$ver" ] || { echo "pyhmmer.sh: venv has pyhmmer $got, not $ver" >&2; exit 1; }

mkdir -p "$out"
cat > "$out/hmmsearch" <<EOF
#!/bin/sh
${ftz:+export PYHMMSEARCH_FTZ=1}
${cpus:+export PYHMMSEARCH_CPUS=$cpus}
exec "$venv/bin/python" -I "$root/dev/perf/pyhmmsearch.py" "\$@"
EOF
sed -i '/^$/d' "$out/hmmsearch"
chmod +x "$out/hmmsearch"
libs=$("$venv/bin/python" -I -c 'import pyhmmer, os; print(os.path.join(os.path.dirname(pyhmmer.__file__), "..", "pyhmmer.libs"))')
{
  printf 'name        %s\n' "$(basename "$out")"
  printf 'date        %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'src         PyPI wheel, venv %s\n' "$venv"
  printf 'packages    %s\n' "$("$venv/bin/pip" freeze --disable-pip-version-check | tr '\n' ' ')"
  printf 'python      %s\n' "$("$venv/bin/python" -I -c 'import sys; print(sys.version.split()[0])')"
  printf 'compiler    %s\n' "$(strings "$libs/liblibhmmer.so" | grep -m1 '^GCC:' || echo unknown)"
  printf 'simd        %s\n' "$(objdump -d "$libs/liblibhmmer.so" | grep -q ymm && echo 'AVX (ymm) present' || echo 'SSE only (no ymm/zmm in libhmmer)')"
  printf 'driver      dev/perf/pyhmmsearch.py %s%s\n' "$(sha256sum "$root/dev/perf/pyhmmsearch.py" | cut -c1-16)" "${ftz:+, PYHMMSEARCH_FTZ=1}${cpus:+, PYHMMSEARCH_CPUS=$cpus}"
  printf 'hmmsearch   %s\n' "$(cat "$out/hmmsearch" "$root/dev/perf/pyhmmsearch.py" "$libs/liblibhmmer.so" | sha256sum | cut -c1-16)"
} > "$out/info.txt"
echo "$out"
