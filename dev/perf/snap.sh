#!/usr/bin/env bash
# dev/perf/snap.sh -- freeze the current build's search binaries as an A/B "before" snapshot.
#
#   dev/perf/snap.sh [--build-dir DIR] [--name NAME]
#
# Rebuilds DIR (default build-dev) so the binaries match the source tree, then copies
# the search programs to bench-results/bins/<name>/ and prints that directory.
# HMMER binaries link libhmmer/libeasel statically (only libc/libm are shared), so a
# copied binary is a complete snapshot. The default name is the short HEAD sha, plus
# -dirty-<hash of the diff> when tracked files are modified.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
BUILD="${HMMER_BUILD_DIR:-$ROOT/build-dev}" NAME=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --build-dir) BUILD=$(realpath "$2"); shift ;;
    --name)      NAME=$2; shift ;;
    -h|--help)   sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done
[[ -f "$BUILD/Makefile" ]] || { echo "no configured build in $BUILD (run dev/check.sh first)" >&2; exit 1; }

SRC=$(cd "$BUILD" && sed -n 's/^srcdir *= *//p' Makefile)
[[ "$SRC" == /* ]] || SRC="$BUILD/$SRC"
SRC=$(cd "$SRC" && pwd)
make -C "$BUILD" -j"$(nproc)" >/dev/null || { echo "build failed in $BUILD" >&2; exit 1; }

if [[ -z "$NAME" ]]; then
  NAME=$(git -C "$SRC" rev-parse --short=8 HEAD)
  dirty=$( { git -C "$SRC" diff HEAD; git -C "$SRC/easel" diff HEAD 2>/dev/null; } | sha256sum | cut -c1-8)
  if [[ -n "$(git -C "$SRC" status --porcelain --untracked-files=no)" || -n "$(git -C "$SRC/easel" status --porcelain --untracked-files=no 2>/dev/null)" ]]; then
    NAME="$NAME-dirty-$dirty"
  fi
fi

OUT="$ROOT/bench-results/bins/$NAME"
mkdir -p "$OUT"
for p in hmmsearch hmmscan phmmer jackhmmer nhmmer nhmmscan hmmpress hmmbuild; do
  [[ -x "$BUILD/src/$p" ]] && install -m 0555 "$BUILD/src/$p" "$OUT/$p.tmp" && mv -f "$OUT/$p.tmp" "$OUT/$p"
done
{
  echo "name        $NAME"
  echo "date        $(date -u +%FT%TZ)"
  echo "src         $SRC"
  echo "git_sha     $(git -C "$SRC" rev-parse HEAD)"
  echo "git_branch  $(git -C "$SRC" rev-parse --abbrev-ref HEAD)"
  echo "easel_sha   $(git -C "$SRC/easel" rev-parse HEAD 2>/dev/null || echo '?')"
  echo "build_dir   $BUILD"
  echo "flags       $(grep -E '^(CC|CFLAGS|SIMD_CFLAGS) *=' "$BUILD/src/Makefile" | tr -s ' ' | tr '\n' ';')"
  echo "hmmsearch   $(sha256sum "$OUT/hmmsearch" | cut -c1-16)"
} > "$OUT/info.txt"
echo "$OUT"
