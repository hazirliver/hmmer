#!/usr/bin/env bash
# dev/perf/snap.sh -- build and freeze hmmsearch & co. as an A/B snapshot.
#
#   dev/perf/snap.sh                          working tree, build-dev (default flags)
#   dev/perf/snap.sh --perf                   working tree, layout-stable perf flags
#   dev/perf/snap.sh --rev HEAD --perf        a commit (temporary git worktree), perf flags
#   dev/perf/snap.sh --cflags "-O3 ..."       any flags; --name NAME to override the name
#
# Copies the search programs to bench-results/bins/<name>/ and prints that dir.
# Binaries link libhmmer/libeasel statically (only libc/libm are shared), so a copy
# is a complete snapshot. Name: short sha, plus -dirty-<diff hash> for a modified
# working tree, plus -<flag tag> for non-default flags.
#
# --perf = CFLAGS "-O3 -falign-functions=64". Every function starts on a cache line,
# so editing one file doesn't shift the alignment of unrelated hot kernels (that
# alone moves cycles by 1-3%). Use it on both sides of a step's A/B.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
PERF_CFLAGS="-O3 -falign-functions=64"
BUILD="" NAME="" REV="" CFLAGS_SET=0 CF="-O3"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --build-dir) BUILD=$(realpath -m "$2"); shift ;;
    --name)      NAME=$2; shift ;;
    --rev)       REV=$2; shift ;;
    --cflags)    CF=$2; CFLAGS_SET=1; shift ;;
    --perf)      CF=$PERF_CFLAGS; CFLAGS_SET=1 ;;
    -h|--help)   sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

tag=""
if [[ $CFLAGS_SET == 1 && "$CF" != "-O3" ]]; then
  if [[ "$CF" == "$PERF_CFLAGS" ]]; then tag="-perf"; else tag="-cf$(printf '%s' "$CF" | sha256sum | cut -c1-6)"; fi
fi

configure_and_make() { # <srcdir> <builddir>
  local src=$1 b=$2 stamp
  mkdir -p "$b"; stamp="$b/.snap-cflags"
  if [[ ! -f "$b/Makefile" || ! -f "$stamp" || "$(cat "$stamp")" != "$CF" ]]; then
    [[ -f "$src/configure" && ! "$src/configure.ac" -nt "$src/configure" ]] || (cd "$src" && autoconf)
    (cd "$b" && "$src/configure" CFLAGS="$CF" > configure.log 2>&1) || { tail -20 "$b/configure.log" >&2; exit 1; }
    echo "$CF" > "$stamp"
  fi
  make -C "$b" -j"$(nproc)" > "$b/make.log" 2>&1 || { tail -30 "$b/make.log" >&2; echo "build failed in $b" >&2; exit 1; }
}

if [[ -n "$REV" ]]; then
  # Build a commit in a temporary worktree that shares the pinned Easel checkout.
  sha=$(git -C "$ROOT" rev-parse --verify "$REV^{commit}")
  WT=$(mktemp -d "${TMPDIR:-/tmp}/hmmer-snap.XXXXXX")
  trap 'git -C "$ROOT" worktree remove --force "$WT" 2>/dev/null || rm -rf "$WT"' EXIT
  git -C "$ROOT" worktree add -q --detach "$WT" "$sha"
  ln -s "$ROOT/easel" "$WT/easel"
  SRC=$WT; BUILD="$WT/build"
  configure_and_make "$SRC" "$BUILD"
  [[ -z "$NAME" ]] && NAME="$(git -C "$ROOT" rev-parse --short=8 "$sha")$tag"
  dirty=""
else
  SRC=$ROOT
  if [[ -z "$BUILD" ]]; then
    if [[ $CFLAGS_SET == 1 && "$CF" != "-O3" ]]; then BUILD="$ROOT/build${tag}"; else BUILD="${HMMER_BUILD_DIR:-$ROOT/build-dev}"; fi
  fi
  if [[ -f "$BUILD/Makefile" && $CFLAGS_SET == 0 ]]; then
    make -C "$BUILD" -j"$(nproc)" >/dev/null || { echo "build failed in $BUILD" >&2; exit 1; }
    CF=$(sed -n 's/^CFLAGS *= *//p' "$BUILD/src/Makefile")
  else
    configure_and_make "$SRC" "$BUILD"
  fi
  if [[ -z "$NAME" ]]; then
    NAME=$(git -C "$SRC" rev-parse --short=8 HEAD)
    if [[ -n "$(git -C "$SRC" status --porcelain --untracked-files=no)" || -n "$(git -C "$SRC/easel" status --porcelain --untracked-files=no 2>/dev/null)" ]]; then
      NAME="$NAME-dirty-$( { git -C "$SRC" diff HEAD; git -C "$SRC/easel" diff HEAD 2>/dev/null; } | sha256sum | cut -c1-8)"
    fi
    NAME="$NAME$tag"
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
  echo "src         ${REV:+$REV @ }$SRC"
  echo "git_sha     $(git -C "$SRC" rev-parse HEAD)"
  echo "git_branch  $(git -C "$ROOT" rev-parse --abbrev-ref HEAD)"
  echo "easel_sha   $(git -C "$ROOT/easel" rev-parse HEAD 2>/dev/null || echo '?')"
  echo "build_dir   $BUILD"
  echo "cflags      $CF"
  echo "hmmsearch   $(sha256sum "$OUT/hmmsearch" | cut -c1-16)"
} > "$OUT/info.txt"
echo "$OUT"
