#!/usr/bin/env bash
# dev/check.sh -- one-command dev test suite for HMMER.
#
#   dev/check.sh                 setup + build + `make check` + golden regression
#   dev/check.sh --quick         skip `make check` (build + golden regression only)
#   dev/check.sh --update-golden regenerate golden outputs (review the git diff!)
#
# See dev/README.md for details.
set -uo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
DEV="$ROOT/dev"
REG="$DEV/regress"
BUILD="${HMMER_BUILD_DIR:-$ROOT/build-dev}"
WORK="$BUILD/regress-work"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"

DO_SETUP=1 DO_BUILD=1 DO_CHECK=1 DO_REGRESS=1 UPDATE=0 SQCLEVEL=2 CLEAN=0 FILTER=""

usage() {
  cat <<EOF
Usage: dev/check.sh [options]
  --quick            skip 'make check' (unit tests + sqc integration tests)
  --no-build         don't (re)configure/build; use existing $BUILD
  --only-regress     only run the golden regression (implies --no-build)
  --no-regress       skip the golden regression
  --update-golden    rewrite dev/regress/golden from the current build
  --cases REGEX      only run regression cases whose name matches REGEX
  --full             'make check' at sqc level 3 (adds valgrind tests; needs valgrind)
  --clean            wipe the build dir and rebuild from scratch
  -j N               parallel jobs (default: $JOBS)
  -h, --help         this help
Env: HMMER_BUILD_DIR (default build-dev), CONFIGURE_FLAGS, CFLAGS
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --quick)         DO_CHECK=0 ;;
    --no-build)      DO_BUILD=0 ;;
    --only-regress)  DO_BUILD=0; DO_CHECK=0 ;;
    --no-regress)    DO_REGRESS=0 ;;
    --update-golden) UPDATE=1 ;;
    --cases)         FILTER="$2"; shift ;;
    --full)          SQCLEVEL=3 ;;
    --clean)         CLEAN=1 ;;
    -j)              JOBS="$2"; shift ;;
    -j*)             JOBS="${1#-j}" ;;
    -h|--help)       usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

if [[ -t 1 ]]; then R=$'\e[31m' G=$'\e[32m' Y=$'\e[33m' B=$'\e[1m' N=$'\e[0m'; else R= G= Y= B= N=; fi
step() { echo; echo "${B}==> $*${N}"; }
die()  { echo "${R}ERROR: $*${N}" >&2; exit 1; }

declare -a SUMMARY=()
FAILED=0
record() { # record <name> <ok|FAIL|skip> [detail]
  local color=$G; [[ $2 == FAIL ]] && { color=$R; FAILED=1; }; [[ $2 == skip ]] && color=$Y
  SUMMARY+=("$(printf '  %-22s %s%-4s%s %s' "$1" "$color" "$2" "$N" "${3:-}")")
}

# --------------------------------------------------------------------------
# 1. Setup: Easel checkout pinned to dev/easel.ref, autoconf.
# --------------------------------------------------------------------------
setup() {
  step "Setup"
  local url ref
  url=$(awk '$1=="url"{print $2}' "$DEV/easel.ref")
  ref=$(awk '$1=="commit"{print $2}' "$DEV/easel.ref")
  if [[ ! -d "$ROOT/easel/.git" ]]; then
    echo "cloning Easel ($ref)"
    git clone -q "$url" "$ROOT/easel" || die "easel clone failed"
    git -C "$ROOT/easel" checkout -q "$ref" || die "easel checkout $ref failed"
  fi
  local have; have=$(git -C "$ROOT/easel" rev-parse HEAD)
  if [[ "$have" != "$ref" ]]; then
    echo "${Y}WARNING: easel/ is at $have, dev/easel.ref pins $ref.${N}"
    echo "${Y}         Golden outputs were generated with the pinned commit.${N}"
  fi
  for t in gcc make autoconf python3 perl; do command -v $t >/dev/null || die "missing required tool: $t"; done
  if [[ ! -f "$ROOT/configure" || "$ROOT/configure.ac" -nt "$ROOT/configure" ]]; then
    (cd "$ROOT" && autoconf) || die "autoconf failed"
  fi
  record setup ok "easel $(git -C "$ROOT/easel" rev-parse --short HEAD)"
}

# --------------------------------------------------------------------------
# 2. Build (out-of-source, in $BUILD).
# --------------------------------------------------------------------------
build() {
  step "Build ($BUILD, -j$JOBS)"
  if [[ $CLEAN == 1 && -d "$BUILD" ]]; then chmod -R u+w "$BUILD"; rm -rf "$BUILD"; fi
  mkdir -p "$BUILD"
  local stamp="$BUILD/.configure-args" args="${CONFIGURE_FLAGS:-} CFLAGS=${CFLAGS:-}"
  if [[ ! -f "$BUILD/Makefile" || ! -f "$stamp" || "$(cat "$stamp")" != "$args" || "$ROOT/configure" -nt "$BUILD/Makefile" ]]; then
    echo "configuring: ${CONFIGURE_FLAGS:-<default>}"
    (cd "$BUILD" && "$ROOT/configure" ${CONFIGURE_FLAGS:-} ${CFLAGS:+CFLAGS="$CFLAGS"} > configure.log 2>&1) \
      || { tail -30 "$BUILD/configure.log"; record build FAIL "configure failed (see $BUILD/configure.log)"; return 1; }
    echo "$args" > "$stamp"
  fi
  if (cd "$BUILD" && make -j"$JOBS" > make.log 2>&1 && make -j"$JOBS" dev > make-dev.log 2>&1); then
    local nw; nw=$(cat "$BUILD/make.log" "$BUILD/make-dev.log" | grep -c 'warning:' || true)
    record build ok "$nw compiler warning(s) in files rebuilt this run (see $BUILD/make*.log)"
  else
    tail -40 "$BUILD/make.log" "$BUILD/make-dev.log" 2>/dev/null
    record build FAIL "see $BUILD/make.log"; return 1
  fi
}

# --------------------------------------------------------------------------
# 3. make check: Easel + HMMER unit tests and sqc integration tests.
# --------------------------------------------------------------------------
make_check() {
  step "make check (SQCLEVEL=$SQCLEVEL)"
  if [[ $SQCLEVEL -ge 3 ]] && ! command -v valgrind >/dev/null; then
    echo "${Y}valgrind not found: level-3 valgrind tests will fail${N}"
  fi
  (cd "$BUILD" && make -j"$JOBS" check SQCLEVEL=$SQCLEVEL > check.log 2>&1)
  local rc=$?
  local summ; summ=$(grep -E 'exercises at level|exercises at level .* FAILED|FAILED' "$BUILD/check.log" | grep -E '^All|FAILED' | tr '\n' ';')
  if [[ $rc -eq 0 ]] && ! grep -q 'FAILED' "$BUILD/check.log"; then
    grep -E '^All [0-9]+ exercises' "$BUILD/check.log"
    record "make check" ok "$(grep -E '^All [0-9]+ exercises' "$BUILD/check.log" | awk '{s+=$2} END{print s" sqc exercises passed"}')"
  else
    grep -E 'FAILED|failed|Error' "$BUILD/check.log" | head -40
    record "make check" FAIL "${summ:-rc=$rc} (see $BUILD/check.log)"
  fi
}

# --------------------------------------------------------------------------
# 4. Golden-output regression.
# --------------------------------------------------------------------------
# Normalize volatile lines in place (text files only). Keep in sync with README.
normalize() {
  local f
  for f in "$@"; do
    [[ -f "$f" && ! -L "$f" ]] || continue
    if grep -Iq . "$f" 2>/dev/null || [[ ! -s "$f" ]]; then
      sed -i -E \
        -e 's/^(DATE  ).{24}$/\1Thu Jan  1 00:00:00 1970/' \
        -e 's/^# Date:.*/# Date:            [normalized]/' \
        -e 's/^# Current dir:.*/# Current dir:     [normalized]/' \
        -e 's/^# CPU time:.*/# CPU time: [normalized]/' \
        -e 's/^# Mc\/sec:.*/# Mc\/sec: [normalized]/' \
        -e 's/^run time: .*/run time: [normalized]/' \
        "$f"
    else
      # binary (e.g. hmmpress .h3m/.h3p, hmmconvert -b): ctime is a 24-char C string
      perl -0777 -pi -e 's/(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun) (?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [ \d]\d \d\d:\d\d:\d\d \d{4}/Thu Jan  1 00:00:00 1970/g' "$f"
    fi
  done
}

# run_one <dir> <name> <command> : run command in dir, capture stdout/stderr/exit
run_one() {
  local dir=$1 name=$2 cmd=$3 rc
  ( cd "$dir" && bash -c "$cmd" > "$dir/stdout" 2> "$dir/stderr" ); rc=$?
  echo "$rc" > "$dir/exitcode"
}

regress() {
  local out="$WORK/actual" golden="$REG/golden"
  step "Golden regression ($([[ $UPDATE == 1 ]] && echo UPDATING golden || echo comparing to dev/regress/golden))"
  [[ -x "$BUILD/src/hmmsearch" ]] || { record regress FAIL "no build in $BUILD"; return 1; }
  [[ -d "$WORK" ]] && chmod -R u+w "$WORK"
  rm -rf "$WORK"; mkdir -p "$WORK/data" "$out/fx" "$out/cases"
  cp -r "$REG/data/." "$WORK/data/" && chmod -R a-w "$WORK/data"

  export PATH="$BUILD/src:$BUILD/easel/miniapps:$PATH"
  unset HMMER_NCPU
  export LC_ALL=C TZ=UTC

  # fixtures (serial)
  local name cmd
  ln -s "$WORK/data" "$out/fx/data"
  mkdir -p "$out/fx/_log"
  while read -r name cmd; do
    run_one "$out/fx" "$name" "$cmd"
    for s in stdout stderr exitcode; do mv "$out/fx/$s" "$out/fx/_log/$name.$s"; done
    if [[ $(cat "$out/fx/_log/$name.exitcode") != 0 ]]; then
      cat "$out/fx/_log/$name.stderr"; record regress FAIL "fixture '$name' failed"; return 1
    fi
  done < <(grep -Ev '^\s*(#|$)' "$REG/fixtures.txt")
  normalize $(find "$out/fx" -type f)
  rm "$out/fx/data"
  chmod -R a-w "$out/fx"

  # cases (parallel)
  local -a names=()
  while read -r name cmd; do
    [[ -n "$FILTER" && ! "$name" =~ $FILTER ]] && continue
    names+=("$name")
    mkdir -p "$out/cases/$name"
    ln -s "$WORK/data" "$out/cases/$name/data"
    ln -s "$out/fx" "$out/cases/$name/fx"
    printf '%s\n' "$cmd" > "$WORK/cmd.$name"
  done < <(grep -Ev '^\s*(#|$)' "$REG/cases.txt")
  export -f run_one
  local t0=$SECONDS
  printf '%s\n' "${names[@]}" | xargs -P "$JOBS" -I{} bash -c 'run_one "$1/{}" {} "$(cat "$2/cmd.{}")"' _ "$out/cases" "$WORK"
  for name in "${names[@]}"; do
    find "$out/cases/$name" -maxdepth 1 -type l -delete
    normalize $(find "$out/cases/$name" -type f)
  done
  echo "ran ${#names[@]} cases in $((SECONDS - t0))s"
  chmod -R u+w "$WORK"   # inputs were read-only only while cases ran

  if [[ $UPDATE == 1 ]]; then
    if [[ -n "$FILTER" ]]; then
      for name in "${names[@]}"; do rm -rf "$golden/cases/$name"; mkdir -p "$golden/cases"; cp -r "$out/cases/$name" "$golden/cases/"; done
    else
      rm -rf "$golden"; mkdir -p "$golden"; cp -r "$out/fx" "$out/cases" "$golden/"
    fi
    chmod -R u+w "$golden"
    record regress ok "golden updated: ${#names[@]} cases (review: git status dev/regress/golden)"
    return 0
  fi

  [[ -d "$golden" ]] || { record regress FAIL "no golden outputs; run with --update-golden first"; return 1; }
  local -a bad=()
  local scope="fx cases"
  [[ -n "$FILTER" ]] && scope="fx $(printf 'cases/%s ' "${names[@]}")"
  local s
  for s in $scope; do
    local d
    d=$(diff -r -q "$golden/$s" "$out/$s" 2>&1) || bad+=("$d")
  done
  if [[ ${#bad[@]} -eq 0 ]]; then
    record regress ok "${#names[@]} cases + fixtures byte-identical to golden"
  else
    echo "${R}Differences from golden:${N}"
    printf '%s\n' "${bad[@]}" | sed "s|$golden/||g; s|$out/||g" | head -50
    local first
    first=$(printf '%s\n' "${bad[@]}" | grep -m1 '^Files ' | awk '{print $2" "$4}')
    if [[ -n "$first" ]]; then
      echo; echo "${B}first diff (golden -> actual):${N}"
      diff -u $first | head -40
    fi
    local nfiles; nfiles=$(printf '%s\n' "${bad[@]}" | grep -c . )
    echo; echo "Actual outputs kept in $out  (compare: diff -r $golden $out)"
    record regress FAIL "$nfiles file(s) differ from golden"
  fi
}

# --------------------------------------------------------------------------
cd "$ROOT"
T0=$SECONDS
[[ $DO_SETUP == 1 ]] && setup
if [[ $DO_BUILD == 1 ]]; then build || { DO_CHECK=0; DO_REGRESS=0; }; fi
if [[ $DO_CHECK == 1 ]];   then make_check; else record "make check" skip; fi
if [[ $DO_REGRESS == 1 ]]; then regress;    else record regress skip; fi

echo; echo "${B}==> Summary ($((SECONDS - T0))s)${N}"
printf '%s\n' "${SUMMARY[@]}"
if [[ $FAILED == 0 ]]; then echo "${G}${B}ALL PASSED${N}"; else echo "${R}${B}FAILURES${N}"; fi
exit $FAILED
