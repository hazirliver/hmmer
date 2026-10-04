# HMMER dev test suite

One command builds HMMER and runs every check we rely on to know nothing broke:

```
dev/check.sh
```

It exits 0 only if everything passes. It takes about 45 s on 16 cores. `dev/check.sh --quick` skips `make check` and takes about 5 s.

To also measure speed, run `dev/check.sh --bench [S|M|L]`. This adds a speed benchmark on large real datasets (Swiss-Prot, Pfam, Dfam, hg38 chr21, Rfam and E. coli). Each benchmark case also has its outputs hash-checked against goldens. Tier S takes about 5 min, M about 30 min and L hours. See [`bench/README.md`](bench/README.md).

For optimization work, `dev/perf/snap.sh` freezes the current binaries and `dev/perf/ab.sh` compares two builds run by run: outputs must be byte-identical on large real inputs, and speed is compared with interleaved runs and CPU counters. See [`perf/README.md`](perf/README.md) and the step plan in [`profiling/PLAN.md`](profiling/PLAN.md).

## What it does

| Step | What it runs | Fails when |
|---|---|---|
| **setup** | Clones Easel into `./easel` at the commit pinned in `dev/easel.ref` if it's missing, and runs `autoconf` if needed | Clone fails or a required tool is missing (warns if `easel/` is at a different commit) |
| **build** | Out-of-source `configure` + `make` + `make dev` in `build-dev/` (examples, benchmarks, utests) | Compile or link error |
| **make check** | Easel and HMMER unit tests plus the sqc integration suites (`testsuite/testsuite.sqc`, `easel/testsuite/testsuite.sqc`). That's 394 exercises at level 2. | Any exercise fails |
| **regress** | The golden-output regression (below) | Any output file isn't byte-identical to `dev/regress/golden/` |
| **bench** (optional) | The `dev/bench/` tier: speed per case, with output hashes checked against `dev/bench/golden.json` | A case fails, its output hash differs, or reps produce different outputs. A slowdown fails the run only with `--fail-on-slowdown` |

Logs go to `build-dev/{configure,make,make-dev,check}.log`.

## Golden-output regression

`dev/regress/` contains these files:

- `data/`: the fixed inputs. They're tutorial and testsuite alignments, plus `protdb.fa`, a synthetic protein db of 2,735 sequences. It contains real globins, planted homologs emitted from 5 profiles, shuffled decoys and i.i.d. random sequences. `make-data.sh` records how `protdb.fa` was made.
- `fixtures.txt`: commands that run serially first. They build the HMMs and query files that the cases use. These outputs are golden-checked too.
- `cases.txt`: 86 commands that run in parallel, each in its own directory. They cover every user-facing program: hmmbuild (weighting, effective-number, prior and `--hand` variants), hmmsearch (all table formats, `-A`, thresholds, `--max`, stdin, multi-query at several thread counts, `--restrictdb_*`, pressed query files, gzip and stdin targets), hmmpress, hmmscan, phmmer, jackhmmer (checkpoints), nhmmer (strands, FM index via makehmmerdb, alignment query), nhmmscan, hmmalign (all output formats, `--mapali`, `--trim`), hmmemit, hmmsim, hmmstat, hmmconvert (including the binary round-trip), hmmfetch, hmmlogo and alimask. The cases also cover error paths and threading. Several searches run with `--cpu 0`, the default and `--cpu 8`.
- `golden/`: the expected results. For each case it holds `stdout`, `stderr`, `exitcode` and every file the case wrote.

**Normalization:** only volatile lines are masked, in both text and binary files:

- `DATE` in HMM files and the embedded ctime in `.h3m` and `-b` files
- `# Date:`
- `# Current dir:`
- `# CPU time:`
- `# Mc/sec:`
- `run time:`

Everything else must match byte for byte, including scores, E-values, alignments, table formatting and binary pressed files.

The goldens are deterministic across runs, across thread counts and between `-O0` and `-O3` builds on x86-64 SSE with gcc 13. A different architecture (NEON or VMX) or compiler might legitimately differ in the last digit of some floats. If that happens, regenerate the goldens on that platform rather than editing them by hand.

## Workflow

```
git checkout dev && git checkout -b feature/foo
# ... hack ...
dev/check.sh                      # must be ALL PASSED before merging back to dev
```

If you change output **on purpose** (new feature, bug fix, format change):

```
dev/check.sh --update-golden            # or: --update-golden --cases 'hmmsearch'
git diff --stat dev/regress/golden      # review EVERY changed file
git add dev/regress/golden              # commit together with the code change
```

Give a golden change its own explanation in the commit message. An unexplained golden diff is a regression.

To add coverage, append a line to `cases.txt` with the format `<name> <shell command>`. Then run `dev/check.sh --quick --update-golden --cases '^name$'` and commit the new `golden/cases/<name>/`.

## Options

```
--quick            skip make check
--no-build         use the existing build as-is
--only-regress     golden regression only (no build, no make check)
--no-regress       skip golden regression
--update-golden    rewrite goldens from the current build
--cases REGEX      restrict the regression to matching cases
--full             make check at sqc level 3 (valgrind tests; needs valgrind)
--clean            rebuild from scratch
-j N               parallelism (default: nproc)
```

You can also set these environment variables:

- `HMMER_BUILD_DIR`: the build directory. The default is `build-dev`.
- `CONFIGURE_FLAGS`: for example `--enable-mpi`.
- `CFLAGS`

For example, `HMMER_BUILD_DIR=$PWD/build-dbg CFLAGS="-O0 -g" dev/check.sh` runs the suite against a separate debug build.

When the regression fails, the actual outputs stay in `build-dev/regress-work/actual/`. To see every difference, run `diff -r dev/regress/golden build-dev/regress-work/actual`.

## Not covered

- `hmmpgmd` and `hmmc2`, the daemon and client. The sqc suite exercises them partly (`hmmpgmd_ga`).
- MPI builds.
- Non-SSE SIMD backends.
