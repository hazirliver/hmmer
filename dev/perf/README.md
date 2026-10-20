# A/B performance tooling for optimization steps

These tools answer two questions for every optimization step in [`../profiling/PLAN.md`](../profiling/PLAN.md):

1. Is the new hmmsearch **byte-identical** to the old one on real data?
2. Is it **faster beyond noise** on the workloads the step targets?

The golden regression (`dev/check.sh`) and the bench hashes (`dev/bench.sh`) check outputs against stored goldens. `ab.sh` instead compares two binaries directly, run by run, on large inputs.

```
dev/perf/snap.sh --rev HEAD --perf                 # A: the commit before the step -> bench-results/bins/<sha>-perf/
# ... change code, then dev/check.sh ...
B=$(dev/perf/snap.sh --perf)                       # B: the working tree, same flags
dev/perf/ab.sh <sha>-perf $B/hmmsearch --set cpu
dev/perf/ab.sh f8bebd8e --set cpu                  # B defaults to build-dev/src/hmmsearch (default flags)
dev/perf/ab.sh <sha> --set scale,query --reps 5
dev/perf/ab.sh <sha> build-xchk/src/hmmsearch --set all --check-only   # outputs only, e.g. a crosscheck build
dev/perf/ab.sh x --list --set all                  # show workloads
```

## `snap.sh`

```
dev/perf/snap.sh --rev HEAD --perf     # a commit, built in a temporary worktree, layout-stable flags
dev/perf/snap.sh --perf                # the working tree, same flags (build dir build-perf/)
dev/perf/snap.sh                       # the working tree as built in build-dev (default -O3)
dev/perf/snap.sh --cflags "-O3 ..."    # any flags; --name NAME overrides the snapshot name
```

- `--perf` means `CFLAGS="-O3 -falign-functions=64"`. **Use it on both sides of a step's A/B.**
  - In a default build, editing one file shifts every later function in the link by a few bytes. That changes the cache-line alignment of the SSE kernels and moves cycles by 1-3% on workloads that never run the edited code.
  - Measured on step 1.1: W5 (`--max`, bias filter not called) was 1.1-1.6% "slower" with an identical instruction count. With both sides aligned it was within noise.
- The working-tree builds reconfigure when `configure.ac`, a `Makefile.in` or a `*config.h.in` is newer than the build's Makefile (HMMER's Makefiles don't regenerate themselves), and a snapshot's `-dirty-` hash covers new untracked files too.
- `--rev` builds a commit in a temporary `git worktree`, with the Easel commit that the commit's own `dev/easel.ref` pins (a detached worktree of the `easel/` checkout, fetched from the ref's url if missing). It needs no stashing, and neither working tree is touched.
- Without `--rev`, it builds the working tree so the binaries match the source.
- Copies `hmmsearch` and the other search programs to `bench-results/bins/<name>/`, together with an `info.txt` (git sha, branch, Easel sha, flags, binary hash).
- `<name>` is the short sha, plus `-dirty-<hash of the diff>` when tracked files are modified, plus `-perf` (or `-cf<hash>` for other flags).
- HMMER links libhmmer and libeasel statically (only libc and libm are shared), so a copied binary is a complete snapshot.

## `ab.sh A [B]`

`A` and `B` can each be a binary, a directory that holds one, or a snapshot name or unique prefix under `bench-results/bins/`.

**Workloads** live in [`workloads.tsv`](workloads.tsv), grouped into sets:

| Set | Workloads | Verdict metric | Targets |
|---|---|---|---|
| `cpu` | W1-W5 at `--cpu 0` | `cycles:u` | Single-thread CPU work: filters and domain stage |
| `reader` | R1 `kern_M50 × sprot` at `--cpu 0` (~55% of its CPU is FASTA parsing) | `cycles_uk` (user + kernel cycles) | Sequence file reading (Phase 5): includes the system calls and page-cache copies, and unlike task-clock doesn't move with the CPU clock rate |
| `scale` | X1 `pfam_100 × sprot_50k` @ 4/8/N (plus @ N without `-Z`), X2 `kern_M200 × sprot` @ N | wall | Reader bottleneck, thread scaling |
| `query` | X3 `pfam_1k × sprot_1k` @ 0/N, X4 Pfam-A (text and `.h3m`) × 5 seqs @ 0/N | wall | Per-query fixed costs |
| `milestone` | X5 `pfam_1k × human --cut_ga` @ N | wall | Realistic annotation run |

A second file, [`annot.tsv`](annot.tsv) (`--file annot.tsv`), holds the **Pfam-A annotation benchmark**: all of Pfam-A with `--cut_ga` at N threads against six protein sets from UniProt (bacteria, phages, fungi, viruses of eukaryotes, eukaryotes other than human, human), each split by length into four groups of 250 (`L1` shortest to `L4` longest). That makes 24 workloads `A_<set>_L<g>` with the pressed Pfam-A (set `annot`), and the same 24 from the text file (`A_<set>_L<g>_txt`, set `annot_txt`). Upstream runs each model on a single worker here, because 250 targets fit in one read block, so these runs show per-model overhead and thread use more than kernel speed. The proteins come from `dev/bench/uniprot.py fetch` (see [`../bench/README.md`](../bench/README.md)); results are in [`../profiling/RESULTS.md`](../profiling/RESULTS.md), as sequences per second.

Data comes from the bench dataset (`dev/bench.sh prepare`), and `U/` in a workload stands for the UniProt sets. Fixtures such as the pressed Pfam-A are built on demand with A's `hmmpress`, in `~/.cache/hmmer-bench/ab-fixtures/`.

**How a comparison runs:**

- **Interleaving.** For each workload, the binaries alternate: A B, B A, A B, and so on. The first pair (`--warmup`, default 1) is untimed and serves as the output reference. Then come `--reps` timed pairs (default 5). Drift on the machine affects both sides equally.
- **Measurements.** Each run is wrapped in GNU `time` (user, sys, peak RSS). Workloads judged by a counter (sets `cpu` and `reader`) also run under `perf stat -e task-clock,instructions:u,cycles:u,cycles`. Workloads judged by wall time don't, unless you pass `--perf-all`. On this VM, counting a threaded run costs wall time: the per-thread counters are switched along with the threads. Pfam-A × 250 proteins at 16 threads measured +41% upstream (27.5 → 38.8 s) and +18% with Phase 7 (1.71 → 2.03 s). Until 2026-10-07 every timed run was counted, so wall-time verdicts on threaded runs carried that cost. Wall time is the parent's clock around a *blocking* wait. (Until 2026-10-05 it used `Popen.wait(timeout=...)`, which polls with sleeps that grow to 50 ms: every wall time was rounded up to the next poll, by up to 50 ms. That inflated the spreads of sub-second runs such as X2, about 5% of its 0.9 s, and is part of the wall noise in the table below. `dev/bench` had the same flaw and the same fix.)
- **argv[0].** The binary always runs as `hmmsearch`, through a symlink on `PATH`. Outputs therefore don't contain its path; the tabular tail echoes argv.
- **Output identity.**
  - Every file a run leaves, plus stdout and stderr, is hashed after the dev/bench normalization: dates, CPU time, Mc/s and thread count are masked.
  - A must equal B. Every rep must also equal its binary's first rep; a mismatch is reported as NONDETERMINISTIC.
  - Files that differ **only in blanks** (runs of spaces collapsed, trailing blanks ignored, text files only) are classified `ws`. They're accepted, printed as a note and shown as `ws` in the summary. Anything else is `DIFF`, which is a failure. This tolerates the upstream threaded column-width flake (`PLAN.md` §2) while still catching any change of content.
  - When files differ, the first lines of a unified diff are printed and the run dirs are kept under `bench-results/ab/work/`.
- **Verdict.**
  - ratio = median(B) / median(A).
  - A ratio counts as `faster` or `SLOWER` only when it differs from 1 by more than both the metric's floor (cycles and cycles_uk 0.5%, instructions 0.1%, wall 2%) and the spread (max−min)/median of either binary's reps.
  - Otherwise it's `~`.
- **Output and exit status.**
  - Results go to `bench-results/ab/<time>_<A>_vs_<B>[_label].json`, with every sample included.
  - Exit status: 0 = outputs identical; 1 = outputs differ, are nondeterministic, or a run failed; 3 = a slowdown, only with `--fail-on-slowdown`.

Options: `--file FILE` (workloads file, default `workloads.tsv`), `--workloads REGEX`, `--reps`, `--warmup`, `--check-only` (one untimed run each, outputs only), `--no-perf`, `--perf-all`, `--keep`, `--label`.

**B-only options.** `--b-args="--dbcache"` (use the `=` form) adds options to B's command line only, for example to compare a new opt-in feature against the old default. The option's own traces (its header echo line and its token in a table's `# Option settings:` line) are masked before comparing, so everything else must still be identical.

**Per-side environment.** `--a-env VAR=VALUE` and `--b-env VAR=VALUE` (repeatable) set environment variables for one side's runs. With the same binary on both sides, `--a-env HMMER_SIMD=sse` compares the SSE filter kernels with the default (widest) ones, with no code-layout difference between A and B. `HMMER_SIMD` and `HMMER_NCPU` are otherwise removed from the runs' environment.

## `series.py` and `plot_series.py`: many snapshots at once

`ab.sh` compares two binaries; `series.py` runs any number of them (for example upstream and the end of every phase) on the same workloads in one session, round-robin with a rotated order each round, after an untimed round that checks every binary's outputs against the first one's. It takes `--file` and `--perf-all` like `ab.sh`. `plot_series.py` turns its JSON into SVG charts (light and dark) and a Markdown table: `--charts` picks the speedup chart, the per-phase chart or both, and `--speedup LABELS` puts several binaries in the speedup chart (e.g. pyhmmer and this branch). [`../profiling/RESULTS.md`](../profiling/RESULTS.md) was made this way; its last section has the commands. `plot_annot.py` draws the annotation benchmark (one grouped bar chart per protein set, sequences per second, one bar per binary) for the same page. Both take `--best-of A,B` (show labels A and B as A, at the faster of the two per threaded workload, `*` in the tables), `--drop LABEL`, `--legend LABEL=TEXT` (a legend name, or the name of a `--best-of` alternative in tooltips) and `--note TEXT` (a footnote under the chart).

## `pyhmmer.sh` and `pyhmmsearch.py`: pyhmmer as a contender

[pyhmmer](https://github.com/althonos/pyhmmer) runs HMMER 3.4's pipeline from Python, with its own threading. To time it with the same tools:

```
dev/perf/pyhmmer.sh                 # venv ~/.cache/hmmer-bench/pyhmmer-0.12.3 + snapshot bench-results/bins/pyhmmer-0.12.3
dev/perf/pyhmmer.sh --ftz           # snapshot pyhmmer-0.12.3-ftz: the same with flush-to-zero on, as in hmmsearch
dev/perf/series.py series-p0 pyhmmer-0.12.3 series-p7e1 --labels upstream,pyhmmer,fork \
    --file vs_pyhmmer.tsv --mask-tail --flush-subnormals --workloads '^(X|W|R)'
dev/perf/ab.sh pyhmmer-0.12.3 pyhmmer-0.12.3 --b-env PYHMMSEARCH_CPUS=8 --file vs_pyhmmer.tsv --set vs   # thread counts
```

- **The snapshot's `hmmsearch`** is a wrapper around `pyhmmsearch.py`, run by the venv's `python -I`. The driver accepts the `hmmsearch` options that [`vs_pyhmmer.tsv`](vs_pyhmmer.tsv) uses, and rejects any other.
- **How it runs pyhmmer.** It reads the targets into memory once and streams the HMMs. With fewer HMMs than threads it passes them as a list, so that pyhmmer splits the targets among its threads. It writes `--tblout` / `--domtblout` with pyhmmer's `TopHits.write()`.
- **No main output.** pyhmmer can't write `hmmsearch`'s report, so every workload in `vs_pyhmmer.tsv` sends it to `-o /dev/null` and writes the two tables instead.
- **Time breakdown.** The driver writes where its time went to `.pyhmmsearch.json`. `ab.py` stores that in the sample, under `pyhmmer`.
- **Two output differences**, each masked by its own option, so that everything else is still compared exactly:
  - `--mask-tail` drops the `# Program: ... # [ok]` block that `hmmsearch` appends to its tables. pyhmmer doesn't write it.
  - `--flush-subnormals` compares numbers below 2.2e-308 as 0. `hmmsearch` sets flush-to-zero and denormals-are-zero (`impl_Init()`). pyhmmer doesn't, so it prints e.g. an E-value of `3.2e-315` where `hmmsearch` prints `0`. `series.py` notes how many such numbers each binary had. The `-ftz` snapshot has none, and its tables equal `hmmsearch`'s exactly.
  - Not masked: an E-value whose P-value is subnormal. `hmmsearch` flushes the P-value to 0 and prints an E-value of `0`; pyhmmer prints P × Z, a normal number such as `3.7e-306`. 2 of the 24 pressed annotation runs have one such hit and report `DIFFER`. `series.py` keeps the run directory of a binary whose files differ (not the reference's), under `bench-results/series/work/`.
- Plan and results: [`../profiling/PYHMMER.md`](../profiling/PYHMMER.md).

## `errpaths.sh A B`: error paths and query batching

Goldens cover few error exits, and a change in *when* output is printed (Phase 7's query batches) changes what is printed before a fatal error. `errpaths.sh` runs both binaries on multi-query searches that end in an error or exercise batching, at several thread counts, each in its own directory with the same file names, and compares every output (dates and timings masked):

```
dev/perf/errpaths.sh 8f490d2f-perf build-dev/src/hmmsearch                 # default: --cpu 0 4 16, sprot_10k
dev/perf/errpaths.sh A B --cpus "0 2" --db ~/.cache/hmmer-bench/derived-*/sprot_1k.fa
```

Cases: `--cut_ga` with GA lines missing on the first / 23rd model (with and without `--dbcache`), an HMM file truncated mid-model, mixed alphabets, a target file with a parse error in its 500th record (with and without `--dbcache`), and plain runs with and without `-Z`, `--cut_ga`, `--dbcache`, `-A`. Exit status 1 on any difference.

## Measured noise (A/A, same binary on both sides)

Machine: 16-vCPU Xeon 8468 VM, idle. Measured on 2026-10-04 with 5 reps per side.

| Workload | Metric | A/A ratio | Spread |
|---|---|---|---|
| W1_hitters_10k | cycles | 1.0015 | 0.75% |
| W2_mix_10k | cycles | 0.9994 | 0.38% |
| W3_m200_sprot | cycles | 0.9981 | 0.60% |
| W4_med_longseqs | cycles | 0.9981 | 0.47% |
| W5_med_1k_max | cycles | 0.9997 | 0.60% |

| X1_p100_50k_c4 | wall | 1.0002 | 5.76% |
| X1_p100_50k_c8 | wall | 0.9884 | 3.56% |
| X1_p100_50k_cN | wall | 1.0000 | 3.72% |
| X2_m200_sprot_cN | wall | 1.0466 | 4.72% |
| X3_p1k_s1k_c0 | wall | 1.0094 | 2.88% |
| X3_p1k_s1k_cN | wall | 0.9999 | 2.28% |
| X4_pfamA_txt_c0 | wall | 0.9944 | 1.87% |
| X4_pfamA_txt_cN | wall | 1.0042 | 1.79% |
| X4_pfamA_h3m_c0 | wall | 1.0001 | 3.03% |
| X4_pfamA_h3m_cN | wall | 1.0081 | 1.44% |

`instructions:u` was identical on both sides (ratio 1.0000) for every workload.

Added 2026-10-05: R1_m50_sprot (reader set) A/A `cycles_uk` 0.9975, spread 0.87% (5 reps). That day the VM's clock rate varied between runs: R1's task-clock spread was 8-14% and X2's wall spread 16-25% (blocking wait already in place), while cycles stayed within 1-2%. On such days only counter metrics resolve small steps.

What this means for verdicts:

- **cpu set (cycles).** A/A ratios stay within ±0.2%, so changes of about 0.75% or more are resolvable.
- **Threaded wall times.** These are noisier: A/A ratios reached ±4.7%, and spreads up to 5.8%. A scaling step therefore needs a clear effect, or more reps (`--reps 9`).
- **Wall time single-threaded vs threaded.** At `--cpu 16`, X4 takes 47 s versus 21 s at `--cpu 0`. That's the per-query overhead Phase 3 targets.

Threaded workloads pin `-Z`. Without it, hmmsearch output widths are nondeterministic upstream (see the comment in `workloads.tsv` and `PLAN.md` §2).

## Crosscheck builds (kernel rewrites)

When a step rewrites a numeric kernel, keep the old code under `#ifdef p7_CROSSCHECK`: run both versions and `abort()` on any bitwise difference. Build that flavor separately, then run it over the golden cases and the real-data workloads:

```
HMMER_BUILD_DIR=$PWD/build-xchk CFLAGS="-O3 -Dp7_CROSSCHECK" dev/check.sh --quick
dev/perf/ab.sh <sha> build-xchk/src/hmmsearch --set all --check-only
```

Remove the crosscheck code when the phase is merged.
