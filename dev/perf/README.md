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
- `--rev` builds a commit in a temporary `git worktree` that shares the pinned `easel/` checkout. It needs no stashing, and the working tree is never touched.
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
| `scale` | X1 `pfam_100 × sprot_50k` @ 4/8/N, X2 `kern_M200 × sprot` @ N | wall | Reader bottleneck, thread scaling |
| `query` | X3 `pfam_1k × sprot_1k` @ 0/N, X4 Pfam-A (text and `.h3m`) × 5 seqs @ 0/N | wall | Per-query fixed costs |
| `milestone` | X5 `pfam_1k × human --cut_ga` @ N | wall | Realistic annotation run |

Data comes from the bench dataset (`dev/bench.sh prepare`). Fixtures such as the pressed Pfam-A are built on demand with A's `hmmpress`, in `~/.cache/hmmer-bench/ab-fixtures/`.

**How a comparison runs:**

- **Interleaving.** For each workload, the binaries alternate: A B, B A, A B, and so on. The first pair (`--warmup`, default 1) is untimed and serves as the output reference. Then come `--reps` timed pairs (default 5). Drift on the machine affects both sides equally.
- **Measurements.** Each run is wrapped in `perf stat -e task-clock,instructions:u,cycles:u` and GNU `time` (wall, user, sys, peak RSS).
- **argv[0].** The binary always runs as `hmmsearch`, through a symlink on `PATH`. Outputs therefore don't contain its path; the tabular tail echoes argv.
- **Output identity.**
  - Every file a run leaves, plus stdout and stderr, is hashed after the dev/bench normalization: dates, CPU time, Mc/s and thread count are masked.
  - A must equal B. Every rep must also equal its binary's first rep; a mismatch is reported as NONDETERMINISTIC.
  - When files differ, the first lines of a unified diff are printed and the run dirs are kept under `bench-results/ab/work/`.
- **Verdict.**
  - ratio = median(B) / median(A).
  - A ratio counts as `faster` or `SLOWER` only when it differs from 1 by more than both the metric's floor (cycles 0.5%, instructions 0.1%, wall 2%) and the spread (max−min)/median of either binary's reps.
  - Otherwise it's `~`.
- **Output and exit status.**
  - Results go to `bench-results/ab/<time>_<A>_vs_<B>[_label].json`, with every sample included.
  - Exit status: 0 = outputs identical; 1 = outputs differ, are nondeterministic, or a run failed; 3 = a slowdown, only with `--fail-on-slowdown`.

Options: `--workloads REGEX`, `--reps`, `--warmup`, `--check-only` (one untimed run each, outputs only), `--no-perf`, `--keep`, `--label`.

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
