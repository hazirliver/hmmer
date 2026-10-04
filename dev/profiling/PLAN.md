# hmmsearch acceleration plan (small, verified steps)

Input: the profile in [`README.md`](README.md) (dev @ 8da3ab2, measured 2026-10-04).
Each step below is one small commit. A step is kept only if it is **byte-identical** to the previous step on every check and **faster beyond noise** on the metric it targets. A pure refactor that later steps need can be kept when it is perf-neutral.

---

## 1. What the profile says, in one paragraph

On realistic Pfam × proteome runs the **serial reader limits wall time**. It re-parses the whole target FASTA for every query at ~220 M res/s. With Pfam-like models (M̄ ≈ 163), hmmsearch can't use more than about 3 cores, and 16 threads are slower than 8. Second come the **per-query fixed costs**: a pthread create/join per query, per-thread clone/alloc, and text HMM parsing. Together these dominate many-model × small-DB runs. CPU work splits into filters S1-S5 (~57%), the domain stage (~15%, but 20-40% on hit-rich, long-model and `--max` runs) and the reader (~26%).

Two consequences for planning:

- **At 16 threads, CPU-only optimizations don't show up in wall time today.** Wall time is reader-bound until the reader is fixed. CPU steps must therefore be measured with CPU counters, not with 16-thread wall time.
- **The order of work matters.** First remove the reader and per-query serial costs (Phases 2-3). Only then do CPU savings (Phases 1, 4, 6) turn into wall time at scale.

---

## 2. How byte-identical testing works here (and what was checked)

| Gate | Command | What it proves | Cost |
|---|---|---|---|
| Golden regression | `dev/check.sh --quick` | 75 cases plus fixtures across all programs, byte-for-byte against `dev/regress/golden`. Only dates and timings are masked. Includes hmmsearch at `--cpu 0`, default and `8`, plus `--max`, `--nobias --nonull2`, all table formats and `-A` | **~3 s** (verified green on f8bebd8) |
| Unit + integration | `dev/check.sh` | Adds `make check` (394 sqc exercises) | ~45 s |
| Real-data output hashes | `dev/bench.sh run S --cases hmmsearch --reps 1 --baseline none` | hmmsearch on Swiss-Prot subsets and Pfam sets. Normalized outputs hash-match `dev/bench/golden.json`, and reps are deterministic | ~1-2 min |
| Milestone | `dev/check.sh --bench M` | Everything, plus 121 bench cases with hashes and speed vs baseline | ~35 min |

**Invariants that make reordering work safe** (checked in the code):

- Domain definition reseeds its RNG at every region (`p7_domaindef.c:612`, `do_reseeding` whenever `--seed` ≠ 0). Results therefore don't depend on which thread sees which sequence, or in what order. `--seed 0` is nondeterministic already and isn't golden-tested.
- Hits are sorted by sortkey, then name, then position (`p7_tophits.c:310`). Merge order doesn't matter.
- On the protein path, `esl_sqio_ReadBlock` does **not** reset block slots: `sqascii_Read` (`esl_sqio_ascii.c:1474`) expects them already `esl_sq_Reuse`'d. Only the DNA/long-target branch resets them itself (:1546). Today the worker's `esl_sq_Reuse(dbsq)` provides that reset. Any change that keeps blocks around (the cache) has to move that responsibility, not just drop it.

**Known upstream nondeterminism (found by the Phase 0 A/A run).** With `-Z` unset and E-value thresholds, each worker sets `pli->Z` to the number of sequences *it* has seen so far (`p7_pipeline.c:580`). Which unreported borderline hits stay in the hit list therefore depends on thread scheduling.

- Reported hits, scores and E-values are unaffected, because `p7_tophits_Threshold` recomputes them with the final Z.
- Column widths are affected. They're computed over *all* kept hits (`p7_tophits.c:1134`, `GetMaxNameLength` for `--tblout`), so the main output and the tables can be a character wider or narrower from run to run.
- Seen once in about 6 runs of `pfam_100 × sprot_50k --cpu 4`. It never showed in 30 runs each of the threaded golden cases on `protdb`.

Consequences:

- Threaded `ab.sh` workloads pass `-Z <true target count>`.
- `dev/bench` cases that use `--cpu N` without `-Z` can, rarely, report NONDETERMINISTIC regardless of the change being tested.
- A fix (each worker uses the sequence's global index as Z, which makes threaded output equal `--cpu 0` output) would change threaded outputs. It's listed in §6 as a decision, not done.

**Gaps to close before touching the code (Step 0.3):**

- No golden case covers multi-query × several worker threads × several 1000-sequence blocks per worker.
- No case covers `-Z` left unset in multi-query threaded mode, `--restrictdb_stkey`/`--restrictdb_n`, binary or pressed HMM input to hmmsearch, or the "isn't rewindable" error for multi-query on `.gz`/stdin.

**Measurement noise** (measured today, same code as the 10-03 baseline):

| Metric | Spread | Use it for |
|---|---|---|
| `instructions:u` (`perf stat -r 5`, `--cpu 0`) | 0.00% | Detecting any change in executed work (deterministic) |
| `cycles:u` | 0.23% | CPU-time effect of single-thread steps |
| `task-clock` | 0.6% | — |
| Wall, within one run (min → median) | 1-4% | Scaling steps, only interleaved A/B |
| Wall, day to day (identical code) | up to 4.5% (`human_pfam1k` 59.5 → 56.9 s) | Never compare against an old baseline file for a small step |

`hmmsearch` links only libc and libm (libhmmer and libeasel are static). A single copied binary is therefore a complete "before" snapshot.

---

## 3. Per-step protocol

```
git checkout -b perf/<id>-<slug> dev        # or one branch per phase, one commit per step
dev/perf/snap.sh --rev HEAD --perf          # A := HEAD built with layout-stable flags -> bench-results/bins/<sha>-perf/
# ... implement the step ...
dev/check.sh                                # build + make check + golden: must be ALL PASSED
dev/bench.sh run S --cases hmmsearch --reps 1 --baseline none      # real-data hashes
# kernel-rewrite steps only: shadow-execute old and new code, abort on any bit difference
HMMER_BUILD_DIR=$PWD/build-xchk CFLAGS="-O3 -Dp7_CROSSCHECK" dev/check.sh --quick
dev/perf/ab.sh <sha> build-xchk/src/hmmsearch --set all --check-only     # crosscheck over real data
# speed + A-vs-B output identity on the step's workload set(s), both sides layout-stable
B=$(dev/perf/snap.sh --perf)                # B := working tree, same flags
dev/perf/ab.sh <sha>-perf $B/hmmsearch --set <cpu|scale|query|milestone>
# record the numbers in PROGRESS.md and the commit message, then commit
```

Tool details are in [`dev/perf/README.md`](../perf/README.md).

**Why `--perf` builds.** The step 1.1 measurement showed that editing one file shifts every later function in the link by a few bytes. That changes the 64-byte alignment of the SSE kernels and moves cycles by 1-3% on workloads that never run the edited code. W5 got 1.1-1.6% "slower" with an identical instruction count. Step A/Bs therefore build both sides with `-O3 -falign-functions=64` (`snap.sh --perf`), which removes the artifact. Default `-O3` builds are compared once per phase, against the phase's start, to catch effects that only show in the shipped configuration.

**Accept a step when all of the following hold:**

1. **Correctness, zero tolerance.**
   - `dev/check.sh` passes.
   - The bench hashes match.
   - `ab.sh` prints `outputs identical`. A and B produce identical normalized `stdout`, `stderr`, `-o`, `--tblout`, `--domtblout` and `-A` files on every workload, and every rep matches its first rep.
   - Kernel steps also pass the crosscheck run.
2. **Speed on the target metric.**
   - `ab.sh` reports `faster` on the step's target workloads.
   - The metric is `cycles:u` for the `cpu` set (`--cpu 0`, deterministic work) and interleaved wall time for the other sets.
   - "Faster" means B/A is below 1 by more than both a floor (cycles 0.5%, wall 2%) and the measured spread (max−min)/median of either binary's reps.
   - Measured A/A spreads are in [`dev/perf/README.md`](../perf/README.md).
3. **No collateral slowdown.**
   - No other workload in the set is `SLOWER`.
   - Peak RSS is reported, and any growth is explained.
4. **A neutral step** is kept only if it's a declared prerequisite (a refactor). Otherwise it's reverted.

**At the end of each phase:**

- Run `dev/check.sh --bench M` against the dev baseline and merge the phase into `dev`.
- Re-save the baseline.
- Re-run the stage profile (`dev/profiling/profile_stage.sh` on human_pfam1k and hitters), because the bottleneck moves.

---

## 4. The steps

Legend:

- **BI**: why the step is byte-identical.
- **Measure**: `ab.sh` set and metric.
- **Expect**: estimate from the profile's stage shares.
- **cpu** = W1-W5 below at `--cpu 0`; **scale** = X1-X2; **query** = X3-X4; **milestone** = X5. They're defined in `dev/perf/workloads.tsv`.

### Phase 0 — measurement and safety net (no `src/` changes) · **done**

| Step | Change | Status |
|---|---|---|
| **0.1** | `dev/perf/snap.sh`: rebuilds, then copies the search binaries to `bench-results/bins/<sha>[-dirty-<diffhash>]/` with an `info.txt` | Done |
| **0.2** | `dev/perf/ab.sh A [B] --set …`: runs A and B interleaved (AB, BA, …). Uses `perf stat` (task-clock, instructions:u, cycles:u) plus GNU time for wall/user/sys/RSS. Compares normalized outputs A vs B and rep vs rep, then prints verdicts and writes JSON to `bench-results/ab/` | Done. A/A runs are identical with ratios at 1.000 ± spread; a negative test (B with `--seed 7`) is flagged as DIFF |
| **0.3** | 11 new golden cases from the unmodified binary: multi-query at `--cpu 1/3/4/8` (3 blocks, `-Z` unset, all tables, `--max -A`); `--restrictdb_stkey/_n` (SSI, threaded and `--cpu 0`); pressed `.h3m` query; `.gz` and stdin targets for a single query, plus the multi-query fatal errors | Done. 86 cases pass on 3 back-to-back runs |
| **0.4** | Workload sets in `dev/perf/workloads.tsv` | Done (table below) |

| Set | Workload | Why |
|---|---|---|
| cpu | **W1** `pfam_hitters × sprot_10k` (with `-A`) | Domain-stage heavy (D2 ≈ 22%) |
| | **W2** `pfam_len_mix × sprot_10k` | Typical filter mix |
| | **W3** `kern_M200 × sprot.fa` (full Swiss-Prot) | Bias filter 11.6%, D2 18% |
| | **W4** `pfam_len_med × sprot_long` | Bias filter 13.5% (long targets) |
| | **W5** `pfam_len_med × sprot_1k --max` | D2 41%, Forward/Backward heavy |
| scale | **X1** `pfam_100 × sprot_50k --noali` @ 4/8/N | Reader-bound scaling |
| | **X2** `kern_M200 × sprot.fa --noali` @ N | Single query, one full parse |
| query | **X3** `pfam_1k × sprot_1k` @ 0/N | Per-query overhead |
| | **X4** `Pfam-A.hmm` (text, and pressed `.h3m`) `× q_pick5.fa` @ 0/N | 22k queries × tiny DB |
| milestone | **X5** `pfam_1k × human --cut_ga` @ N | Realistic annotation (57 s today) |

### Phase 1 — bias filter (S2: 4-12% of CPU) · **done** (1.1 measured and dropped, 1.2 kept; see PROGRESS.md)

| Step | Change (files) | BI | Measure / Expect |
|---|---|---|---|
| **1.1** | Keep the filter DP matrix in `P7_BG`, grow-only, instead of `esl_hmx_Create`/`Destroy` on every `p7_bg_FilterScore` call (`src/p7_bg.c`: Create/Clone/Destroy/FilterScore) | Same arithmetic, different memory | cpu W3, W4 `cycles:u`. Expect a small gain (≤ 1%): this proves the pipeline end to end |
| **1.2** | Replace `esl_hmm_Forward` with a static 2-state streaming forward in `p7_bg.c` that needs no matrix. It must keep the exact operation order: `0.0 + dp0·t0k + dp1·t1k`, then `· eo`, divide by the row max, and store `(float) log(max)`. The float running sum must be accumulated in `i` order (identical to summing `sc[1..L+1]` afterwards). The old call stays under `p7_CROSSCHECK` | Same IEEE ops in the same order. No FMA on the SSE build. Crosscheck asserts bitwise equality | W3, W4 `cycles:u`. Expect S2 to drop ~2× (the `log()` stays), i.e. 2-6% of total CPU |

Out of scope, because it isn't byte-identical: a cheaper `log` (`logf`, or a periodic log of a running product).

Outcome:

- **1.1:** identical outputs; −0.15% to −0.5% instructions, below the cycles floor. Dropped, because 1.2 removes the matrix altogether.
- **1.2:** filter CPU went from 9.4% to 2.8% of W3. Cycles fell 1.2-4.3% (layout-stable builds) and 1.8-5.0% (default `-O3`) on every workload that runs the filter.
- **What's left:** the per-row `log()` and a division-bound dependency chain (row i needs row i−1's max). Both are inherent to doing the same operations in the same order.

**Candidate (unscheduled):** `-falign-functions=64` as a default build flag. The aligned base was 0.9-1.9% faster than the default base on W2-W4 and 0.5% slower on W5. It doesn't affect outputs. It needs its own A/B across all sets before it's proposed.

### Phase 2 — parse the target DB once (Bottleneck #1: 3.5-4× wall at scale)

| Step | Change (files) | BI | Measure / Expect |
|---|---|---|---|
| **2.1** | Refactor `thread_loop`/`serial_loop` (`src/hmmsearch.c`) behind a small "block source" interface (`next_block()`/`release_block()`). No behavior change | Pure refactor | cpu `instructions:u` within ±0.1%, X1 neutral. Kept as a prerequisite |
| **2.2** | Move the `esl_sq_Reuse` of block slots from the workers (`pipeline_thread`/`serial_loop`) to the block source, just before it refills a block. Workers then never write target sequences. Make the worker-side `ESL_SQ` pointer `const` where the API allows | Same inputs to the pipeline; the reset still happens before every refill | Neutral. Check X1: the reset now runs on the reader, which is the bottleneck thread (it's cheap, but verify) |
| **2.3** | **DB cache, opt-in.** Add a `--dbcache` option, off by default; while the step is being built, `HMMER_DBCACHE=1` forces it on for testing. Query 1 reads blocks from a growing pool and keeps them. Queries ≥ 2 enqueue the cached blocks and parse nothing. EOF and empty-block handshakes stay unchanged | Same `ESL_SQ` contents (name, acc, desc, dsq, n, idx) fed to the same pipeline. **Extra gate:** `HMMER_DBCACHE=1 dev/check.sh --only-regress` must match the *same* goldens | X1 @ 4/8/16, X3. Expect X1 @ 16: 9.8 s → ~2-3 s (scales past 4 threads). Record RSS on X2 |
| **2.4** | Compact cache storage. Shrink each cached `ESL_SQ` to fit, or pack into an arena: default `ESL_SQ` buffers cost ~1 KB per sequence, which is ~0.5 GB of overhead for Swiss-Prot | Storage only | RSS for Swiss-Prot multi-query: expect ≈ residues + headers (~0.3 GB) |
| **2.5** | **Memory cap, sized from the host.** Add `--dbcache_max <MB>`. Its default is computed at run time from the machine (a fraction of physical/available RAM via `sysconf`/`/proc/meminfo`). If the DB exceeds the cap, hmmsearch falls back to re-reading per query and prints a note. The cache stays **off unless `--dbcache` is given**. The frontend stays backward compatible: every existing command line behaves and prints exactly as before | Default command lines unchanged. With `--dbcache` the only output difference is the echoed option in the header | `scale` + `milestone` sets with `--dbcache`. Expect X5: 57 s → ~15-20 s |
| **2.6** *(intentional, new capability)* | With `--dbcache`, multi-query search works on `.gz` and stdin targets; this is a fatal error today. Without the flag the error is unchanged | Golden cases added for the new path; the existing error-case goldens stay as they are | Functional |

Decided 2026-10-04: the cache is opt-in, its cap defaults from the host, and `.gz`/stdin multi-query works only under `--dbcache`.

### Phase 3 — per-query fixed costs (many models × small DB)

| Step | Change | BI | Measure / Expect |
|---|---|---|---|
| **3.1** | Reuse per-thread `P7_PIPELINE` and `P7_TOPHITS` across queries. Add a reset that zeroes stats and keeps the grow-only matrices (`p7_pipeline.c`, `hmmsearch.c`) | The RNG reseeds per region. Stats are zeroed exactly as by `Create` | X3, X4 @ 16: wall and **sys** time (44 s sys today) |
| **3.2** | Reuse per-thread `P7_OPROFILE` buffers: copy into a preallocated, grow-only profile instead of `p7_oprofile_Clone`/`Destroy` | Same profile bytes | X3, X4 |
| **3.3** | Persistent worker threads. Create them once and park them on a condvar between queries. Removes per-query `pthread_create`/join and `impl_Init` | Same work, different scheduling. Results are independent of scheduling (§2) | X4 @ 16: expect 39 s → ~cpu0 level (~20 s text, ~4-5 s pressed). X3 |
| **3.4** | While workers search query q, the main thread reads and converts query q+1 (`p7_hmmfile_Read`, `ProfileConfig`, `oprofile_Convert`) | Output order unchanged | X3, X4, X5 |
| 3.5 *(measure first)* | Overlap query q's merge and output with query q+1's search, using double-buffered per-thread results. Do it only if the post-3.4 profile shows a barrier stall | Output order unchanged | X5 |

### Phase 4 — stochastic-traceback domain stage (D2: 12-41% of CPU on realistic and hit-rich runs)

Sub-profile of D2 on W1 (perf, `--cpu 1`, % of all CPU):

| Function | % of CPU |
|---|---|
| `esl_rnd_FChoose` | 5.2 |
| `p7_Null2_ByTrace` | 4.8 |
| `esl_vec_FNorm` | 3.7 |
| `link_spsamples` | 2.5 |
| `p7_StochasticTrace` + `select_*` | 2.2 |
| `p7_trace_Reverse` | 2.0 |
| `p7_trace_Append` | 1.4 |
| `p7_trace_Index` | 0.9 |
| `esl_cluster_SingleLinkage` | 0.4 |
| **Total** | **≈ 23** |

Every step must keep the RNG stream exactly the same: the same number and order of `esl_random()` draws, and the same comparisons.

| Step | Change (`src/impl_sse/stotrace.c` unless noted) | BI | Measure / Expect |
|---|---|---|---|
| **4.1** | `select_c/j/b/d/i`: replace `esl_vec_FNorm(path,2)` + `esl_rnd_FChoose(path,2)` with one static inline 2-way chooser. It does the same float Kahan sum (which is a plain sum for n=2), the same float divisions, the same double `norm`/`sum` and the same `roll < sum/norm` test. Crosscheck against the Easel pair | Identical op sequence | W1, W5 `cycles:u`. Expect 2-4% on W1 |
| **4.2** | `select_m`: the same for the 4-way choice, with the 4-term Kahan sum replicated exactly | Identical op sequence | W1, W5 |
| **4.3** | Remove the per-step `(k-1) % Q` and `/ Q` divisions. Track `(q, r)` incrementally as k decrements | Integer only | W1, W5 |
| **4.4** | Trace building: pre-grow the trace once per call (`p7_trace_GrowTo`), then use an unchecked append fast path in the traceback loop | Same trace contents | W1, W5. Expect ~1-2% |
| **4.5** | `p7_Null2_ByTrace` (`src/impl_sse/null2.c`): delete the insert-state term. Insert counts are never incremented (line 163 bins every k>0 state into the M slot), so `sv += +0.0` on `sv ≥ +0` is a bitwise no-op | Exact: x + (+0) = x for x ≥ +0 | W1, W5. Expect ~1% of W1 |
| **4.6** | Same function: hoist the per-q count load out of the K loop (q outer, x inner, in two halves of ~10 accumulators to avoid spills). The per-x accumulation stays in ascending q | Same per-lane add order | W1, W5. Keep only if faster |
| **4.7** | `p7_spensemble.c`: a specialized single-linkage with the link test inlined (no callback), producing **the same cluster numbering** as `esl_cluster_SingleLinkage`. Add a unit test against the Easel routine on random sample sets | Same assignments, checked by the test | W1, W5. Expect ~1.5% |

Out of scope, because it isn't byte-identical: fewer than 200 samples, a different RNG, or a deterministic replacement for sampling.

### Phase 5 — FASTA parsing speed (first pass, and single-query runs such as X2)

These steps change **Easel**, which is pinned in `dev/easel.ref` and lives outside this repo. **Decided:** Easel changes go to a fork of Easel in the user's GitHub. Each step commits on a branch there and bumps `dev/easel.ref` (url + commit).

| Step | Change (`easel/esl_sqio_ascii.c`) | BI | Measure / Expect |
|---|---|---|---|
| **5.1** | Track the file offset arithmetically and stop calling `ftello` (`lseek`) on every 4 KB read | Offsets are compared by an Easel utest | X2, cpu W2. `strace -c` shows no `lseek` |
| **5.2** | Raise `eslREADBUFSIZE` from 4 KB to 64-256 KB for sequence files | Buffer size only | X2 |
| **5.3** | `seebuf`/`addbuf`: scan lines with `memchr` and validate and digitize in one fused table-driven pass | A utest compares parsed blocks to the old parser on the bench FASTA files | X2. Expect reader 220 → ≥ 500 M res/s |
| **5.4** | `header_fasta`: `memchr`/`strpbrk` instead of byte-at-a-time `nextchar()` | Same name/desc strings | X2 |
| 5.5 *(optional, new dependency)* | In-process zlib or libdeflate for `.gz` instead of a `gzip -dc` pipe | Same bytes | `hmmsearch_in_gzip_sprot` (2.1× slower today) |

### Phase 6 — wider SIMD for the integer filters (S1 25-55%, S3 9-15% of worker CPU)

SSV/MSV (uint8/int8, saturating) and the Viterbi filter (int16) compute exact integer DP values. The lane count changes only the layout, so **bit-exact AVX2/AVX-512 versions are possible**, including identical overflow detection.

| Step | Change | BI | Measure / Expect |
|---|---|---|---|
| **6.0** | Build plumbing: configure check for AVX2/AVX-512BW, runtime dispatch (`__builtin_cpu_supports`), and an env override `HMMER_SIMD=sse` so one binary can be A/B-tested and golden-tested on both paths | No kernel change | Neutral |
| **6.1** | AVX2 SSV (`calc_band_*`) using a 32-lane striped byte profile, added to `P7_OPROFILE` at `p7_oprofile_Convert`. The crosscheck build compares the score against the SSE path | Integer-exact. Crosscheck plus golden under both `HMMER_SIMD` values | `k_msv_*` Mc/s, W2 `cycles:u`. Expect SSV up to ~1.7-2× |
| **6.2** | AVX2 full MSV fallback | Same | Kernel bench and W2 |
| **6.3** | AVX2 Viterbi filter (16 × int16) | Same | `k_vit_*`, W2 |
| **6.4** | AVX-512BW variants of 6.1-6.3. This box has them. Watch for clock-down | Same | Same as above |

Out of scope, because it isn't byte-identical: wider Forward/Backward and the posterior kernels. Lane count changes the float summation order (the D-chain serialization and the E-state horizontal sum). They need an explicit, reviewed golden update, which is a separate decision.

---

## 5. Order and expected end state

```
0.1-0.4  ->  1.1-1.2 (protocol warm-up)  ->  2.1-2.5  ->  3.1-3.4  ->  4.1-4.7  ->  6.0-6.4
                                              \-> 5.x can run in parallel once the Easel decision is made
```

| Workload | Today | After Phase 2-3 | After Phase 4 + 6 (CPU-bound by then) |
|---|---|---|---|
| X5 pfam_1k × human, `--cpu 16` | 57 s | ~15-20 s (≈ CPU/cores) | further ~1.3-1.6× |
| X1 p100 × sprot_50k, `--cpu 16` | 9.8 s | ~2-3 s | — |
| X4 Pfam-A × 5 seqs, `--cpu 16` | 39 s | ~20 s text / ~4-5 s pressed | — |
| Full Pfam-A × Swiss-Prot (projection) | ≥ 6.1 h reader floor | ~1.3-1.6 h on 8 cores | ~1 h or less |

These numbers come from the stage shares in `README.md` §4-§6. Each phase's milestone run replaces the estimates with measurements.

---

## 6. Decisions (made 2026-10-04)

1. **DB cache (2.3-2.5):** opt-in `--dbcache`, never on by default; the frontend stays backward compatible. The memory cap defaults to a value computed from the host.
2. **`.gz`/stdin multi-query (2.6):** yes, under `--dbcache` (it doesn't work today).
3. **Easel (Phase 5):** fork into the user's GitHub and pin through `dev/easel.ref`.
4. **Non-byte-identical speedups:** not for now. Every step in this plan must be byte-identical.

**Open:**

5. **Threaded determinism (§2).** Should threaded runs without `-Z` be made deterministic, equal to `--cpu 0`? This fixes an upstream flake in column widths. It's not a speedup, and it can change threaded outputs, so it would get its own golden review.

Track measured results per step in `dev/profiling/PROGRESS.md`. Use one row per step: step, commit, workload, metric, A → B, ratio, and the A/A spread.
