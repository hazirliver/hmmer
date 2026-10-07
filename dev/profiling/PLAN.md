# hmmsearch acceleration plan (small, verified steps)

Results, upstream vs every phase: [`RESULTS.md`](RESULTS.md). Re-profile after Phase 6 and proposed Phase 7: [`REPROFILE.md`](REPROFILE.md). Lessons learned so far: [`LEARNINGS.md`](LEARNINGS.md). Input: the profile in [`README.md`](README.md) (dev @ 8da3ab2, measured 2026-10-04).
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

### Phase 2 — parse the target DB once (Bottleneck #1: 3.5-4× wall at scale) · **done**

| Step | Change (files) | BI | Measure / Expect |
|---|---|---|---|
| **2.1** ✔ | **`--dbcache`, opt-in** (`src/hmmsearch.c`). The first query's pass keeps the blocks it reads (exactly what `ReadBlock` reads, `--restrictdb_*` included). Later queries are served from memory and parse nothing; no rewind is needed, so multi-query `.gz`/stdin works (was 2.6). Threaded mode keeps the work-queue protocol: pool blocks act as carriers whose `list` points at cached sequences while in flight. Workers skip `esl_sq_Reuse` on cached targets | Same `ESL_SQ` contents, same pipeline calls. **Extra gate:** `dev/check.sh --force-opt --dbcache` reruns all goldens with `--dbcache` on every hmmsearch call (84 cases byte-identical; the two multi-query `.gz`/stdin error cases are expected to differ). 4 new golden cases cover the option explicitly | Measured: X1 @ 4/8/16 = 0.50/0.33/0.28, X3 0.74-0.76. But single query × Swiss-Prot is 1.39× slower and RSS 114 → 653 MB (see 2.2) |
| **2.2** ✔ | Compact cache storage. The first pass reads into the usual reused block (exactly the default path) and appends a compact copy: the block, its `ESL_SQ` array, the strings and the `dsq`, all packed into 64 MB arena chunks that are 2 MB-aligned with `madvise(MADV_HUGEPAGE)`. This replaces fresh 1000-sequence blocks with Easel's default growable buffers. Workers reset the reused blocks in the first pass (per-pass `reuse_seqs`) and never touch cached targets | Same field values (struct copy, then exact-size copies of every buffer). The same gates as 2.1, plus AddressSanitizer runs | Measured: X2 (single query) 1.39 → 1.15× slower; RSS 653 → 547 MB, which is the data itself (206 MB residues, 120 MB headers, ~130 MB `ESL_SQ` structs). X1/X3 gains kept. Default path neutral (instructions 1.0000) |
| **2.3** ✔ | **Memory limit from the host.** `--dbcache_max <MB>`; the default is half of the memory available to the process (physical RAM, or a lower cgroup v1/v2 limit). If the cache outgrows the limit during the first pass, its copies are freed, a warning goes to stderr, and later queries re-read the file exactly as without the option. The current query is unaffected, because the first pass reads as the default path does (2.2). A non-rewindable target then fails with an explanatory message. Documented in `-h` and the man page | Default command lines unchanged. 2 new golden cases: over-limit fallback (hits = uncached golden) and the `.gz` error | Functional |

Decided 2026-10-04: the cache is opt-in, its cap defaults from the host, and `.gz`/stdin multi-query works only under `--dbcache`. The original 2.1 refactor and 2.2 reuse move weren't needed: a per-worker `reuse_seqs` flag leaves the default path's code unchanged.

### Phase 3 — per-query fixed costs (many models × small DB) · **done** (3.2 deferred; 3.3b, 3.3c, 3.4b added; see PROGRESS.md)

| Step | Change | BI | Measure / Expect |
|---|---|---|---|
| **3.1** | Reuse per-thread `P7_PIPELINE` and `P7_TOPHITS` across queries. Add a reset that zeroes stats and keeps the grow-only matrices (`p7_pipeline.c`, `hmmsearch.c`) | The RNG reseeds per region. Stats are zeroed exactly as by `Create` | X3, X4 @ 16: wall and **sys** time (44 s sys today) |
| ~~3.2~~ *(deferred: profile shows allocation < 1%; see PROGRESS)* | Reuse per-thread `P7_OPROFILE` buffers: copy into a preallocated, grow-only profile instead of `p7_oprofile_Clone`/`Destroy` | Same profile bytes | X3, X4 |
| **3.3** | Persistent worker threads. Create them once and park them on a condvar between queries. Removes per-query `pthread_create`/join and `impl_Init` | Same work, different scheduling. Results are independent of scheduling (§2) | X4 @ 16: expect 39 s → ~cpu0 level (~20 s text, ~4-5 s pressed). X3 |
| **3.3b** *(added)* | Start only as many workers per query as the previous pass had blocks, each woken on its own condvar. On a 5-target db, 15 of 16 workers used to wake every query just to collect an end-of-pass block (kernel runqueue spinlocks: 17% of cycles) | Which worker searches which target doesn't matter (§2) | X4 @ 16: sys time |
| **3.3c** *(added)* | With a complete `--dbcache`, split a small db into ~4 blocks per worker (≥ 20 targets each) by handing out slices of the cached blocks. `sprot_1k` is a single 1000-sequence block, so X3 @ 16 searched on one worker. Uncached passes keep full blocks: the serial reader sets their pace, and tiny blocks only added futex churn (X3 sys 0.1 → 21 s) | Same as 3.3b | X3 @ 16 with `--dbcache` |
| **3.4** | While workers search query q, the main thread reads and converts query q+1 (`p7_hmmfile_Read`, `ProfileConfig`, `oprofile_Convert`) | Output order unchanged | X3, X4, X5 |
| **3.4b** *(added)* | Several helpers prepare the next queries in parallel: reading stays serial and in file order (under a mutex), profile building runs concurrently into a ring of slots that the master takes in order. After 3.4 the single helper was the critical path on X4 | Same calls per query; output order unchanged | X4 @ 16 |
| **3.5** *(measured: barrier stall)* | Query-parallel mode with `--dbcache` and a small cached db (≤ 2M residues): each worker searches whole queries against all cached targets, in file order, with its own pipeline; the master prints jobs in query order (3.5a: output moved into `output_results()`). After 3.4b, X3 `--dbcache` @ 16 kept only ~4.7 of 16 cores busy: per-query hand-offs and the barrier cost ~1.4 ms for ~5 ms of CPU work | Each query is searched serially in file order, so its output equals `--cpu 0` even without `-Z`. Model errors (`--cut_ga`) and read errors are reported in order | X3, X4 `--dbcache`; X5 `--dbcache` must not regress |
| ~~3.5 (original)~~ | Overlap query q's merge and output with query q+1's search, using double-buffered per-thread results. Do it only if the post-3.4 profile shows a barrier stall | Output order unchanged | X5 |

### Phase 4 — stochastic-traceback domain stage (D2: 12-41% of CPU on realistic and hit-rich runs) · **done** (4.3 measured and dropped; 4.8 added and dropped; see PROGRESS.md)

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
| 4.8 *(added)* | `p7_spensemble.c`: skip the overlap divide when `nov <= 0` (exact for `n > 0`, `min_overlap > 0`) | Same predicate, checked by the utest | W5. Dropped: rarely fires, 2% slower |

Out of scope, because it isn't byte-identical: fewer than 200 samples, a different RNG, or a deterministic replacement for sampling.

### Phase 5 — FASTA parsing speed (first pass, and single-query runs such as X2) · **done** (5.1 measured and dropped; 5.0/5.0b test harness and 5.4b added; 5.5 not started; see PROGRESS.md)

These steps change **Easel**, which is pinned in `dev/easel.ref` and lives outside this repo. **Decided:** Easel changes go to a fork of Easel in the user's GitHub. Each step commits on a branch there and bumps `dev/easel.ref` (url + commit). Fork: `hazirliver/easel`, branch `hmmer-perf/phase5-fasta` (easel-0.49 + one commit per step).

| Step | Change (`easel/esl_sqio_ascii.c`) | BI | Measure / Expect |
|---|---|---|---|
| **5.0** ✔ *(added)* | Test harness: `esl_sqio_ascii_utest` includes upstream 0.49 `esl_sqio_ascii.c` verbatim as `ref_*` and compares both parsers call by call through the esl_sqio API, at read sizes down to 1 byte (test build only) | Test only | 6 of 6 mutants caught |
| **5.0b** ✔ *(added)* | Edge-case files (record/file ends) at every read size | Test only | Caught a mutant the random files missed |
| ~~5.1~~ *(dropped)* | Track the file offset arithmetically and stop calling `ftello` (`lseek`) on every 4 KB read | Offsets are compared by the utest | `lseek` 80,846 → 3 on X2, kernel cycles −9% on R1, but net −0.5% at best under ±3% loop-alignment noise; 5.2 cut the calls 4× anyway |
| **5.2** ✔ | Read files 16 KB at a time instead of 4 KB (measured; the plan's 64-256 KB was slower: the old two-pass parser wants the buffer in L1). Pipes and daemon format keep 4 KB, where results depend on where reads end | Read size only; the utest reads at 1 B-64 KB | R1 −4.8% user+kernel cycles |
| **5.3** ✔ | `seeaddbuf()`: validate and digitize in one pass, lines and record end found with `memchr`; exact hand-over to `seebuf`/`addbuf` where it can't stand in | Utest, incl. a mismatched alphabet that forces the hand-over | X2 wall −29%, X1 −34 to −40%, R1 −24% |
| **5.4** ✔ | `header_fasta`: run scans within the buffer (`header_run()`), `memcpy`, `nextchar()` only across buffers | Same strings, same allocation growth | R1 −6.1%, X1 −7 to −9% |
| **5.4b** ✔ *(added)* | Description end with `memchr('\n')`, then `\r`/`^A` before it (the byte loop was ~85% of `header_fasta()` after 5.4) | Same end position | R1 −11%, X2 −11%, X1 −11 to −17% |
| 5.5 *(optional, new dependency; not started)* | In-process zlib or libdeflate for `.gz` instead of a `gzip -dc` pipe | Same bytes | `hmmsearch_in_gzip_sprot` (2.1× slower today) |

Outcome (milestone, default `-O3` vs 62d1c72, outputs identical everywhere): X5 `pfam_1k × human --cut_ga` @ 16 uncached 48.2 → 25.7 s (−47%); X1 @ 4/8/16 −41 to −51%; X2 0.86 → 0.50 s; R1 −39% user+kernel cycles; W1-W4 −11 to −22% cycles; X3 −16 to −18%; `--dbcache` runs, X4 and W5 unchanged (little FASTA parsing). The reader alone (Swiss-Prot, 1000-sequence blocks, digital) is 2.6× faster, ~270 → ~700 M residues/s (cycles −61%, instructions −58%), with an identical hash over every parsed field on four bench files.

Found on the way (upstream behaviours, documented in the utest): on pipes, record offsets are -1 plus the position in the current read; `GuessAlphabet()` leaves the file mid-recording when it fails; daemon-format files read a byte past the buffer at record ends; `skip_whitespace()` indexes the input map with a negative `char` for bytes > 127; with an alphabet that disagrees with the format's input map (no standard alphabet does), which residues survive depends on where reads end, so 5.2 changed those (already garbage) results.

**Candidates:** ~~a 256-entry combined classify-and-digitize table in `seeaddbuf()`~~ done as **5.6** after Phase 6 (Easel 2be6bcc: R1 −10.8%, X1 @ 16 −18%, X5 uncached −17% wall; see PROGRESS.md); a faster name scan in `header_fasta()`; reads of 32-128 KB (~1% in the reader alone, not resolvable on R1 today).

### Phase 6 — wider SIMD for the integer filters (S1 25-55%, S3 9-15% of worker CPU) · **done** (6.3b, 6.4b, 6.4c added; AVX-512 used only for long models; see PROGRESS.md)

SSV/MSV (uint8/int8, saturating) and the Viterbi filter (int16) compute exact integer DP values. The lane count changes only the layout, so **bit-exact AVX2/AVX-512 versions are possible**, including identical overflow detection.

Profile at the start (f9dcb6b, pfam_1k × human `--cut_ga`, `--cpu 1`): SSV/MSV 38.6%, Viterbi filter 12.9% of all CPU (51.5% together; the reader 14.7%, down from 26.2% before Phase 5).

| Step | Change | BI | Measure / Outcome |
|---|---|---|---|
| **6.0** ✔ | Build plumbing: configure checks for AVX2 (`-mavx2`) and AVX-512BW, run-time choice `p7_simd_Level()` (`__builtin_cpu_supports`), env override `HMMER_SIMD=sse\|avx2\|avx512` (fatal on a bad value), `dev/check.sh --simd LEVEL`, `ab.sh --a-env/--b-env` | No kernel change | Neutral (instructions 0.9995-1.0000). Also fixed: `src/*.o` didn't depend on `impl_sse.h`; `check.sh`/`snap.sh` didn't reconfigure on `Makefile.in` changes |
| **6.1** ✔ | AVX2 SSV (`ssvfilter_avx2.c`) on a 32-lane copy of the SSV scores, made in `p7_oprofile_Convert()` (`p7_oprofile_SetSIMD()`, `simd.c`) | Integer-exact (argument in the file header). `simd_utest` bitwise vs SSE; crosscheck; goldens under every `HMMER_SIMD` | W2 −16.8%, W4 −16.1%, W3 −12.4% cycles; X5 −6.1% wall (`--dbcache` −11.3%) |
| **6.2** ✔ | AVX2 full MSV (`msvfilter_avx2.c`); row max kept in vector registers (a scalar round trip was 25% slower than SSE) | Same; plus an overflow-boundary scan in the utest | W3 −1.2% cycles |
| **6.3** ✔ | AVX2 Viterbi filter (`vitfilter_avx2.c`, 16 × int16) | Same; plus a deletion scan (M->D->M, D->D at every node) | W2 −7.8%, W4 −6.5%, W3 −4.4% cycles; X5 `--dbcache` −4.1% wall; X4 pressed serial +1.9% cycles (per-query copies) |
| **6.3b** ✔ *(added)* | The Viterbi filter's copies made on first use (lock-free, SSE fallback while another thread builds them) | Same; re-conversion and racing-clone utests; TSan | X4 pressed serial −1.8% cycles (back to the 6.2 level) |
| **6.4** ✔ | AVX-512BW variants of the three kernels, used only for models longer than one AVX2 SSV pass (M > 448; `p7_simd_LevelFor()`) | Same, for both wide levels | AVX-512 for every model was slower here (W3/W4 +3.5% cycles, X5 @ 16 +1.3/+2.6% wall, clock −1 to −1.6% at 16 threads); for M > 448: W2 −7.1% cycles, X5 neutral |
| **6.4b** ✔ *(added)* | Copy only the extra SSV vectors a band reads (a run-time vector width had turned 16 `memcpy`s per row into `rep movsq`) | Same | X4 pressed serial −2.8% cycles |
| **6.4c** ✔ *(added)* | SSE SSV for models of ≤ 2 wide vectors (AVX2: M ≤ 64), where the AVX2 cross-lane shift latency loses | Chooses between two exact kernels | R1 −6.7% cycles |

Outcome (milestone, default `-O3` vs f9dcb6b, outputs identical everywhere): W2 −29.5%, W4 −21.8%, W3 −17.2%, W1 −5.8% cycles; X5 `pfam_1k × human --cut_ga` @ 16 28.2 → 26.6 s, with `--dbcache` 14.4 → 12.0 s (−16%); X3 −14 to −18%; X1 −4 to −15%; X4, W5 unchanged; nothing slower. Instructions −26 to −52% on filter-heavy runs. Uncached X5 is held back by the serial reader at 16 threads. Profiles read from pressed files (hmmscan) get no wide copies and stay on SSE: per model and target sequence, restriping costs about what the AVX2 SSV saves.

Out of scope, because it isn't byte-identical: wider Forward/Backward and the posterior kernels. Lane count changes the float summation order (the D-chain serialization and the E-state horizontal sum). They need an explicit, reviewed golden update, which is a separate decision.

### Phase 7 — after the re-profile · **done** (7A.0-7A.4, 7B, 7C.1-7C.4, 7D.1; 7D.2 measured and dropped; see PROGRESS.md and [`REPROFILE.md`](REPROFILE.md))

The re-profile at 432ee21 found that at 16 threads the **end-of-query barrier** now costs more than any kernel. On X5 with `--dbcache`, workers are busy 61% of each pass, and 4.8 s of the 10.2 s run is per-query tail. Uncached, the reader is busy only 11.1 s of 18.9 s: it waits out each tail before it can rewind.

| Step | Change | BI | Measure / Outcome |
|---|---|---|---|
| **7B** ✔ | Exact `atof()` fast path in `read_asc30hmm()` (Clinger: ≤ 15 digits, ≤ 22 decimals, one division), `atof()` otherwise | Correctly rounded = glibc `strtod()`; utest bitwise on 1M strings + every token of Pfam-A and Dfam | Pfam-A text × 5 seqs: serial 20.6 → 9.5 s, @ 16 17.7 → 7.8 s |
| **7A.0** ✔ | 16 slices per worker with `--dbcache` when splitting can't change results (`-Z` or bit-score thresholds) | Which hits a worker keeps depends on its running Z only with E-values and no `-Z` | X5 `--dbcache` −7% |
| **7A.1** ✔ | Per-worker query slots sharing one DP workspace | Refactor | Neutral |
| **7A.2** ✔ | Up to 16 queries per pass over the targets (`MAXQ`; hit lists of a batch are held until printed) | Same targets and pipeline calls per query; output in query order | X5 uncached 0.43, X3 0.29, X1 @ 16 0.37 |
| **7A.3** ✔ | The same batching with `--cpu 0` | Same | W4 −8%, W2 −4%, X3 serial −8% |
| **7A.4** ✔ | Consecutive batches overlap ("streaming": two slot sets, tagged blocks); slots configured from a template | Scheduling only | X5 −11%, X3 −24%, X4 pressed −28% |
| fixes ✔ | Two error paths changed by printing headers after a pass (`--cut_ga` failing on the first model; a target-file parse error); goldens added, `dev/perf/errpaths.sh` | Same output as upstream, checked | — |
| **7C.1** ✔ | `p7_Pipeline()` split into `p7_pli_Filter()`, the parsers, `p7_pli_ForwardPass()`, `p7_pli_Domains()` | Refactor | Neutral |
| **7C.2** ✔ | Forward/Backward parsers on two targets at once (one per 128-bit half of an AVX2 register) | Per half, the SSE code's operations; bitwise utest, crosscheck | 1.4-1.5× Forward, 1.7-1.9× Backward per pair |
| **7C.3/7C.4** ✔ | hmmsearch pairs targets that pass the filters; `--cpu 0` keeps waiting targets in a small pool | Per target, `p7_Pipeline()`'s calls and state | W5 −9%, W1 −7%, W3 −4%; few targets reach the parsers on the realistic run |
| **7D.1** ✔ | Bias filter chain without `0.0 +` / `max(0, ·)` | Exact for values ≥ +0 (crosscheck vs `esl_hmm_Forward()`) | W3 −2%, W4 −3% |
| 7D.2 *(dropped)* | `choose()` with multiply + margin before the divide | Exact by a rounding argument | No gain |

Outcome (milestone, default `-O3` vs 7d637cb, outputs identical everywhere): X5 `pfam_1k × human --cut_ga` @ 16 20.8 → 7.85 s (2.6×), with `--dbcache` 11.5 → 7.6 s; X3 @ 16 0.216; X1 @ 16 0.384; Pfam-A text × 5 seqs @ 16 0.422, pressed 0.624; W1-W5 0.86-0.95 cycles. The default and `--dbcache` realistic runs are now equally fast: no re-parsing bound is left at 16 threads. Measured and dropped in the re-profile: an exact inline `log()` for the bias filter (exhaustively verified, −1%: the libm samples were stalls), LTO (−0.5 to −1%), the query-parallel mode for big dbs (1.53× slower on X5).

---

## 5. Order and expected end state

```
0.1-0.4  ->  1.1-1.2 (protocol warm-up)  ->  2.1-2.5  ->  3.1-3.4  ->  4.1-4.7  ->  6.0-6.4  ->  7B, 7A, 7C, 7D (all done)
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

5. **Whitespace-only differences (decided 2026-10-05):** acceptable. When results are 100% the same but a column is padded one character wider or narrower (the upstream threaded `-Z` flake, §2), the step still counts as correct. `ab.sh` reports these as `ws` (accepted with a note) instead of `DIFF`. No fix to upstream threaded determinism is planned. Golden regression cases stay strict, since they are deterministic.

Track measured results per step in `dev/profiling/PROGRESS.md`. Use one row per step: step, commit, workload, metric, A → B, ratio, and the A/A spread.
