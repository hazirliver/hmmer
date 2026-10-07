# hmmsearch re-profile after Phases 0-6 (dev @ 432ee21, Easel 2be6bcc)

The first profile ([`README.md`](README.md), 2026-10-04) found a serial reader, per-query fixed costs and 128-bit kernels. Phases 0-6 ([`PLAN.md`](PLAN.md), [`RESULTS.md`](RESULTS.md)) fixed those. This is the second pass (2026-10-06): the same stage profiles, plus instrumented 16-thread runs and some cheap experiments. Each experiment checks one hypothesis before it goes into a plan. Lessons from earlier phases ([`LEARNINGS.md`](LEARNINGS.md)) were applied throughout, and two of them changed conclusions here (§4).

Machine, data and build: as in README.md (16 vCPU = 8 cores × 2 HT, Sapphire Rapids, gcc 13.3 `-O3`). Profiling build: `build-prof6` (frame pointers, 92 goldens byte-identical).

**Outcome:** the proposals in §5 became Phase 7 ([`PLAN.md`](PLAN.md), step numbers in [`PROGRESS.md`](PROGRESS.md)). The estimates held: the realistic annotation run at 16 threads went 20.8 → 7.85 s by default (§5 estimated ~7-8 s) and 11.5 → 7.6 s with `--dbcache`; Pfam-A text × 5 seqs 18.1 → 7.7 s. Pairing targets in Forward/Backward (7C) gave the hit-rich and `--max` runs 7-9%, as estimated, but little on the realistic run, where few targets reach the parsers. [`RESULTS.md`](RESULTS.md) has upstream vs every phase.

---

## 1. TL;DR

| # | Finding | Evidence | What to do |
|---|---|---|---|
| 1 | **At 16 threads, the end-of-query barrier costs more than anything else.** Per-target cost is heavy-tailed: a target with hits goes through domain definition, which takes milliseconds. So every query waits for its slowest block while the other workers idle. | X5 `--dbcache` @ 16: workers busy only **61%** of each pass; summed over 1000 queries, **4.8 s of a 10.2 s run** is tail (first worker done → last worker done). Master serial time between queries is only 0.1 s. Uncached X5: the "serial reader" is busy only **11.1 s of 18.9 s**. It waits out each query's tail (6.7 s) before it can rewind. | **Batch several queries per database pass** (§5, 7A). Estimate: X5 uncached 19 → ~7-8 s, `--dbcache` 10 → ~7 s, and the reader floor at scale divided by the batch size |
| 2 | **Uncached runs are still reader-bound where the reader isn't hidden by the barrier.** | X1 (100 queries × Swiss-Prot 50k) @ 16: reading takes 2.35 s of the 2.87 s run. `seeaddbuf` is 7.2% of all CPU on the realistic run | The same batching: one parse serves K queries |
| 3 | **Forward/Backward are now the top functions, and they're latency-bound.** | `forward_engine` + `backward_engine` take 15.6% of the realistic run and **39.6%** of the hit-rich one (largest of any stage). `perf annotate`: about 2/3 of `forward_engine` is the serial D→D passes (one add→mul or mul chain per vector, Q long); the main DP loop is ~20%. Two copies on hyperthread siblings each keep 72% of their solo speed (SSV-heavy: 63%) | Interleave two targets per call, one per 128-bit half (§5, 7C). Exact, but needs the pipeline split per block. ≤ 1.4× on S4/S5 |
| 4 | **Text HMM parsing is 54% `atof()`.** | Pfam-A text × 5 seqs `--cpu 0`: `read_asc30hmm` 78% of the run, `atof` 54%. The serial read bounds the 16-thread run (15.4 s) | **Exact fast path for decimal tokens** (§5, 7B). Prototype: 0 mismatches on 230 M real tokens + 50 M random ones, 7.5× faster per token |
| 5 | The 6.5% of samples in libm, which the stage table filed under "pipeline glue", are the bias filter's per-residue `log()`. But they're **stalls on the filter's division chain, not `log()`'s own cost.** | An exact inline log (exhaustively verified on all 2^32 inputs) gained only 0.9-1.2% on W3/W4 | Dropped. LEARNINGS "profile shares can mislead", again |
| 6 | SSV/MSV and the Viterbi filter are still the largest CPU stages on typical runs (S1 35.9%, S3 10.4%), but they run at the vector ports' limit since Phase 6 | 71 G cells/s in situ, as in the kernel benchmark | Nothing exact left beyond Phase 6 |
| 7 | Build flags don't help | LTO: −0.5 to −1.3% cycles (within layout noise), outputs identical | Not worth a build-system change |

**Ranking for a Phase 7:** 7A (batched multi-query scheduler) is by far the largest lever on the realistic workloads, at every database size. 7B is small and safe with a large factor on one scenario. 7C has the largest CPU gain left on hit-rich runs, at high complexity. §5 has the details.

---

## 2. Where the CPU goes now (`--cpu 1`, % of samples)

Realistic annotation: `pfam_1k × human --cut_ga`. Wall time at `--cpu 1`: 109.4 s at upstream → 94.3 s at the start of Phase 6 → **75.3 s now**.

| Stage | upstream (8da3ab2) | start of Phase 6 (f9dcb6b) | **now (432ee21)** |
|---|---:|---:|---:|
| S1 SSV/MSV filter | 31.4 | 38.6 | **35.9** |
| R target read/parse (main thread) | 26.2 | 14.7 | **13.3** |
| D2 stochastic traceback clustering | 12.1 | 10.0 | **11.7** |
| S3 Viterbi filter | 9.6 | 12.9 | **10.4** |
| S4 Forward parser | 6.3 | 8.4 | **10.1** |
| "P pipeline glue": actually the bias filter's `log()` in libm, see §4 | ~2 | 5.6 | **6.9** |
| S5 Backward parser | 1.7 | 2.4 | **2.8** |
| S2 bias filter (own code) | 7.7 | 1.9 | **2.4** |
| D3-D7 envelope F/B, OA, decoding, null2 | 2.6 | 4.3 | **3.9** |

Hit-rich: `pfam_hitters × Swiss-Prot`, `-A`, tables. 40.3 s at the start of Phase 6 → **37.3 s now**. D2 18.7, S4 15.8, S1 13.7, S5 12.9, R 9.3, D3 10.0, S3 4.8, D6 4.0.

Top functions (self time):

| Function | realistic | hit-rich | Note |
|---|---:|---:|---|
| `calc_band_*` (SSE + AVX2 + AVX-512 SSV) | ~30 | ~11 | At the port limit (Phase 6) |
| `forward_engine` | **11.6** | **21.2** | Float, 4 lanes; latency-bound DD passes (§3) |
| `p7_ViterbiFilter_avx2` + `_avx512` | 10.3 | 4.7 | |
| `seeaddbuf` + `header_fasta` (main thread) | 7.9 | 5.3 | Serial at 16 threads |
| libm `log` (from `filter_forward`) | 6.5 | ~1.5 | Stall sink, §4 |
| `backward_engine` | 4.0 | **18.4** | Same structure as Forward |
| `p7_StochasticTrace` | 4.1 | 8.0 | Hottest instruction: `comisd` waiting for `roll < sum/norm`'s `divsd` |
| `p7_spensemble_Cluster` | 3.4 | 2.2 | `divss` of the overlap test |
| `p7_Null2_ByTrace` | 2.5 | 4.4 | |
| `p7_bg_FilterScore` | 2.4 | 1.0 | |
| `p7_OptimalAccuracy` | 1.0 | 4.2 | |

`stages.py` was updated for the AVX kernels and the new parser functions. Caveat: libm has no frame pointers, so a libm leaf loses its direct caller and is charged to the caller's caller. `--call-graph dwarf` on a 60-model run showed that the libm samples come from `filter_forward()` (the bias filter).

---

## 3. What bounds the 16-thread runs now

A throwaway build of `hmmsearch.c` (patch kept out of the tree) timestamped every pass: start, each worker's end, the master's wait, time inside `esl_sqio_ReadBlock`, and time blocked on a full queue.

| Run (@ 16 unless noted) | wall | worker busy / (workers × pass time) | tail: first worker done → last, summed | master serial between passes | reader inside `ReadBlock` |
|---|---:|---:|---:|---:|---:|
| X5 `--dbcache` | 10.2 s | **61%** | **4.8 s** | 0.11 s | — |
| X5 `--dbcache`, `--cpu 8` | 11.6 s | 72% | 4.3 s | 0.07 s | — |
| X5 uncached | 18.9 s | (72%; includes waiting for blocks) | 6.7 s | 0.12 s | **11.1 s** |
| X5 uncached, `--cpu 8` | 17.9 s | | 5.4 s | 0.09 s | 10.4 s |
| X1 `--dbcache` (100 queries × 50k) | 1.38 s | 73% | 0.48 s | 0.01 s | — |
| X1 uncached | 2.87 s | | 0.41 s | 0.02 s | **2.35 s** |

Reading:

- With `--dbcache`, the per-query tail takes about half of the realistic run. That's not the master: its serial share is 1%. It's that every query ends with a few workers finishing their last, slowest slices. Hits are rare, but each costs milliseconds (stochastic traceback: 3-6 ms per multi-domain region, README §3).
- Finer slices help only partly (one-line change to `thread_loop()`, same outputs):

  | slices per worker | 4 (now) | 8 | 16 |
  |---|---:|---:|---:|
  | X5 `--dbcache` @ 16 wall | 10.0 s | 9.7 s | **9.3 s** |
  | tail | 4.7 s | 3.9 s | 3.3 s |

  So most of the tail is a few expensive targets, not slice granularity.
- Without the cache, the reader works 11.1 s of 18.9 s and waits the rest: it can't rewind for query q+1 before query q's tail is done. A run without the barrier would be reader-bound at ~11-12 s. Parsing each block once for several queries would also shrink those 11 s.
- **The existing query-parallel mode is not the answer for big databases.** Forced on X5 (`QP_MAXRES` raised in a throwaway build), it was **1.53× slower** (16.5 vs 10.8 s, outputs identical). Results are printed in query order, so one long model blocks the 32-job window while the other workers run dry. It stays right for small databases (X3/X4).
- Compute floor for X5 at 16 threads: the busy time at `--cpu 8` (66 s on 8 cores) divided by 8 cores × ~1.27 (hyperthread gain measured below) ≈ **6.5 s**. Today: 18.9 s uncached, 10.2 s cached.

---

## 4. Experiments in this pass (all outputs byte-identical)

| Experiment | Result | Verdict |
|---|---|---|
| **Exact inline `log()` in the bias filter.** fdlibm-style double log (< 1 ulp), rounded to float only when the double is ≥ 4096 double-ulps from a float midpoint, else libm (fallback rate 1.5e-5). | **0 mismatches against glibc over all 2^32 float inputs.** In place: W3 −0.9%, W4 −1.2%, W2 0.0% cycles (`ab.sh` vs `-O3` dev; W5, which never runs the filter, −0.5%: layout) | **Dropped.** The 6.5% libm share was stall samples: the filter's critical path is the per-residue `mul → add → add → mul → max → div` chain (~30 cycles), and `log` runs in its shadow. Removing the always-redundant `0.0 +` and `max(0, ·)` from that chain is exact for non-negative values but worth ~1% at most |
| **Exact `atof()` fast path** (Clinger: `[+-]digits[.digits]`, ≤ 15 significant and ≤ 22 fraction digits → `(double) m / 10^e`, one correctly rounded division = glibc `strtod`; everything else → `atof`) | 0 mismatches on all 230.5 M whitespace tokens of Pfam-A.hmm + dfam_human.hmm and 50 M random decimal strings; 7 vs 53 ns per value | **Candidate 7B** |
| **LTO** (`-O3 -flto=auto`, `gcc-ar`) | 92 goldens identical; cycles W1-W5 0.991-1.001, R1 0.987 | Not worth it. AVX2/AVX-512 flags don't enable FMA, so LTO can't introduce contraction; checked with `gcc -dM -E` |
| **Hyperthread probe** (two `--cpu 0` copies pinned to sibling CPUs vs one alone) | Per-copy slowdown: W5 `--max` (F/B, traceback) 1.38×, W1 hit-rich 1.49×, W3 SSV-heavy 1.58× | F/B and traceback are more latency-bound than the integer filters. Interleaving two independent targets in one thread can gain up to ~1.4× there |
| Forced query-parallel mode on X5 | 1.53× slower | §3 |
| 8 / 16 slices per worker with `--dbcache` | X5 −3% / −7% | Cheap, but subsumed by 7A |

Not tried, and why:

- **PGO**: gains would be in branchy code (parser, traceback), and it needs training runs in the build. Same category as LTO.
- **Wider float SIMD for Forward/Backward**: it changes the summation order, so it isn't byte-identical (decision 4). It remains the largest non-identical lever: 2-4× on 16-40% of CPU.

---

## 5. Proposed Phase 7 (ranked)

### 7A. Batched multi-query scheduler: K queries per database pass, no per-query barrier

What changes (`src/hmmsearch.c`, threaded path only):

- Each worker keeps a small array of per-query slots: a profile clone, a `P7_BG` (its filter HMM is per model), pipeline counters and a hit list. The DP scratch matrices can stay one per worker.
- A work item is (block, set of active queries). A block read or sliced once is searched by every query in the batch, then released (refcount).
- Uncached: the reader parses each block once per batch, not once per query.
- With `--dbcache`: each slice is searched by K queries, so one query's slow tail overlaps with the others' work.
- Queries are merged and printed in file order when their last block is done. The next batch starts while the current one drains, so no worker waits at a barrier.

Why it's byte-identical: each query sees the same targets, parsed by the same code, through the same `p7_Pipeline()` calls with its own profile, null model and counters. Hits are sorted by key, name and position, and domain definition reseeds per region (PLAN §2), so neither assignment nor order matters. Targets shared between queries are read-only: the `--dbcache` path's `reuse_seqs = FALSE` already covers this.

Behaviour that must not change:

- Running Z without `-Z`: per (worker, query), exactly as today. It's the accepted `ws` class; check how often it shows, as after step 3.3c.
- The multi-query `.gz`/stdin error must stay where it is.
- Read errors are reported after the preceding query's output.
- Memory: K slots per worker.

Expected (from §3):

- X5 uncached @ 16: 18.9 s → **~7-8 s** (the compute floor; the reader drops to ~11 s / K).
- X5 `--dbcache`: 10.2 s → **~7 s**.
- X1 uncached: 2.9 s → ~1.2 s.
- Full Pfam × databases too big to cache: the reader floor (README §1, now ~1 G residues/s) is divided by K. That's the at-scale case, where `--dbcache` can't help.

Steps:

| Step | Change |
|---|---|
| 7A.0 | Optional now: 16 slices per worker with `--dbcache` (one line, −7% on X5 `--dbcache`) |
| 7A.1 | Per-worker query slots, K = 1. A refactor, must be neutral (instructions 1.0000) |
| 7A.2 | `--dbcache`: slices searched by K queries, batches overlapped |
| 7A.3 | Uncached: blocks shared by K queries, refcounted |
| 7A.4 | Choose K from memory and model size; fold the query-parallel mode in for small databases |

Checks: TSan, ASan, stress loop, error paths (LEARNINGS "Threading changes"), all goldens with `--cpu` 1-16, and the `ws` frequency count.

### 7B. Exact fast path for `atof()` in the text HMM reader

`p7_hmmfile.c` `read_asc30hmm()` reads ~50 numbers per node with `atof()`, as `expf(-1.0 * atof(tok))` or plain `atof(tok)`. A static `fast_atof()` with the prototype's rule replaces those calls; other tokens still go to `atof()`.

- Byte-identical by construction, and checked on every token of Pfam-A and Dfam.
- Unit test: the same comparison over the test HMM files plus random and edge-case strings (`-0.0`, `1.`, `.5`, 15/16/17 significant digits, 22/23 fraction digits, `1e5`, `inf`, trailing junk).
- Mutation test: drop the 15-digit limit, so the division is no longer exact.
- Expected: Pfam-A text × 5 seqs serial 19.2 s → ~10 s, @ 16 15.4 s → ~5 s (the read is the serial path). Text HMMs read by any program (hmmpress, hmmstat, hmmconvert, hmmfetch) get the same gain. Realistic runs: none (models are a small part of their time).

### 7C. Forward/Backward on two targets at once

The parser kernels are SSE `__m128` code. In AVX2 (no FMA), a `__m256` op does the same op on two independent 128-bit halves. So target A in the low half and target B in the high half run each target's exact SSE operation sequence:

- in-lane shifts (`vpslldq`) and in-lane shuffles for the horizontal sums;
- per-half scalar special states;
- rescaling by `(1/xE or 1.0)` per half (`x * 1.0 == x`);
- a finished half's DD passes masked by zeroing its `dcv` (`x + 0 == x` for the non-negative values here).

The two DD chains then overlap. Needs `p7_Pipeline()` split per block into a filter stage (MSV, bias, VF) and a post-Forward stage, so Forward/Backward candidates can be paired.

- Expected ≤ 1.4× on S4/S5 (the HT probe): realistic −3 to −4%, hit-rich −8 to −10% CPU.
- High complexity; do it after 7A, which makes the realistic runs CPU-bound.

### 7D. Small exact steps (each ~1-2%, only if cheap to verify)

- Bias filter chain: drop `0.0 +` and `max(0, ·)` (values are ≥ +0).
- `choose()`: decide `roll < sum/norm` from `roll * norm` with a relative margin, and divide only when it's too close.
- `p7_spensemble` overlap test: the same trick for its `divss`.

---

## 6. Scaling out: more threads, or several processes?

**Measured: query sharding across processes on the same 16 vCPUs.** `reprofile/shard_queries.sh` splits `pfam_1k.hmm` into P contiguous shards of equal total model length and runs P processes with `--cpu T` at once. X5, two rounds, interleaved; the rounds agreed within 0.5 s:

| processes × threads | uncached | `--dbcache` |
|---|---:|---:|
| 1 × 16 (today) | 22.1 s | 12.0 s |
| 2 × 8 | 13.2 s | 9.2 s |
| 4 × 4 | 9.6 s | **8.2 s** |
| 8 × 2 | **9.3 s** | 8.5 s |
| 16 × 1 (`--cpu 0`) | 10.8 s | 9.7 s |

The shards' outputs concatenated in order equal the single run byte for byte: tables, domain tables and the main output (per-query sections, timing lines masked), in all 10 configurations. Each process searches every target, so Z, domZ and the hit lists of a query don't depend on the split.

Why it's faster: separate processes are separate barriers. While one waits for its query's slowest target, the others keep the cores busy, and without the cache each process has its own reader. That's the effect 7A aims for inside one process. These numbers back its estimate of ~7-8 s; 7A would also balance load dynamically instead of by static shards.

**Splitting the target database is not exact with text outputs.** `domZ` defaults to the number of targets reported across the whole search (`p7_tophits.c:967`), so each shard would compute its own domain E-values (and E-value-based domain reporting). Without `-Z`, each shard's Z is its own target count. Exact target sharding has to merge the raw hit lists and pipeline counters before `p7_tophits_Threshold()`, as upstream `--mpi` does. It's needed only when one node can't hold or read the database fast enough.

**Beyond 16 threads in one process** (not measurable on this 16-vCPU VM; inferred from §3):

- The per-query tail is a fixed cost, so efficiency falls as threads grow: 72% busy at 8 workers, 61% at 16.
- Without the cache, a query has only as many blocks as `n_targets / 1000`. Human has 20, so no more than 20 workers can ever be busy on a query.
- The single reader parses about 1 G residues/s.
- A single mutex work queue with a master handing out slices one by one; probably fine at 16, untested at 100+.

7A removes the first three. On multi-socket hosts, one process per socket (query sharding, `numactl`) will likely stay useful anyway: the `--dbcache` arena is first-touched on one NUMA node.

**A shared database cache across processes:** a pre-digitized, mmap-able target file. ESL_SQ views would point into the mapping: NUL-terminated names and descriptions, `dsq` with its sentinels. That gives one copy in the page cache for any number of processes and no parsing at all, even on the first pass.

- Byte-identical by the `--dbcache` argument: the same `ESL_SQ` contents, read-only targets, the same `--force-opt` gate.
- It matters when the database is large (UniRef-size: P per-process `--dbcache` copies won't fit) or the run is short (the first parse dominates).
- For human-size databases, per-process `--dbcache` already does the job (4 × 4 above).

---

## 7. Reproduce

```
HMMER_BUILD_DIR=$PWD/build-prof6 CFLAGS="-O3 -g -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer" dev/check.sh --quick
D=~/.cache/hmmer-bench/derived-b36ad365c6d1; export HMMSEARCH=$PWD/build-prof6/src/hmmsearch OUT=$PWD/bench-results/profiling/p7
dev/profiling/profile_stage.sh human_pfam1k 1 --cut_ga --domtblout h.dom $D/pfam_1k.hmm $D/human.fa
dev/profiling/profile_stage.sh hitters_sprot 1 --tblout hs.tbl --domtblout hs.dom -A hs.sto $D/pfam_hitters.hmm $D/sprot.fa
# per-thread cycles at 16 threads
perf record -F 499 -o x5.data $HMMSEARCH --cpu 16 [--dbcache] --cut_ga --tblout /dev/null --domtblout /dev/null -o /dev/null $D/pfam_1k.hmm $D/human.fa
perf script -i x5.data -F pid,tid,period | awk '...sum period by tid...'
# who calls libm (frame pointers lose it): DWARF unwinding on a small run
perf record --call-graph dwarf,16384 -F 499 $HMMSEARCH --cpu 0 --cut_ga -o /dev/null p60.hmm $D/human.fa
```

Throwaway experiments, kept in [`reprofile/`](reprofile/) (never committed to `src/`):

- `barrier-stats.patch`: per-pass timestamps and reader time (§3). Prints `XSTAT` lines to stderr. Use `git apply`, build in a scratch build dir, then `git checkout src/hmmsearch.c`.
- `bias-fastlog.patch` and `fastlog.c`: the exact inline log (`gcc -O3 fastlog.c -lm && ./fastlog` checks all 2^32 inputs in about a minute).
- `fastatof.c`: the exact `atof()` fast path (`./fastatof Pfam-A.hmm ...` compares every token).
- `shard_queries.sh`: the processes × threads runs (§6).
- The forced query-parallel and slice variants: `QP_MAXRES` and the `4 * gate->nworkers` split in `thread_loop()`, edited by hand.

Raw profiles are in `bench-results/profiling/p7/` (gitignored).
