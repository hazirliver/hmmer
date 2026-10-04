# hmmsearch end-to-end profile (HMMER 3.4, `dev` @ 8da3ab2)

This is a step-by-step breakdown of where an `hmmsearch` run spends its time on the protein datasets of the benchmark (`dev/bench`). Each step is marked as **expected** (running at the speed its design allows), **inefficient** (slower than it needs to be at the same algorithmic cost) or **bottleneck** (limits wall time or scaling at production scale). Nothing in `src/` or `easel/` was changed. This directory holds only the measurement tooling and this report.

Machine: 16 vCPU Xeon Platinum 8468 KVM guest. That is **8 physical cores × 2 hyperthreads**, ~3.2 GHz under load, 196 GB RAM, gcc 13.3 at `-O3`, SSE code path (HMMER 3.4 has no AVX/AVX-512 kernels).
Data: `~/.cache/hmmer-bench/derived-b36ad365c6d1`. Swiss-Prot is 570,420 seqs and 206.3 M residues; Pfam-A 37.0 is 21,979 models with ΣM = 3.59 M (mean M ≈ 163).

---

## 1. TL;DR — what limits hmmsearch at scale

| # | Finding | Evidence | Verdict |
|---|---|---|---|
| 1 | **A single reader thread re-parses the whole target FASTA for every query model.** Above about 2-3 worker threads, wall time is set by this serial parse, not by the search. | The reader runs at 98% of wall at `--cpu 16` (p100×50k). CPU utilization plateaus at 2.9× whatever the thread count. Wall time is flat from 4 to 16 threads and gets *worse* at 12-16. Full Swiss-Prot parse is ~1.0 s per query (~220 M res/s); 22k Pfam queries need ≥ 6 h of serial parsing. | **Bottleneck #1** |
| 2 | Per-query fork/join. Every query creates and joins N pthreads, clones the profile and allocates a pipeline per thread, and reads the next HMM as text. | Full Pfam-A against 3 seqs: 20.2 s at `--cpu 0`, 39.0 s at `--cpu 16` (44 s sys). 69% of the serial overhead is text HMM parsing (`read_asc30hmm`, 0.72 ms/model). A pressed `.h3m` brings it to 4.3 s. | **Bottleneck** for many models × small DB (genome annotation) |
| 3 | SSV/MSV filter, the dominant compute stage (25-55% of worker CPU). | 35-45 G cells/s at M ≥ 200 and 22 G cells/s at M = 50. IPC 3.6-3.8, L1-resident. | **Expected** (at the SSE peak). The biggest lever is wider SIMD (AVX2 2×, AVX-512 4× lanes) |
| 4 | Stochastic-traceback domain clustering (200 sampled traces per multi-domain region). | 12% of all CPU on Pfam-1k × human `--cut_ga`, 17-24% on long models or hit-rich sets, 41% with `--max`. Costs 14-30 µs per sample and 3-6 ms per region, all scalar. | **Inefficient / bottleneck** on hit-rich workloads |
| 5 | Bias filter (`p7_bg_FilterScore` → `esl_hmm_Forward`). | ~20 ns/residue (~50 M res/s). Each call mallocs an L×2 matrix and does a `log()` per residue. Per residue it costs as much as an MSV pass of an M ≈ 700-900 model. Takes 4-12% of CPU. | **Inefficient** |
| 6 | Input formats. gzip goes through a single `gzip -dc` pipe (2.1× slower). EMBL/UniProt text is 3× slower. 4 KB `fread`s each come with an `lseek` (`ftello`). gzip and stdin can't be used with more than one query at all. | §5 | **Inefficient**, and a functional limit |
| 7 | Viterbi filter, Forward/Backward parsers, envelope DP, OA, output, memory. | §3 | **Expected** for 128-bit SIMD. Memory is not a concern (≤ 335 MB at 16 threads) |

**Rule of thumb from the measurements.** Per query and per target residue, the reader costs about 1/220e6 s, serial. The workers cost about M̄/18e9 s of CPU, end to end. The search is therefore **reader-bound as soon as the worker core count exceeds M̄/82**. For Pfam (M̄ ≈ 163) that's about 2 cores, which matches the observed plateau at ~2.9× utilization (1 reader + ~2 busy workers). At today's speeds, hmmsearch can't use more than about 3 cores on Pfam-like query sets, however large the machine.

**Projection: full Pfam-A × Swiss-Prot on this box.**

- Compute: ΣM×R = 7.4e14 cells at ~16-20 G cells/s per core is ~10-13 core-hours, about 1.3-1.6 h on 8 cores if it scaled.
- Reader floor: 21,979 × ~1.0 s ≈ **6.1 h wall**, whatever `--cpu` is.

Everything scales linearly with database size. For a database N× Swiss-Prot, multiply both numbers by N.

---

## 2. What an end-to-end run does (and where the code is)

```
main thread (serial)                                     worker threads (--cpu N; per query!)
-------------------------------------------------------  ------------------------------------------------
open seq db, open HMM file
for each query HMM:                                       
  [Q] p7_hmmfile_Read (text parse of HMM)                 
  [Q] p7_ProfileConfig, p7_oprofile_Convert               
  [Q] per thread: p7_oprofile_Clone, p7_pipeline_Create,  
      esl_threads_AddThread -> pthread_create             
  rewind db (esl_sqfile_Position 0)                       
  [R] thread_loop: esl_sqio_ReadBlock (1000 seqs/block,   pipeline_thread: for each seq in block
      queue depth 2N) -- parses+digitizes FASTA again        [P] p7_pli_NewSeq, p7_bg_SetLength, ReconfigLength
                                                             [S0] p7_bg_NullOne
                                                             [S1] p7_MSVFilter (-> p7_SSVFilter, calc_band_*)   P<0.02
                                                             [S2] p7_bg_FilterScore (bias)                     P<0.02
                                                             [S3] p7_ViterbiFilter                             P<1e-3
                                                             [S4] p7_ForwardParser                             P<1e-5
                                                             [S5] p7_BackwardParser
                                                             [D]  p7_domaindef_ByPosteriorHeuristics:
                                                                  D1 p7_DomainDecoding (B/E posteriors)
                                                                  per region: if multidomain ->
                                                                  D2 p7_Forward(region) + 200x p7_StochasticTrace,
                                                                     p7_Null2_ByTrace, p7_spensemble_Cluster
                                                                  per envelope (rescore_isolated_domain):
                                                                  D3 p7_Forward + p7_Backward (envelope)
                                                                  D4 p7_Decoding   D5 p7_Null2_ByExpectation
                                                                  D6 p7_OptimalAccuracy + p7_OATrace
                                                                  D7 p7_alidisplay_Create
  pthread_join all (barrier)                              
  [O] p7_tophits_Merge/Sort/Threshold, Targets, Domains,  
      tabular outputs, -A alignment, p7_pli_Statistics    
```

Files: `src/hmmsearch.c` (main loop at :462, `thread_loop` :1321, `pipeline_thread` :1376), `src/p7_pipeline.c:697` (`p7_Pipeline`), `src/p7_domaindef.c` (domain definition), `src/impl_sse/*` (SIMD kernels), `easel/esl_sqio_ascii.c` (FASTA reader: `seebuf` :2250, `addbuf` :2336, `header_fasta` :3009, 4 KB `loadmem` :2057), `src/p7_bg.c` (`p7_bg_FilterScore`), `easel/esl_hmm.c` (`esl_hmm_Forward`), `easel/esl_threads.c` (`pthread_create` per `AddThread`).

---

## 3. Every step, annotated

Throughputs are in situ, from `--cpu 1` runs with perf at 997 Hz. Cell counts come from bpftrace uprobes (`cells.sh`) on single-model × full Swiss-Prot runs. One "cell" is one model state × one residue, M × L.

| Step | Measured cost | Share of CPU (typical → extreme) | Verdict | Notes / acceleration direction (not implemented) |
|---|---|---|---|---|
| **R: target DB read + parse + digitize** (main thread) | ~220-245 M res/s, i.e. ~1.0 s per full Swiss-Prot pass. 331 MB/s of FASTA | 15-35% of total CPU (54% at M=50). **~98-100% of wall at ≥ 4 threads** | **Bottleneck #1** | Repeated **once per query**. Serial. Two scalar byte-at-a-time passes (`seebuf` validates, `addbuf` digitizes through `inmap[]`). Headers go one `nextchar()` per byte. 4 KB `fread` + `lseek` per block (7,371 reads + 7,360 lseeks for 30 MB). Options: parse once and keep a digitized DB in RAM across queries (`src/cachedb.c` / hmmpgmd already does this); batch several queries per DB pass; parallel or SIMD parsing; a binary pre-digitized DB (mmap) |
| R-gz: `.gz` input | `gzip -dc` child process costs 2× the parse itself | full Swiss-Prot: 2.05 s vs 0.99 s wall | **Inefficient** | Single-threaded external gunzip through popen. In-process zlib/isa-l/libdeflate or pigz would help. **Multi-query on .gz/stdin is refused** ("isn't rewindable") |
| R-embl: UniProt `.dat` | 52 M res/s (3.0-3.6× slower than FASTA per residue) | 51% of CPU on 1 model | **Inefficient** (12× more bytes than FASTA) | Convert to FASTA or a binary DB up front |
| **Q: per-query model setup** (main thread) | Text HMM parse 0.72 ms/model (`read_asc30hmm` + libc `strtod`); with a pressed `.h3m` 0.19 ms/model total | 1-2% of typical runs. **100% when the DB is tiny** | **Bottleneck** for many-model × small-DB runs | Pfam-A × 3 seqs: 20.2 s serial (text) → 4.3 s (pressed). Use `hmmpress`'ed input, or pre-convert oprofiles once |
| **Q/T: thread fork/join per query** | ~0.8-0.9 ms/query extra at `--cpu 16` + 44 s sys over 22k queries | 0 for single-query runs. Dominant for 22k-query × small-DB | **Bottleneck** (grows with threads) | `esl_threads_AddThread` calls `pthread_create` per query and `WaitForFinish` joins. Plus per-thread `p7_oprofile_Clone`, `p7_pipeline_Create` (omx alloc + page faults) and `impl_Init`. pfam_1k × sprot_1k: 5.1 s at cpu0 → 6.25 s at cpu16. A persistent pool with reused per-thread pipelines would fix it. The end-of-query barrier also idles workers while main does output and setup |
| T: work-queue sync | `esl_workqueue_*` < 0.3%, 1000-seq blocks | negligible | Expected | The mutex/condvar queue isn't contended. The coarse 1000-seq blocks are fine |
| P: per-seq setup (`NewSeq`, `bg_SetLength`, `oprofile_ReconfigLength`, `omx_GrowTo`) + pipeline glue (`esl_gumbel_surv`, `exp`) | ~0.3-1 µs/seq | 1-3% (≈10% with ≤60-aa sequences) | Expected / minor | Only matters for very short targets |
| S0: null1 (`p7_bg_NullOne`) | O(1) per seq | < 0.1% | Expected | — |
| **S1: SSV/MSV filter** (`p7_SSVFilter` → `calc_band_*`; full `p7_MSVFilter` fallback is rare) | **35.0 G cells/s (M=205), 45.2 (M=397), 41.7 (M=992), 21.8 (M=48)**. IPC 3.6-3.8, ~14 cells/cycle | 25-55% of worker CPU (31% of all CPU on Pfam-1k × human) | **Expected** (at the 128-bit SIMD peak; matches the kernel microbench 46.7 G cells/s) | Not memory-bound. The only lever is wider vectors (AVX2 32 lanes, AVX-512 64 lanes of uint8 → up to 2-4×) or a smarter prefilter. Short models (M < 64) lose ~half their speed to per-row overhead (`Q` ≥ 2 vectors minimum, horizontal max per row). `-march=native` gives nothing (§7) |
| **S2: bias filter** (`p7_bg_FilterScore` → `esl_hmm_Forward`) | **~50 M res/s ≈ 20 ns (~64 cycles) per residue**, 8 µs/call | 3-12% (7.7% on Pfam-1k × human, 11.6% at M=205, 13.5% on long seqs) | **Inefficient** | Runs on every MSV survivor (~2-11%), so cost is proportional to the survivors' length. It's a 2-state HMM, scalar, with a `malloc` of an L×M matrix (`esl_hmx_Create`) per call, a division per state and a `log()` per residue. The filter HMM depends on the model's composition (`p7_bg_SetFilter` per model), so it can't simply be cached per sequence. A 2-row/streaming SIMD version with a periodic log of the scale factor is ~10× cheaper. Keep float evaluation order in mind for byte-identical outputs |
| **S3: Viterbi filter** (`p7_ViterbiFilter`, 8×int16) | **2.3-3.1 G cells/s** | 9-15% of worker CPU | **Expected** (kernel microbench 2.8 G cells/s) | AVX2/AVX-512 would give 2-4× the lanes. ~1-3% of targets reach it |
| **S4: Forward parser** (`p7_ForwardParser` → `forward_engine`, 4×float) | **~0.5 G cells/s** (~6.4 cycles/cell) | 3-8% typical. **23% at M=1000**, 12-14% on hit-rich sets | **Expected** for SSE float | Share grows with M and hit density. Wider float SIMD helps directly. Only ~0.1-1% of targets |
| S5: Backward parser | 0.4-0.5 G cells/s | 1-2% typical, 7-11% hit-rich, 26% `--max` | Expected | Only for targets passing F3 |
| D1: domain decoding (`p7_DomainDecoding`) | — | < 0.3% | Expected | — |
| **D2: stochastic traceback clustering** (`region_trace_ensemble`: 200 × `p7_StochasticTrace`, `p7_Null2_ByTrace` per sampled domain, `p7_spensemble_Cluster` → `esl_cluster_SingleLinkage`/`link_spsamples`) | **14 µs/sample at M=205, 30 µs/sample at M=992. 2.8-6 ms per multi-domain region** | 5% (p100) → **12% (Pfam-1k × human `--cut_ga`) → 17-24% (M=1000, hitters × Swiss-Prot) → 41% (`--max`)** | **Inefficient / bottleneck on hit-rich and long-model workloads** | `nsamples = 200` is hard-coded (`p7_domaindef.c:110`). Per sample, `p7_Null2_ByTrace` is O(M·K) and the traceback walks with scalar `esl_vec_FNorm` + `esl_rnd_FChoose` per step, then `p7_trace_Reverse`. Clustering is O(n²) callbacks (n = samples × domains). Constraint: outputs depend on the RNG stream (reseeded per region), so byte-identical speedups must keep the exact sampling sequence: vectorize null2-by-trace, cheapen per-step choice, sort-based clustering. Fewer samples would change outputs |
| D3: envelope / region Forward + Backward (`p7_Forward`, `p7_Backward` full matrices) | 0.38-0.5 G cells/s | 1-4% typical, 7-10% at M=1000 or hit-rich | Expected | Re-runs F/B on each envelope (unihit) and on multi-domain regions (multihit) |
| D4: envelope posterior decoding | — | < 1% | Expected | — |
| D5: null2 by expectation | — | < 1% | Expected | — |
| D6: optimal accuracy DP + OA traceback | — | 0.3-3.5% | Expected | — |
| D7: `p7_alidisplay_Create` | — | < 0.3% | Expected | — |
| **O: merge / sort / threshold / text output** (main thread, serial after each query) | 60 MB text + 48 MB Stockholm + 21 MB domtbl for hitters × Swiss-Prot in ~0.7 s | ≈ 1% CPU, ~5% of wall in that output-heavy run | Expected (minor serial section) | Serial and outside the worker region. It adds to the per-query barrier stall. Becomes relevant only once the reader bottleneck is gone |
| Memory | 20 MB/thread base. Max seen 335 MB (M=1000 × Swiss-Prot, 16 threads). Titin target with 472 domains: 9 MB | — | **Not a concern** | Region/envelope DP matrices are sized to the region, not the full sequence |

---

## 4. Stage mix per workload (`--cpu 1`, % of all sampled CPU; main = reader/output thread)

| Workload (models × targets) | R read/parse | S1 MSV | S2 bias | S3 Vit | S4 Fwd | S5 Bck | D2 stotrace | D3-D7 | other |
|---|---|---|---|---|---|---|---|---|---|
| p100 × sprot_50k (100 Pfam models) | **30.9** | 36.3 | 6.1 | 11.0 | 5.7 | 0.7 | 5.1 | 1.1 | 3 |
| pfam_1k × human `--cut_ga --domtblout` (realistic annotation) | **26.2** | 31.4 | 7.7 | 9.6 | 6.3 | 1.7 | **12.1** | 2.6 | 2 |
| kern_M50 × Swiss-Prot | **54.1** | 26.9 | 5.6 | 3.6 | 1.4 | 0.2 | 0 | 0 | 7.3 (per-seq) |
| kern_M200 × Swiss-Prot | 20.4 | 26.8 | **11.6** | 9.3 | 7.2 | 1.1 | **18.1** | 3.9 | 1.5 |
| kern_M400 × Swiss-Prot | 25.1 | 53.8 | 5.0 | 9.6 | 4.2 | 0.1 | 0 | 0 | 1.6 |
| kern_M1000 × Swiss-Prot | 4.5 | 24.1 | 2.8 | 10.8 | **23.3** | 6.9 | **17.0** | 9.9 | 0.8 |
| len_mix (40 models) × sprot_10k | 14.5 | 41.0 | 4.4 | 12.8 | 12.7 | 1.4 | 2.7 | 4.7 | 2.5 |
| len_xlong (10 models, M > 1000) × sprot_10k | 6.2 | 40.8 | 3.1 | 14.4 | 16.6 | 2.2 | 4.2 | 7.4 | 3 |
| hitters × sprot_10k `--noali` | 14.4 | 12.6 | 5.0 | 5.5 | 11.7 | 10.9 | **21.7** | 13.3 | 3 |
| hitters × Swiss-Prot, full output + `-A` (`--cpu 16`) | 15.8 | 13.1 | 3.9 | 5.2 | 11.7 | 9.3 | **23.6** | 13.0 | 3 (O 1.2) |
| pfam_10 × random_100k (null DB) | **34.0** | 43.8 | 4.6 | 11.1 | 3.1 | 0 | 0 | 0 | 2.3 |
| len_med × sprot_long (≥ 3,000 aa) | **30.6** | 34.2 | **13.5** | 12.2 | 5.3 | 0 | 0 | 0 | 3 |
| len_mix × sprot_short (≤ 60 aa) | **24.8** (+12 unresolved main) | 36.5 | 3.6 | 3.8 | 0.7 | 0 | 2.7 | 0 | 10 (per-seq) |
| len_med × sprot_1k `--max` | 0.5 | 0.3 | 0 | 0 | 19.9 | 26.1 | **41.2** | 8.8 | 2 |
| kern_M200 × sprot_50k.dat (EMBL) | **50.8** | 14.7 | 7.4 | 5.1 | 4.3 | 0.5 | 12.8 | 1.8 | 2 |
| kern_M200 × sprot_50k.fa.gz | 15.8 + **33.8 gzip child** | 15.8 | 7.8 | 5.3 | 2.9 | 0.6 | 12.5 | 3.3 | 2 |

Filter pass rates are at their nominal values in all workloads (MSV 2-11%, bias 1-4%, Vit 0.1-1%, Fwd ≤ 0.6%), so none of the filters is leaking.

---

## 5. Thread scaling (wall seconds, build-dev, no perf)

| workload | cpu0 | 1 | 2 | 4 | 8 | 12 | 16 | max CPU util | limited by |
|---|---|---|---|---|---|---|---|---|---|
| p100 × sprot_50k | 23.2 | 17.6 | 10.6 | 8.83 | **8.74** | 9.41 | 9.85 | 2.9× | reader (98% busy at cpu16) |
| pfam_1k × human `--cut_ga` | 140.2 | 109.4 | 66.0 | 52.5 | **52.1** | 55.9 | 56.7 | 3.0× | reader (~42 s of parsing) + per-query serial |
| kern_M200 × Swiss-Prot | 4.04 | 3.40 | 1.71 | 1.05 | **1.03** | 1.03 | 1.04 | 4.4× | reader (one 1.0 s parse) |
| len_mix × sprot_10k | — | — | 1.80 | 1.29 | **1.23** | 1.27 | 1.29 | 3.6× | reader |
| hitters × sprot_50k | 4.59 | 4.04 | 2.20 | 1.40 | 1.24 | **1.23** | 1.26 | 4.5× | reader |
| kern_M1000 × Swiss-Prot | 19.6 | 18.8 | 9.47 | 4.85 | 2.58 | 2.21 | **2.01** | 14.3× | compute: linear to 8 cores, +28% from HT |

- Compute alone scales well, linearly up to the physical core count (M1000 case). Hyperthreads add about 25%, as expected for SIMD-saturated kernels with IPC ≈ 3.7.
- When the reader is the limit, adding threads beyond about 8 *hurts*. Per-query thread creation, HT contention with the reader, and more blocks in flight all grow with N.
- `--cpu 1` already beats `--cpu 0` by 1.2-1.4×, because parsing overlaps one worker.

## 6. Many queries × small DB (per-query fixed cost)

| run | cpu0 | cpu1 | cpu4 | cpu16 |
|---|---|---|---|---|
| Pfam-A (21,979 models, text) × 3 seqs | 20.2 s | 22.0 | 25.0 | **39.0** (34.7 user + 44.7 sys) |
| Pfam-A (hmmpress'ed `.h3m`) × 3 seqs | 4.3 s | — | — | **22.5** |
| pfam_1k × sprot_1k | 5.09 | 5.52 | — | **6.25** |

In this regime more threads always make it slower.

## 7. Compiler flags

A `-O3 -march=native` build (AVX-512 available) ran within ±3% of the default build on M1000 × Swiss-Prot, p100 × 50k and hitters × 50k. It also **changed outputs in the last digit** of some E-values (11 lines in p100 × 50k), presumably from FMA contraction. So flags give no speedup, and they break the byte-identical golden regime. Real SIMD gains need hand-written 256/512-bit kernels.

## 8. Prioritized acceleration directions (analysis only)

Amdahl framing for the typical Pfam × proteome case (pfam_1k × human), from the CPU split: reader 26%, filters S1-S5 57%, domain stage 15%, other 2%.

1. **Remove the per-query re-parse and the serial reader.** Keep the digitized DB resident and stream it to workers for every query, or swap the loop order and batch many models per DB pass. Expected effect: wall drops from ~52 s to roughly CPU/cores ≈ 105/8 ≈ 13-15 s on this box (3.5-4×). Unlocks scaling to more cores. Prerequisite for everything below.
2. **A persistent thread pool and per-thread state reused across queries**, plus pressed or pre-converted models. Removes the ~1 ms/query fork/join and the text HMM parse. 4.7× on many-model × small-DB from pressing alone. Also needed for the multi-query `.gz`/stdin limitation.
3. **Wider SIMD kernels** (AVX2/AVX-512) for SSV/MSV, Viterbi filter and Forward/Backward. These are already at the 128-bit peak, so lane count is the lever: up to 2-4× on ~60-70% of worker CPU.
4. **Stochastic traceback clustering (D2).** Optimize under an RNG-stream-preserving constraint. 12-24% of CPU on realistic/hit-rich sets.
5. **Bias filter.** Remove per-call malloc and per-residue `log()`, and vectorize. 4-12% of CPU, potentially ~10× cheaper.
6. **Input.** In-process or parallel decompression, larger reads, no per-read `lseek`, faster header scan. Matters only while #1 isn't done, or for one-pass single-query runs.

Correctness constraint for all of these: `dev/check.sh` requires byte-identical outputs. #3, #4 and #5 change floating-point evaluation order and RNG use, so each needs either an exact-order implementation or an explicit, reviewed golden update.

## 9. Method and caveats

- `build-prof/` uses the same `-O3` plus `-g -fno-omit-frame-pointer`. Its outputs are byte-identical to build-dev and its speed is within noise.
- Stage attribution: `perf record -g` at 997 Hz. Each sample is assigned to the innermost known pipeline function on its call chain (`stages.py`). `forward_engine`/`backward_engine` are tail-called, so they're split into filter vs domain stage by whether `p7_domaindef_*` is on the stack. libc leaves are attributed through their callers. Unresolved main-thread frames ("X other/libc: [unknown]") are ≤ 2% except in `shortseqs` (12%, almost certainly parse/libc).
- `perf record` adds a ~1.2 s wall-time floor to small runs. All wall times in §5-§7 come from runs without perf.
- Throughput = stage seconds (perf) ÷ cells (bpftrace counts of ΣL per kernel call × M). Viterbi/Forward/Backward cells use the exact L of each call.
- This is a KVM guest with hyperthreading. Absolute numbers will shift on other hardware, but the ratios (reader vs compute, per-stage cost per cell) carry over.

## 10. Reproduce

```
dev/profiling/build-prof.sh                                   # profiling build in build-prof/
D=~/.cache/hmmer-bench/derived-b36ad365c6d1
dev/profiling/profile_stage.sh human_pfam1k 1 --cut_ga --domtblout h.dom $D/pfam_1k.hmm $D/human.fa
dev/profiling/scale.sh p100_50k 1 --noali $D/pfam_100.hmm $D/sprot_50k.fa
dev/profiling/cells.sh --noali $D/kern_M1000.hmm $D/sprot.fa  # sudo bpftrace: calls + ΣL per stage
```

Raw artifacts (perf.data, `.stages`, outputs) go to `bench-results/profiling/`, which is gitignored.

## 11. Open follow-up: full Pfam × millions of proteins (not yet measured)

Extrapolated linearly from pfam_1k × human (`--cut_ga`; 1/22 of Pfam-A, human proteins average ~560 aa), one hmmsearch process with full Pfam-A costs, per target protein:

- ≈ 0.12 core-seconds of search work (parallel)
- ≈ 0.045 s of serial re-parsing (the reader floor)

For 1 M proteins that is ~34 core-hours of search, but a ≥ ~12.5 h wall-time floor from parsing.

**Workaround with no code change:** split the targets into P files and run P processes, each with `--cpu 2-3`. Give every process the same `-Z` (total number of targets) so E-values are comparable, and read models from an hmmpress'ed Pfam file. This divides the parsing floor by P.

Still to measure; the existing tooling is enough, apart from small additions to `stages.py`:

1. All of Pfam-A with `--cut_ga` against 20k and 100k samples of the real target set. Confirm cost grows linearly, then extrapolate.
2. hmmscan on the same sample. It re-reads the pressed Pfam MSV data (~0.4 GB) for every target sequence and is known to use few threads. Add `p7_oprofile_ReadMSV`/`ReadRest` and the hmmscan loop to `stages.py`.
3. Memory and the serial merge/sort/output time at 10k → 100k → 1M targets, for hit-rich families (hit lists are held in memory per query).
4. The stage mix on the real target data. Short fragments carry ~10% per-sequence overhead.
