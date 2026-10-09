# Agent prompt: an opt-in GPU path for `hmmsearch` (`--gpu`)

*How to use: give everything below the line to the agent as its task. It runs on this machine, in `/home/asokolov/Projects/hmmer`. Written 2026-10-09, with `dev` at 5935bd9.*

---

## 0. Your task

You are adding an **opt-in GPU path to HMMER's `hmmsearch`**: a run-time flag `--gpu`, built only when CUDA is enabled at configure time. With `--gpu`, the filter stages of the search pipeline run on an NVIDIA GPU. Everything else stays on the CPU.

The output must be **byte-identical** to the CPU path's. The only exception is the option's own echo line in the header.

This repository is the user's fork of HMMER 3.4. Over nine days it was made **7-16× faster on CPU**, with byte-identical outputs throughout (Phases 0-7). That work produced:

- a strict test suite;
- an A/B performance protocol;
- a long list of measured lessons.

You inherit all three. **Most of the risk in this task is in measuring and verifying, not in writing kernels.** Read §3 before writing any code.

The user's deployment targets are the **NVIDIA H100** (this machine) and the **RTX PRO 6000 Blackwell**. You measure on the H100 and project to the RTX PRO 6000 (§9.5).

Work in small, verified steps. Each step must be:

- byte-identical (proved by the gates in §7);
- measured with the protocol in §8;
- recorded as described in §9.

Stop at the decision gates in §6 and ask the user.

---

## 1. Ground rules (non-negotiable)

1. **Byte-identical output.**
   - `hmmsearch --gpu ...` must produce the same stdout, stderr, `-o`, `--tblout`, `--domtblout`, `--pfamtblout` and `-A` files as `hmmsearch ...`, after the existing normalization (dates and timings masked).
   - The "Internal pipeline statistics summary" counts are part of the output and must match. These are `n_past_msv`, `n_past_bias`, `n_past_vit`, `n_past_fwd` and the residue counts.
   - The only accepted differences:
     - the `--gpu` option's own header echo line and its token in a table's `# Option settings:` line;
     - whitespace-only column-width differences in threaded runs **without `-Z`**. This is a known upstream flake that `ab.sh` classifies as `ws`; see `dev/profiling/LEARNINGS.md` §Determinism.
2. **The default path stays unchanged.**
   - Without `--enable-cuda` at configure time, the build must not change behaviour or speed.
   - In a CUDA build run without `--gpu`, the CPU code path must be literally unchanged.
   - Prove both with instruction counts (`instructions:u` 1.0000 on the `cpu` and `reader` sets) against a snapshot of `dev` (§8.1).
3. **Opt-in at two levels.**
   - `./configure --enable-cuda` (default **off**).
   - `hmmsearch --gpu` (default off).
   - A binary built without CUDA must fail `--gpu` with a clear error.
   - No CUDA device, or an unsupported case (v1 is the amino-acid alphabet only), must also fail with a clear error.
   - **Never fall back to the CPU silently.**
4. **No non-byte-identical shortcuts.** Some stages are float code where exactness on a GPU is expensive (the bias filter, the Forward/Backward parsers). Keep them on the CPU, or reproduce their exact operation order. Relaxing exactness, even behind another flag, is the user's decision. Escalate (§1.6).
5. **Git.**
   - Branch from `dev`: one branch per phase (`perf/gpu-g0`, `perf/gpu-g1`, ...) with one commit per step, or one branch `perf/gpu` with one commit per step.
   - No global git identity is configured on this machine. Commit with `git -c user.name=asokolov -c user.email=sokolov.ars.a@nebius.com commit ...`.
   - Messages follow the existing style (see `git log -5`): a `component: what (step Gx.y)` subject, then a body with what changed, why it is byte-identical, and the measured numbers. End with the trailer the user's tooling expects (`Co-Authored-By: ...`, as in recent commits).
   - **Never push.** This shell has no GitHub credentials, and the user does all pushes.
   - Fast-forward `dev` only after a phase passes its milestone (§6).
6. **Stop and ask the user** before any of these:
   - (a) relaxing byte-identity anywhere;
   - (b) making CUDA or `--gpu` a default;
   - (c) adding a dependency beyond the CUDA toolkit (CUB and Thrust ship with it and are fine);
   - (d) any Easel change. Easel is pinned through `dev/easel.ref` to the user's fork, `hazirliver/easel` branch `hmmer-perf/phase5-fasta`; the user's clone is `~/Projects/easel` and the user pushes it;
   - (e) changing machine or GPU settings: clock locking, persistence mode, MIG, anything needing `sudo` beyond read-only probes;
   - (f) going past a decision gate in §6;
   - (g) choosing the chart type or metric for a new results section (§9.4).
   - Also: **don't re-propose** the query-parallel output-window fix for small databases. The user declined it on 2026-10-09.
7. **Report faithfully.**
   - Measured on the H100 means measured. RTX PRO 6000 numbers are projections and are labelled as such.
   - A skipped check is reported as skipped.
   - A slower workload is reported as `SLOWER`, not hidden.

---

## 2. Machine and environment

| | |
|---|---|
| CPU | 16-vCPU KVM guest: Intel Xeon Platinum 8468 (Sapphire Rapids), **8 cores × 2 hyperthreads**, ~3.2 GHz under load, AVX-512 |
| RAM / disk | 196 GB / 1.3 TB free on `/` |
| GPU | **1× NVIDIA H100 80GB HBM3** (SXM, 700 W cap), compute capability **9.0**, 132 SMs, max SM clock 1980 MHz, driver 580.173.02 |
| CUDA | Toolkit **13.0** in `/usr/local/cuda-13.0` (`nvcc` on PATH); `compute-sanitizer`, `nsys`, `ncu`, `cuobjdump`, `nvdisasm` available |
| Host tools | gcc 13.3, `perf`, `bpftrace` (needs sudo), python3 **without matplotlib** (the plot scripts write SVG by hand). No `gdb`, no `valgrind` |
| Benchmark data | `~/.cache/hmmer-bench/derived-b36ad365c6d1/`: Swiss-Prot `sprot.fa` (570,420 seqs, 206 M residues) and subsets `sprot_{1k,10k,50k}`, `sprot_long`, `sprot_short`; `human.fa` (20,428 seqs, 11.4 M res); `random_100k.fa`; Pfam 37.0 `Pfam-A.hmm` (21,979 models, ΣM 3.59 M) and subsets `pfam_{10,100,1k}`, `pfam_hitters`, `pfam_len_*`, `kern_M{50,200,400,1000}`. Pressed Pfam-A: `~/.cache/hmmer-bench/ab-fixtures/`. UniProt annotation sets: `~/.cache/hmmer-bench/uniprot-ccf4a8601189/` |

Deployment targets (the RTX PRO 6000 isn't on this machine; check its specs before relying on them):

| GPU | Arch / CC | SMs | Memory | Host link | Notes |
|---|---|---|---|---|---|
| H100 (here) | Hopper, sm_90 | 132 | 80 GB HBM3, 3.35 TB/s | — | Hardware DPX; TMA, clusters, wgmma (Hopper-only) |
| RTX PRO 6000 Blackwell (Server) | GB202, sm_120 | 188 | 96 GB GDDR7, ~1.8 TB/s | PCIe 5 x16 | Unified FP32/INT32 cores |

Consequences:

- Build a fatbin for `sm_90` and `sm_120`. Every build must compile both.
- **Don't put Hopper-only features in the main path** (TMA, thread-block clusters, distributed shared memory, wgmma). The RTX PRO 6000 doesn't have them.
- DPX intrinsics (`__vimax3_s16x2` and friends) compile everywhere but are hardware-accelerated on Hopper. Whether GB202 accelerates them is unverified, so an H100 number from a DPX-heavy kernel may overstate the RTX PRO 6000. Measure a DPX-free variant as well (§9.5).
- Design for **80 GB** of device memory (the H100, the smaller of the two).

Machine quirks:

- The VM's speed drifts. On 2026-10-09 it ran ~12% slower for hours. Compare only interleaved runs.
- Bash tool commands in the foreground stop after 10 minutes, and background tasks are capped at 2 h. Start long runs detached (`setsid nohup script.sh > log 2>&1 &`, with the script appending its exit status) and wait for the completion notice; don't poll.
- Before every timed GPU session, check that the GPU is idle (`nvidia-smi`). It may be shared.

---

## 3. Read before writing code (in this order)

1. `dev/README.md`: the test suite, the golden regression and its normalization.
2. `dev/profiling/PLAN.md`:
   - §2: how byte-identical testing works, and the invariants that make reordering safe (RNG reseeded per region; hits sorted by key, name and position);
   - §3: the per-step protocol and acceptance rules;
   - §6: the user's decisions.
3. `dev/profiling/LEARNINGS.md`, **all of it**. Most relevant here:
   - Measuring 1-3, 6, 9, 15, 16, 19, 20;
   - Byte-identical FP 1-3;
   - Phase 6 items 1, 5-9, 11 (exact integer DP at other widths: the layout-independence argument, boundary scans, mutants, "one binary, two paths");
   - Phase 7 items 1-3, 6-8, 10 (barriers, batch memory, error paths, what deferred stages read, read granularity);
   - Determinism;
   - Design patterns 1, 2, 7;
   - Verification habits;
   - Tooling gotchas.
4. `dev/perf/README.md`: `snap.sh`, `ab.sh` (`--b-args=`, `--a-env/--b-env`, `--check-only`, `--file`), `series.py`, `plot_series.py`, `plot_annot.py`, `errpaths.sh`, crosscheck builds, measured noise.
5. `dev/profiling/README.md` §2-§4 (pipeline anatomy, cost per stage, stage mix per workload) and `dev/profiling/REPROFILE.md` §1-§3 and §6 (what bounds threaded runs; scaling out).
6. `dev/profiling/PROGRESS.md` for the format of step records, and `dev/profiling/RESULTS.md` for the style and method of the final report.
7. Code:
   - `src/p7_pipeline.c`: `p7_Pipeline()` = `p7_pli_Filter()` → `p7_ForwardParser()` → `p7_pli_ForwardPass()` → `p7_BackwardParser()` → `p7_pli_Domains()`. Read the contract in `p7_pli_Filter()`'s header comment. **This split is your seam.**
   - `src/hmmsearch.c`:
     - `thread_loop`, `pipeline_thread`, `worker_Search`;
     - `WORKER_QUERY` slots, with `MAXQ` = 16 queries per pass and two slot sets for overlapping batches;
     - `pair_feed`/`pair_flush`, which pair targets for the AVX2 two-target parsers;
     - `split_is_free`, `dbcache_*`, `output_batch`, `worker_Trim`.
   - `src/impl_sse/`:
     - `ssvfilter*.c`, `msvfilter*.c`, `vitfilter*.c` (SSE / AVX2 / AVX-512 kernels, chosen in `simd.c`, `HMMER_SIMD`);
     - `fwdback.c`, `fwdback_pair.c`;
     - `p7_oprofile.c` (striped profile layouts);
     - `impl_sse.h` (`impl_Init()` sets FTZ and DAZ at lines ~405-420; GPU float code must match).
   - `src/p7_bg.c`: the null model and the bias filter (`p7_bg_FilterScore`, a streaming 2-state forward with a per-residue `log`).

---

## 4. What is already known (starting estimates; verify them)

### 4.1 Current CPU state (`dev` @ 5935bd9, 16 threads, `RESULTS.md`)

- Realistic annotation: Pfam 1k × human proteome, `--cut_ga` (workload `X5`). **7.9 s**; upstream 56 s.
- All of Pfam-A × 5 sequences: 1.5 s pressed, 7.3 s from text.
- Pfam-A annotation of 6 × 1,000 UniProt proteins: **71 seq/s** pressed; upstream 5.2.
- Already done on CPU, so don't redo it:
  - AVX2 SSV/MSV/Viterbi filter, with AVX-512 for M > 448;
  - query batches (16 per pass over the targets, overlapped);
  - an exact streaming bias filter;
  - Forward/Backward on two targets per AVX2 register;
  - a fast FASTA reader (~1 G residues/s, single thread);
  - opt-in `--dbcache`.

### 4.2 Where CPU time goes (`--cpu 1`, `REPROFILE.md` §2, parsing excluded; Phase 7 amortizes it over 16 queries)

| Stage | Realistic (X5) | Hit-rich | Fits a GPU | Exact on a GPU? |
|---|---:|---:|---|---|
| SSV/MSV filter (int8, saturating max-plus) | ~41% | ~15% | Very well | **Yes**: order-independent integer DP, if the saturation, bias and overflow semantics are reproduced |
| Viterbi filter (int16) | ~12% | ~5% | Well; published GPU results are weakest here | **Yes**, same argument (see LEARNINGS Phase 6 item 1 on the DD bound) |
| Bias filter (float, scalar per target, per-residue `log`) | ~11% | ~2% | One target per thread | Only with an op-for-op copy, no FMA, FTZ matched, and a `log` bit-identical to glibc's |
| Forward/Backward parsers (float, striped SSE) | ~15% | ~32% | Well | Only by emulating the SSE striped operation order exactly (D→D passes, horizontal sums, rescaling) |
| Domain stage + glue: stochastic traceback (200 samples per region, serial RNG stream), envelope F/B, OA, null2, alignment display | **~21%** | **~45%** | Poorly | Practically no. **Stays on the CPU** |

Filter pass rates: MSV 2-11%, bias 1-4%, Viterbi 0.1-1%, Forward ≤ 0.6% of targets.

### 4.3 Amdahl ceilings (if the GPU did SSV, MSV, bias, Viterbi and the parsers for free)

- Realistic: ~4.8×.
- Hit-rich: ~2.2×.
- `--max`: ~1.9×.
- Few hits (e.g. against `random_100k`): 20× or more; then the target feed is the bound.

With only SSV/MSV + Viterbi on the GPU, the realistic ceiling is ~2.5-3×.

### 4.4 Scale

- X5's SSV work is ΣM 161,888 × 11.4 M residues = **1.85e12 cells**. On the CPU that is ~26 core-seconds at **71 G cells/s per core** (in situ, AVX2). On a GPU at ~5-15 TCUPS it is 0.1-0.4 s.
- Expected end to end: X5 7.9 s → ~2 s, with the **GPU idle most of the time**. One such GPU has the filter throughput of roughly 150-200 of these cores (an estimate from published GPU filter rates), so it balances ~40-50 host cores on annotation work.

### 4.5 New bounds a GPU exposes

- **Feed rate.** A GPU consumes residues at (GPU CUPS ÷ ΣM of the query batch). At 12 TCUPS with 16 Pfam models (ΣM ~2,600) that is ~4.6 G res/s, while the reader parses ~1 G/s. The fixes are larger query batches per pass (more memory: a batch's hit lists are held until printed), a device-resident or pinned `--dbcache`, or a pre-digitized database (`REPROFILE.md` §6). A single-model search against a big database needs the database resident on the device.
- **Domain-stage tails.** On long proteins (annotation group human L4), a few (query, target) domain-stage items of 0.1-0.7 s already set the pace. A GPU doesn't touch them.
- **CPU per-pair glue.** Turning a raw filter score into a P-value is ~tens of ns per (query, target) on the host.

### 4.6 Published references

| Work | Result |
|---|---|
| CUDAMPF (K40, 2016) | SSV 440, MSV 277, **Viterbi only 14.3** GCUPS |
| CUDASW++4.0 (Smith-Waterman with affine gaps; similar to the Viterbi filter) | **H100 5.71 TCUPS** |
| MMseqs2-GPU (gapless prefilter, similar to SSV) | ~100 TCUPS on 8 GPUs. One GPU was 177× faster than JackHMMER on 128 cores, but **~6×** vs CPU MMseqs2 on 2×64 cores. That like-for-like factor is the realistic expectation here |

---

## 5. Architecture guidance

These are starting points. Change them when measurements say so, and record why.

1. **Where to cut.**
   - v1: the GPU computes the **raw integer filter scores** for each (query, target): the MSV/SSV byte score (and overflow), later the Viterbi filter's int16 score.
   - The host converts raw scores to P-values and pass/fail **with the existing float code** (`p7_pli_Filter()`'s arithmetic, `esl_gumbel_surv`). That makes it exact by construction.
   - Survivors continue on the CPU workers through bias, VF (until it moves), parsers and domains, in each target's normal call order. The pipeline counters are updated exactly as `p7_pli_Filter()` updates them.
   - Any host shortcut must be proved exact by evaluating the same code, and covered by the crosscheck. An example is a per-(query, length) table of passing raw scores. Count first how often it would fire (LEARNINGS Measuring 10).
2. **What each later stage reads.**
   - Before a target's later stages run, restore everything they read (LEARNINGS Phase 7 item 6): `om`/`bg` length configuration (`p7_bg_SetLength`, `p7_oprofile_ReconfigLength`), `nullsc`, `filtersc`, and the pipeline's running `Z` when there is no `-Z`.
   - The parsers' outputs consumed by the domain stage are only `xmx`, `L`, `has_own_scales` and `totscale`.
3. **Batching.**
   - Reuse the Phase 7 scheduler: batches of queries per pass, overlapped, printed in query order.
   - The GPU path will want far more than 16 queries per pass, or the whole database resident. Size this by memory: device memory and the host hit lists (LEARNINGS Phase 7 item 2).
   - Decide modes once per run from fixed facts such as database size (LEARNINGS Design 7, Measuring 19).
   - Change the payload, not the work-queue protocol (Design 2).
   - Don't add a per-query barrier (Phase 7 item 1).
4. **Kernels.**
   - Write the layout-independence argument in each kernel's header, as for the AVX2 kernels: padding nodes, post-overflow work, per-layout shortcuts.
   - Handle short models (M ≤ 64) and long ones (M > 1000, up to 2,000+). Model length decides which kernel wins (Phase 6 item 11).
   - Bin targets by length to limit divergence. Results are keyed by target index, so binning can't change outputs.
   - Check the emitted SASS (`cuobjdump -sass`). Saturating SIMD-in-word intrinsics (`__vaddss4`, `__vmaxs4`, ...) are emulated with several instructions on recent architectures.
5. **Float device code, if you get there (bias, parsers).**
   - Compile with `-fmad=false`, FTZ matched to `impl_Init()` (`-ftz=true` only if the CPU path really flushes in that code), `-prec-div=true -prec-sqrt=true`, no fast-math, no `__expf`/`__logf` intrinsics.
   - A `log` bit-identical to glibc's exists as a starting point: `dev/profiling/reprofile/fastlog.c`, checked on all 2^32 float inputs. It falls back to libm on ~1.5e-5 of inputs, and on the device those can come from a table.
   - The parsers must reproduce the SSE striped evaluation tree. `fwdback_pair.c` and its utest show how this was done for two targets per AVX2 register.
6. **Host side.**
   - CUDA context creation costs ~0.1-0.5 s, which is visible on sub-second runs (X2, X4). Overlap it with reading the HMMs and targets, and report it.
   - Use pinned buffers with double buffering, and `cudaHostRegister` for the `--dbcache` arena.
   - One process-wide context; either a dedicated GPU feeder thread or streams per worker. Measure which.
   - Device errors are asynchronous. Check after syncs, and run debug builds under `compute-sanitizer`.
7. **Build.**
   - `configure.ac`: `AC_ARG_ENABLE(cuda)`, default no; find `nvcc`; `CUDA_ARCHS` default `90;120`.
   - Put the code in a new `src/impl_cuda/` (or similar) with its own `Makefile.in`.
   - Link `cudart_static`, so a copied binary stays a complete snapshot (`snap.sh` relies on that). Only `libcuda.so` from the driver is loaded at run time.
   - HMMER's Makefiles don't regenerate themselves. New utests go into `UTESTS` in `src/Makefile.in` and `testsuite/testsuite.sqc`, and existing build dirs need `./config.status`.
   - Prove every snapshot contains the new code: `nm` for a new symbol, plus a negative control (LEARNINGS Phase 6 item 8).
   - Teach `snap.sh` a way to pass configure flags (`--enable-cuda`), and confirm `snap.sh --rev` still builds the Easel commit each revision pins.

---

## 6. Work plan and decision gates

Each step: implement → quality gates (§7) → performance A/B (§8) → record (§9) → commit. A neutral step is kept only if it's a declared prerequisite.

**G0. Ceiling and baselines (no GPU code in `src/`).**

- Snapshot the reference: `dev/perf/snap.sh --rev 5935bd9 --perf` and a default `-O3` snapshot.
- Measure the real Amdahl ceiling instead of inferring it from profile shares (LEARNINGS Measuring 15, Phase 6 item 10). Use a throwaway patch, kept under `dev/profiling/gpu/` and never in `src/`:
  - run 1 records every (query, target) filter outcome: raw scores and pass flags;
  - run 2 replays them and skips SSV/MSV/VF (variant A), then also the bias filter and parsers (variant B);
  - the outputs must stay identical, and the replay's wall time at 16 threads is the floor for a perfect GPU.
- Run it on X5, X1, X2, W1 (hit-rich), W5 (`--max`), X4 pressed, and the annotation groups `A_human_L1`, `A_human_L4`, `A_bacteria_L2`.
- Measure per-pass idle time with timestamps (LEARNINGS Measuring 16; `dev/profiling/reprofile/barrier-stats.patch` no longer applies after 7A, so adapt the idea).
- **Gate:** report the ceilings to the user. If the realistic (X5) and annotation ceilings are below 2× with the GPU free, stop and ask whether to continue.

**G1. Exact SSV/MSV kernel, standalone.**

- The kernel plus a permanent bitwise utest against `p7_SSVFilter`/`p7_MSVFilter`:
  - models from `p7_hmm_SampleUniform` (Pfam-like) and `p7_hmm_Sample`;
  - target lengths 0, 1 and 2; all residue codes, including degenerate ones;
  - scans that hit exact boundaries: consensus prefixes so the score climbs through the overflow threshold, and consensus with 1-8 nodes deleted at every position;
  - assert that every outcome (ok, overflow) is reached;
  - the utest skips cleanly when there is no device.
- Mutation tests that the utest must catch, for example an off-by-one in the overflow test, or a wrong bias.
- Real data: every (query, target) raw score of `pfam_1k × human` and `kern_M* × sprot` equals the CPU's.
- Microbenchmark: GCUPS by M (40, 50, 100, 200, 400, 1000, 2000) and by length mix, on the H100, with and without DPX.
- **Gate:** report the H100 GCUPS and the RTX PRO 6000 projection. If the H100 measures below ~2 TCUPS (~30 Sapphire Rapids cores), stop and ask.

**G2. `--gpu` integration, with SSV/MSV on the GPU and the rest on the CPU.**

- The configure flag, the option (help text, header echo line, error messages), the host pipeline, and the survivor hand-off.
- All of §7. Performance per §8.2 on every set.
- Re-profile the host side (per-thread `perf script -F tid,period` sums, `nsys` timeline) and name the new bound.

**G3. Viterbi filter on the GPU** (int16, exact). Same verification as G1. It runs on only 2-11% of pairs, so measure whether launch and transfer overheads let it beat the CPU kernels, per model length.

**G4. Feed and scheduling.**

- Device-resident targets where they fit.
- GPU-sized query batches.
- Pinned double-buffering.
- Several hmmsearch processes sharing the GPU (query sharding was the best CPU configuration, `REPROFILE.md` §6). Measure 1×16, 2×8 and 4×4.
- Goal: keep host workers busy on the domain stage while the GPU filters the next work. Report GPU busy % and host cores busy.

**G5 (optional, gated by the residual measured after G4).** The bias filter on the GPU (exact `log`), then the Forward/Backward parsers on the GPU (exact SSE-order emulation). Each needs the full kernel evidence (§7.3). If exactness proves infeasible or too slow, **escalate with numbers**; don't relax it.

**Milestone after each phase.**

- All workload sets, both files (`workloads.tsv`, `annot.tsv`), default `-O3` builds.
- `dev/check.sh --bench M` on the default build.
- `PROGRESS.md` milestone rows.
- Fast-forward `dev`.

**End: report** (§9).

---

## 7. Quality gates (zero tolerance)

### 7.1 Every step

1. `dev/check.sh` on the **default build** (no CUDA): ALL PASSED (build, `make check`, all ~97 golden cases).
2. `HMMER_BUILD_DIR=$PWD/build-cuda CONFIGURE_FLAGS=--enable-cuda dev/check.sh`: ALL PASSED. Your new utests run in `make check`.
3. `HMMER_BUILD_DIR=$PWD/build-cuda dev/check.sh --force-opt --gpu`: every golden hmmsearch case with `--gpu` forced equals its golden, with the option's echo lines masked.
   - Check that `check.sh`'s mask and `ab.sh`'s `--b-args` mask actually cover the `--gpu` echo line; extend them if not.
   - Cases that legitimately differ go into `dev/regress/force-opt.expected-diff` with a reason. Expect few or none, since DNA and other unsupported cases error out.
4. Real data, same binary: `dev/perf/ab.sh <cuda-bin> <cuda-bin> --b-args=--gpu --set all --check-only`, plus `--file annot.tsv --set annot,annot_txt --check-only` at milestones. It must print `outputs identical`.
5. Default path unchanged: `ab.sh` A = the `5935bd9` `--perf` snapshot, B = the new tree's `--perf` build without CUDA, `--set cpu,reader`. Instructions 1.0000, cycles `~`. Repeat with B = the CUDA build run without `--gpu`.

### 7.2 Threading or scheduling changes

- `dev/perf/errpaths.sh A B` with B running `--gpu`. Extend `errpaths.sh` to take per-side arguments, or use a wrapper. The cases cover:
  - GA missing on the first and the 23rd model;
  - a truncated HMM file;
  - mixed alphabets;
  - a target parse error in record 500;
  - batching with and without `-Z`, `--cut_ga`, `--dbcache` and `-A`.
- ThreadSanitizer: configure, build and run under `setarch -R`. Check first that it flags a toy race.
- ASan with leak detection: confirm libasan is linked with `ldd`.
- `compute-sanitizer` with `memcheck`, `racecheck`, `synccheck` and `initcheck` on the utests and a small search.
- A stress loop of 120+ runs at `--cpu` 0, 1, 2, 4, 8, 16 and 24 with `--gpu` under `timeout`. Every output must equal `--cpu 0` without `--gpu`.
- The `ws` frequency: 10 threaded runs without `-Z` (e.g. X1 `defZ`), before and after (LEARNINGS Determinism 3).

### 7.3 Every kernel (exact-rewrite evidence, LEARNINGS "Byte-identical rewrites" item 3)

1. A permanent bitwise utest against the **original CPU routine** itself, not a wrapper around the new code.
2. Mutation tests showing the utest catches near-misses. Write down which mutants survive and why (Phase 7 item 9).
3. A `-Dp7_CROSSCHECK` build (`HMMER_BUILD_DIR=$PWD/build-xchk-cuda CFLAGS="-O3 -Dp7_CROSSCHECK" CONFIGURE_FLAGS=--enable-cuda`). It reruns the CPU filter for every pair and `abort()`s on any difference. Run it over the goldens with `--force-opt --gpu`, plus `ab.sh ... --check-only --set all` and the annotation sets.
   - Configure crosscheck builds through `CFLAGS`. Never edit generated Makefiles.
   - Show that the crosscheck fires on a mutant.

### 7.4 Negative controls for every new checker

- `ab.sh` with a B that adds `--seed 7` must report DIFF.
- `--force-opt --no-such-opt` must fail.
- A deliberately broken kernel must fail its utest and abort the crosscheck.
- `hmmsearch --gpu` from a non-CUDA build must fail with the intended message.

---

## 8. Performance protocol

### 8.1 Rules carried over (LEARNINGS "Measuring")

- **Interleaved A/B in one session only** (`ab.sh`, `series.py`). Never compare against a stored baseline or another day's file.
- **Choose the metric that resolves the effect.**
  - CPU-path neutrality: `instructions:u` (0.00% noise) and `cycles:u` at `--cpu 0`, from `--perf` (`-falign-functions=64`) builds on both sides.
  - Everything with `--gpu`: wall time.
- **Never run threaded or GPU runs under `perf stat` when judging wall time.** It added +18% to +41% on this VM. `ab.sh` already counts only the `cpu`/`reader` sets; don't pass `--perf-all` for wall verdicts.
- **"One binary, two paths"** is the main GPU comparison. The same CUDA binary is A (no option) and B (`--b-args=--gpu`), so there is no code-layout confound.
- **A verdict counts only beyond noise.** It must exceed the floor (wall 2%, cycles 0.5%) and both sides' spread. Measured A/A spreads: cycles ±0.2%; threaded wall up to ±4.7% (use `--reps 9` for small effects). Runs under 0.5 s are noisy.
- **Know the bound before predicting a gain** (Measuring 6, 9, 16). After every GPU step, find what limits the run:
  - GPU busy % from an `nsys` timeline, in a separate untimed run;
  - host per-thread busy time;
  - the reader;
  - per-pass idle time.

### 8.2 Workload sets

`dev/perf/workloads.tsv` (`ab.sh --list --set all` shows them all):

| Set | Workloads | Metric (`--gpu` comparisons use wall time) |
|---|---|---|
| `cpu` | W1 hit-rich ×10k, W2 length mix ×10k, W3 M200 × Swiss-Prot, W4 long sequences, W5 `--max` (all `--cpu 0`) | cycles (CPU neutrality) |
| `reader` | R1 M50 × Swiss-Prot, `--cpu 0` | cycles_uk |
| `scale` | X1 Pfam 100 × Swiss-Prot 50k @ 4/8/16 (+ no `-Z`), X2 M200 × Swiss-Prot @ 16 | wall |
| `query` | X3 Pfam 1k × Swiss-Prot 1k @ 0/16 (+ `--dbcache`), X4 Pfam-A text / `.h3m` × 5 seqs @ 0/16 (+ `--dbcache`) | wall |
| `milestone` | X5 Pfam 1k × human `--cut_ga` @ 16 (+ `--dbcache`) | wall |
| `annot`, `annot_txt` (`--file annot.tsv`) | All of Pfam-A `--cut_ga` × 6 UniProt taxa × 4 length groups of 250 (`A_<set>_L<g>`) | wall → **sequences/s** |

Add a `gpu` set for the cases a GPU is meant for, which today's sets under-represent. Pin `-Z` wherever E-value thresholds are used, so threaded outputs are deterministic. Suggested:

- `G1`: Pfam 1k × all of Swiss-Prot, `--cut_ga` @ 16 (big database, realistic hits);
- `G2`: pressed Pfam-A × Swiss-Prot 50k, `--cut_ga` @ 16 (many models × medium database);
- `G3`: Pfam 100 × `random_100k` @ 16, `-Z 100000` (few hits: the feed-bound case);
- `G4`: M200 × Swiss-Prot @ 16 = X2 (one model × a big database).

Commit the new workloads in a separate "coverage only" commit, with their CPU outputs checked for determinism first.

### 8.3 Comparisons to run at each milestone

1. **CPU neutrality.** §7.1 item 5.
2. **GPU gain.** Same binary, `ab.sh <bin> <bin> --b-args=--gpu --set scale,query,milestone,gpu --reps 5`, and `--file annot.tsv`.
3. **Host-core dependence.** X5, G1 and `A_human_L2` with `--gpu` at `--cpu` 2, 4, 8 and 16. This shows how many host cores one GPU can use: the balance point from §4.
4. **Kernel microbenchmarks.** GCUPS per kernel by M and L, with and without DPX, plus achieved occupancy and the limiting pipe from `ncu`.
5. **Resources.**
   - Host peak RSS (`ab.sh` records it).
   - Peak device memory. Print it under a debug environment variable to stderr only in non-compared runs, or read it from `nvidia-smi --query-compute-apps`.
   - CUDA init time.
   - GPU clocks, power and temperature during the runs: `nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv -l 1`.
6. **Process sharing.** 1×16, 2×8 and 4×4 processes on one GPU, with X5 query shards (`dev/profiling/reprofile/shard_queries.sh`).

### 8.4 Acceptance per step (`PLAN.md` §3, adapted)

- **Correctness:** all of §7 passes.
- **Speed:** `faster` on the step's target workloads, and nothing `SLOWER` elsewhere in its sets. Tiny workloads (X2, X4, X3 @ 0) may be slower with `--gpu` because of context creation. Report that honestly and keep it as small as possible; don't hide it with stderr notes, because stderr is compared.
- **Memory:** report any growth in host RSS and device memory, and explain it.

---

## 9. Recording and reporting

1. **`dev/profiling/PROGRESS.md`.** One row per step and per milestone, in the existing columns: step, commit, set, workload, metric, A → B, ratio, A/A spread, notes. Also record dropped steps with their numbers. The existing rows are full of "measured and dropped" entries, and that's expected.
2. **`dev/profiling/GPU.md`**: the design document and full report. It covers:
   - the G0 ceilings;
   - the architecture as built;
   - per-kernel exactness arguments and their evidence (utests, mutants, crosscheck);
   - kernel GCUPS tables;
   - end-to-end results;
   - host/GPU balance;
   - projections;
   - limitations and open questions;
   - exact commands to reproduce.
3. **`dev/profiling/LEARNINGS.md`**: a new section "GPU path (`--gpu`)". Use the same style: numbered items, each a concrete lesson with the number that taught it. Add machine or tool gotchas to "Tooling gotchas".
4. **`dev/profiling/RESULTS.md`**: a new section "GPU offload (`--gpu`)". RESULTS.md is the user's **single** results page.
   - Follow the page's existing structure and tone:
     - headline bullets with numbers;
     - an SVG chart in light and dark (`<picture>` with `prefers-color-scheme`), made with `series.py` + `plot_series.py`/`plot_annot.py`, extended as needed (stdlib only, no matplotlib);
     - a full table of medians;
     - "How this was measured", with machine, binaries (commit SHAs), rounds and noise;
     - reproduce commands.
   - Say "outputs identical" only for what was actually compared.
   - `series.py` runs binaries, not options. Add per-binary arguments, or a wrapper that runs as `hmmsearch` with `--gpu` (argv[0] must stay `hmmsearch`).
   - **Ask the user which chart type and metric before drawing.** Known preference: annotation results go as one grouped bar chart per taxon, length groups as the groups, throughput in sequences/s. Propose three bars per group: upstream 3.4, CPU `dev`, CPU + `--gpu`.
5. **Projections to the RTX PRO 6000.**
   - Label them as projections. State the method:
     - find the limiting pipe of each kernel with `ncu` on the H100;
     - scale by SM count × clock × that pipe's per-SM rate on the target architecture;
     - give a range.
   - Cross-check against the DPX-free build.
   - End-to-end projections must also account for the host: the target nodes' CPU core counts set the CPU share, and when targets aren't resident on the device, the PCIe 5 host link bounds the feed.
   - Recommend a validation run on a real RTX PRO 6000 before any deployment decision. The user can arrange access.
6. **Final message to the user.**
   - What was built.
   - Measured results (H100), with the ceiling next to the achieved speedup.
   - Projections, labelled.
   - Every check that ran, and anything skipped.
   - The commits and branches waiting for the user's push.
   - Open decisions (§1.6).

Writing style for all documents: plain, short declarative sentences, numbers with units, tables for data, no marketing language. Match the existing pages.

---

## 10. Definition of done

- [ ] `./configure --enable-cuda` builds a fatbin for sm_90 and sm_120; the default build is unchanged (instructions 1.0000).
- [ ] `hmmsearch --gpu` runs at least SSV/MSV on the GPU. The Viterbi filter and further stages are there if they pay, by measurement.
- [ ] Outputs are byte-identical on:
  - [ ] all goldens (`--force-opt --gpu`);
  - [ ] all workload sets, including `annot`, `annot_txt` and the new `gpu` set;
  - [ ] the crosscheck build;
  - [ ] the error paths;
  - [ ] the stress loop.
- [ ] Kernel utests with mutants and negative controls are in `make check`.
- [ ] TSan, ASan and `compute-sanitizer` are clean.
- [ ] PROGRESS.md rows, GPU.md, the LEARNINGS section and the RESULTS.md section are written; the chart type was confirmed with the user.
- [ ] Measured H100 results and labelled RTX PRO 6000 projections, each next to the measured G0 ceiling.
- [ ] Everything is committed on `perf/gpu*` branches; `dev` is fast-forwarded after the milestones; nothing is pushed; the user has the list of what to push.
