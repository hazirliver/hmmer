# GPU path for hmmsearch: G0 findings and decision (2026-10-09)

Status: **stopped at the G0 gate. No GPU code was written.** The task prompt is `GPU-AGENT-PROMPT.md`; tools and raw results are in `gpu/`. The user decided to test larger CPU VMs instead.

## What was measured

A throwaway record/replay patch (`gpu/g0-replay.patch`, never in `src/`) records every (query, target) filter outcome at `--cpu 0`, then replays it at 16 threads with stages skipped:

- **A**: SSV/MSV and the Viterbi filter skipped (what a GPU doing the integer filters would remove);
- **B**: also the bias filter and the Forward parser for F3 failures.

The speedup of A or B over the reference is the ceiling for a free GPU. Outputs equal the reference in every run; a lookup-only control (nothing skipped) was within 1% of the reference on the small sets and 3-6% slower on the big ones.

| Workload (16 threads) | ref | ceiling A | ceiling B |
|---|---:|---:|---:|
| X5 Pfam 1k × human `--cut_ga` | 8.03 s | 1.76× | 2.03× |
| X1 Pfam 100 × Swiss-Prot 50k | 1.26 s | 2.09× | 2.41× |
| X2 M200 × Swiss-Prot | 0.45 s | 0.95× | 0.99× |
| X4 Pfam-A × 5 sequences | 1.46 s | 1.03× | 1.03× |
| W1 / W5 (16 thr) | 0.18 / 2.81 s | 1.11× / 1.02× | 1.14× / 1.01× |
| A_bacteria_L2 / A_human_L1 / A_human_L4 | 2.14 / 2.22 / 14.84 s | 1.14 / 1.09 / 1.15× | 1.15 / 1.11 / 1.23× |
| **Pfam-A × human 20,428** | 193.5 s | 1.62× | 1.84× |
| **Pfam-A × bacteria 50,000** (15.7 M res) | 177.9 s | 2.62× | 3.34× |

## Lessons

1. The profile-share estimate of the ceiling (4.8× realistic) was about 2.4× too high. Measure the ceiling by replay (LEARNINGS Measuring 15).
2. Small annotation runs (250 proteins × Pfam-A) are bound by per-model overhead and the domain stage: ceiling 1.1-1.25×. Long, hit-rich proteins (human) keep the ceiling near 1.8×, because workers wait on domain-stage tails (workers busy 52% in B). Short bacterial proteins reach 3.3× (60% busy).
3. X2 (one model × big database) and X4 are reader/setup bound: a GPU can't help.
4. Replay needs target identity. `ESL_SQ.idx` is not set by the ASCII reader, so the patch keys on a 64-bit hash of the name (fatal on a collision).
5. Bug seen while building it: `esl_exp_surv()` used in `hmmsearch.c` without `#include "esl_exponential.h"` compiled with only an implicit-declaration warning and returned garbage (int return). It dropped half the hits. Check warnings in throwaway patches.
6. A GPU adds transfer and launch cost to these ceilings, and the feed rate (reader ~1 G residues/s) becomes a bound. The realistic gain is lower than the ceiling.

## Extrapolation to all of GenBank nr vs Pfam-A (estimate, not measured)

Assuming ~6e8 proteins, mostly bacterial-like: 281 proteins/s at 16 vCPU (8 cores × 2 HT) gives ~25 days. A human-like mix (106 proteins/s) would take ~65 days. Scaling by vCPU, assuming 85-95% efficiency and similar cores: 32 vCPU ~12.5 d, 64 ~6.5 d, 128 ~3.3-3.8 d. One GPU (ceiling 3.3×) with a 16-vCPU host would reach ~7.5 d at best, so a 128-vCPU CPU-only machine should beat it.

Not measured: scaling beyond 16 threads. Candidates to limit one big process: the master thread, the in-order output window, the reader, domain-stage tails on long proteins. The best layout is probably several processes, each sharding the database (with `--cut_ga` thresholds don't depend on Z, so shard outputs just concatenate), as the CPU query-sharding results in REPROFILE.md §6 suggest.

## Suggested test on a larger CPU VM

Pfam-A (pressed, `--cut_ga`) × `bacteria_50k.fa` (`gpu/mk_bact50k.py`), with 1×N threads and k processes of N/k threads each, N = 16, 32, 64, 128; compare against 177.9 s at 16 vCPU here. Pin `-Z` only if E-value thresholds are used.

## Not done

Everything after G0 in the prompt (kernels, `--gpu`, build system, utests, RESULTS.md section, RTX PRO 6000 projections).
