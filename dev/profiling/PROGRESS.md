# hmmsearch acceleration — measured progress

One row per step of [`PLAN.md`](PLAN.md). A is the snapshot before the step, B is the step. Numbers come from `dev/perf/ab.sh` (median of 5 interleaved reps). "Spread" is the larger (max−min)/median of A's and B's reps. Every row also passed `dev/check.sh` and the bench hashes, and `ab.sh` reported identical outputs.

| Step | Commit | Set | Workload | Metric | A → B | B/A | Spread | Notes |
|---|---|---|---|---|---|---|---|---|
| 0 (A/A) | f8bebd8 | cpu | W1-W5 | cycles:u | — | 0.998-1.002 | 0.38-0.75% | Noise calibration; instructions:u ratio 1.0000 |
| 0 (A/A) | f8bebd8 | scale, query | X1-X4 | wall | — | 0.988-1.047 | 1.4-5.8% | Noise calibration; outputs identical with `-Z` pinned on threaded runs |
| 1.1 | (dropped) | cpu | W1-W4 | cycles:u (`--perf`) | — | 0.995-0.999 | 0.2-1.1% | Grow-only filter matrix in `P7_BG`: instructions −0.15% to −0.54%, below the floor. Superseded by 1.2. In default builds the same change *looked* +1-3% slower from code layout alone, which is why `snap.sh --perf` exists |
| 1.2 | "Bias filter: streaming 2-state forward" | cpu | W3 m200 × sprot | cycles:u (`--perf`) | 12.642G → 12.102G | **0.957** | 0.73% | Streaming 2-state filter forward. Filter share of CPU 9.4% → 2.8% |
| 1.2 | | cpu | W4 long seqs | cycles:u (`--perf`) | 1.087G → 1.056G | **0.972** | 1.72% | |
| 1.2 | | cpu | W1 hitters | cycles:u (`--perf`) | 2.853G → 2.816G | **0.987** | 0.87% | |
| 1.2 | | cpu | W2 mix | cycles:u (`--perf`) | 10.857G → 10.725G | **0.988** | 0.59% | |
| 1.2 | | cpu | W5 `--max` | cycles:u (`--perf`) | 11.850G → 11.818G | 0.997 (~) | 1.06% | Bias filter isn't used with `--max`; instructions B/A 1.0000 |
| 1.2 | | cpu | W1-W5 | cycles:u (default `-O3` vs f8bebd8) | — | 0.986 / 0.982 / **0.950** / 0.971 / 1.010 | 0.2-1.6% | W5's +1.0% has identical instructions: a layout artifact |
| Phase 1 milestone | 2547ac9 | milestone | X5 pfam_1k × human `--cpu 16` | wall (default `-O3` vs f8bebd8, 3 reps) | 58.73 s → 57.78 s | 0.984 (~) | 2.39% | Instructions −2.7%, but at 16 threads wall time is reader-bound (Phase 2's target). `dev/check.sh --bench M`: all 121 hashes match. Its 0.920 overall ratio vs the 10-03 baseline is mostly machine drift (untouched kernels 0.877), so the baseline was re-saved from this run |
| 2.1 | "hmmsearch --dbcache" | cpu | W1-W5, no `--dbcache` | cycles:u (`--perf`) | — | 0.990-1.013 | 0.2-1.2% | Default path: instructions B/A 1.0000 on all five; W3's +1.3% cycles with identical instructions is a microarchitectural artifact (rechecked in 2.2) |
| 2.1 | | scale | X1 p100 × 50k @ 4 / 8 / 16 / 16 no `-Z` | wall, B with `--dbcache` | 8.28 / 7.98 / 8.98 / 9.48 s → 4.12 / 2.62 / 2.54 / 2.49 s | **0.50 / 0.33 / 0.28 / 0.26** | 3.6-12.7% | Reader re-parse gone: scales past 4 threads. RSS 166 → 130 MB @ 16 |
| 2.1 | | scale | X2 1 query × Swiss-Prot @ 16 | wall | 1.02 s → 1.42 s | **1.39 (slower)** | 5.0% | Single query, so the cache is pure overhead: 570 fresh blocks with default `ESL_SQ` buffers. RSS 114 → 653 MB. Target of 2.2 |
| 2.1 | | query | X3 pfam_1k × sprot_1k @ 0 / 16 | wall | 4.93 / 6.52 s → 3.62 / 4.97 s | 0.74 / 0.76 | 2.3-3.1% | |
| 2.1 | | query | X4 Pfam-A × 5 seqs | wall | — | 0.95-0.99 (~) | 3.7-7.4% | DB is 5 sequences; per-query overhead dominates (Phase 3) |
| 2.2 | "dbcache: compact arena storage" | cpu | W1-W5, no `--dbcache` | cycles:u (`--perf`) | — | 0.993-1.001 (~) | 0.3-1.5% | Default path neutral; W3 now 1.0000 (2.1's +1.3% was noise) |
| 2.2 | | scale | X2 1 query × Swiss-Prot @ 16, `--dbcache` | wall | 1.02 s → 1.17 s | 1.15 (was 1.39) | 5.0% | Arena copy: instructions B/A 1.007. Remaining cost: sys time from first-touching ~430 MB, lowered by huge pages (1.26 → 1.15). RSS 547 MB (was 653) |
| 2.2 | | scale | X1 p100 × 50k @ 8 / 16, `--dbcache` | wall | 8.43 / 9.43 s → 2.72 / 2.60 s | 0.32 / 0.28 | 5.5-11.7% | Unchanged from 2.1 |
| 2.2 | | query | X3 @ 0 / 16, `--dbcache` | wall | 4.92 / 6.47 s → 3.67 / 4.97 s | 0.75 / 0.77 | 3.1-4.1% | |
| 2.3 | "dbcache: --dbcache_max" | — | protdb, 1 MB limit | functional | — | — | — | Over the limit: warning, results identical to uncached (`--cpu` 0 and 3); a `.gz` target gives a clear fatal error; with 2 MB the `.gz` multi-query works. The help column for "Other expert options" is 4 characters wider (longer option name) |
| Phase 2 milestone | b7e123a | milestone | X5 pfam_1k × human `--cut_ga` @ 16, `--dbcache` | wall (default `-O3` vs 176cb7e, 3 reps) | 57.33 s → 20.28 s | **0.354** | 3.59% | 2.8× faster on the realistic annotation run, outputs identical; RSS 327 → 253 MB |
| Phase 2 milestone | | scale | X1 p100 × 50k @ 16, default `-Z`, `--dbcache` | wall | 8.93 s → 2.44 s | **0.273** | 5.99% | |
| Phase 2 milestone | | cpu | W1-W5, no option | cycles:u (default `-O3`) | — | 0.988-1.014 | 0.2-1.1% | Instructions identical (1.0000); ±1.4% cycles in both directions is code layout in unaligned builds |
| 3.1 | "hmmsearch: reuse per-thread pipelines and hit lists across queries" | query | X3 pfam_1k × sprot_1k / X4 Pfam-A text / X4 Pfam-A pressed, × 5 seqs, @ 16 | wall (`--perf`) | 6.42 / 47.64 / 31.98 s → 6.27 / 44.68 / 28.68 s | **0.977 / 0.938 / 0.897** | 1.9-2.5% | `p7_pipeline_Reinit()` instead of Create/Destroy per query and thread. Instructions −0.9% / −4.6% / −18%; sys time still 4.2 / 85 / 87 s (thread creation, step 3.3) |
| 3.1 | | query | X3 / X4 text / X4 pressed @ 0 | cycles:u (`--perf`) | — | 1.012 / 1.001 / 0.986 | 0.5-5.8% | X3's +1.2% (+1.9% before the row-stride fix below) is code layout: a control build with B's library and reuse switched off is +1.5% vs A, so reuse itself is −0.5%. B has 13% fewer L1 misses, 58% fewer dTLB misses, 74% fewer page faults |
| 3.1 | | cpu | W1-W5 (single query, reuse never runs) | cycles:u (`--perf`) | — | 0.998-1.009 | 0.2-1.1% | Instructions 1.0000; the same layout shift |
| 3.2 | (deferred) | — | — | profile | — | — | — | In the serial X4 pressed profile, allocation doesn't appear among the top 40 symbols; building each profile does: `p7_oprofile_Convert` 20%, `esl_sse_expf` 19%, `roundf` 9%, `p7_ProfileConfig` 6%. Reusing `gm`/`om` buffers would save < 1% |
| 3.3 | "hmmsearch: persistent worker threads" | query | X3 / X4 text / X4 pressed @ 16 | wall (`--perf`) | 6.43 / 45.36 / 28.73 s → 5.47 / 25.35 / 9.18 s | **0.851 / 0.559 / 0.320** | 0.9-1.6% | Workers created once, parked at a mutex/condvar gate between queries. Sys time 4.2 / 86 / 87 s → 1.1 / 16 / 16 s. TSan clean (`setarch -R`), 100-run stress at 2-32 threads |
| 3.3 | | query | X3 / X4 text / X4 pressed @ 0 | cycles:u (`--perf`) | — | 0.985 / 0.990 / 1.001 | 1.5-3.0% | Instructions 1.0000: the serial path is untouched |
| 3.3b | "hmmsearch: wake only the workers a pass needs" | query | X3 / X4 text / X4 pressed @ 16 | wall (`--perf`) | 5.47 / 25.31 / 9.23 s → 5.47 / 22.72 / 6.78 s | 1.000 / **0.898 / 0.735** | 0.7-2.7% | Workers per pass = blocks in the previous pass, one condvar per worker. Sys time 0.97 / 15.9 / 16.2 s → 0.17 / 1.1 / 0.85 s. X3 is one block, so it's unchanged (3.3c). TSan clean on 1-block and 3-block dbs; 100-run stress |
