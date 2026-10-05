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
