# hmmsearch acceleration — measured progress

One row per step of [`PLAN.md`](PLAN.md). A is the snapshot before the step, B is the step. Numbers come from `dev/perf/ab.sh` (median of 5 interleaved reps). "Spread" is the larger (max−min)/median of A's and B's reps. Every row also passed `dev/check.sh` and the bench hashes, and `ab.sh` reported identical outputs.

| Step | Commit | Set | Workload | Metric | A → B | B/A | Spread | Notes |
|---|---|---|---|---|---|---|---|---|
| 0 (A/A) | f8bebd8 | cpu | W1-W5 | cycles:u | — | 0.998-1.002 | 0.38-0.75% | Noise calibration; instructions:u ratio 1.0000 |
| 0 (A/A) | f8bebd8 | scale, query | X1-X4 | wall | — | 0.988-1.047 | 1.4-5.8% | Noise calibration; outputs identical with `-Z` pinned on threaded runs |
