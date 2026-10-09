# GPU work (throwaway G0 tools)

`g0-replay.patch` (against 5935bd9, never applied to `src/`) records every (query, target) filter outcome at `--cpu 0` and replays it: variant A skips SSV/MSV and the Viterbi filter, variant B also skips the bias filter and the Forward parser for F3 failures. `g0.py` interleaves ref / ctl (lookups only) / A / B at 16 threads and checks outputs against ref. `g0-run.sh` runs the task's workload list; `results/` holds its JSON.

Measured 2026-10-09 (wall medians, 3 reps, 16 threads, outputs identical in every run):

| Workload | ref | A | B | ceiling A | ceiling B |
|---|---:|---:|---:|---:|---:|
| X5 Pfam 1k × human `--cut_ga` | 8.03 s | 4.55 | 3.96 | 1.76× | 2.03× |
| X1 Pfam 100 × Swiss-Prot 50k | 1.26 | 0.60 | 0.52 | 2.09× | 2.41× |
| X2 M200 × Swiss-Prot | 0.45 | 0.47 | 0.46 | 0.95× | 0.99× |
| X4 Pfam-A × 5 seqs (pressed) | 1.46 | 1.42 | 1.41 | 1.03× | 1.03× |
| W1 hit-rich ×10k (16 thr) | 0.18 | 0.16 | 0.15 | 1.11× | 1.14× |
| W5 `--max` (16 thr) | 2.81 | 2.76 | 2.80 | 1.02× | 1.01× |
| A_bacteria_L2 | 2.14 | 1.89 | 1.87 | 1.14× | 1.15× |
| A_human_L1 | 2.22 | 2.04 | 2.01 | 1.09× | 1.11× |
| A_human_L4 | 14.84 | 12.86 | 12.05 | 1.15× | 1.23× |

ctl (lookups, nothing skipped) was within 1% of ref everywhere, so the harness itself adds no visible cost.
