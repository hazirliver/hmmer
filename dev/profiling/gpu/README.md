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

## Large-database annotation (all of Pfam-A, pressed, `--cut_ga`, 16 threads, 2 reps)

`g0-big.tsv`, `g0-run-big.sh`; the bacterial set is made by `mk_bact50k.py` (50,000 random non-fragment bacterial Swiss-Prot proteins, seed 1, 15.7 M residues). Outputs identical in every variant. The lookup-only control (ctl) was 3-6% slower than ref here, so the ceilings are slightly conservative.

| Workload | ref | ctl | A | B | ceiling A | ceiling B |
|---|---:|---:|---:|---:|---:|---:|
| Pfam-A × human 20,428 | 193.5 s | 199.8 | 119.4 | 105.1 | 1.62× | 1.84× |
| Pfam-A × bacteria 50,000 | 177.9 s | 188.2 | 68.0 | 53.2 | 2.62× | 3.34× |

In B, the 16 workers' CPU time is 52% (human) and 60% (bacteria) of 16 × wall; the rest is waiting for work.
