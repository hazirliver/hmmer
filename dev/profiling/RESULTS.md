# hmmsearch speedups: results

What the work in [`PLAN.md`](PLAN.md) achieved. Upstream HMMER 3.4 and the end of every phase ran on the same workloads in one session (2026-10-06). **Every speedup comes with byte-identical output**: in each scenario, every binary's output files were compared with upstream's (dates and timings masked), and all matched. Per-step details are in [`PROGRESS.md`](PROGRESS.md), lessons in [`LEARNINGS.md`](LEARNINGS.md).

**Headline:**

- Realistic annotation, Pfam (1,000 families) against the human proteome with `--cut_ga`, 16 threads: **59 s → 21 s (2.8×)** by default, **11 s (5.4×)** with the new opt-in `--dbcache`.
- All of Pfam-A (22k pressed HMMs) against a few sequences, 16 threads: **35 s → 2.9 s (12×)**, 1.5 s (23×) with `--dbcache`.
- Single-threaded runs: **1.2-2.0×**, depending on where the time goes (filters and parsing gained most, the domain stage least).

## Speedup over upstream, per scenario

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/speedup-dark.svg">
  <img alt="Horizontal bar chart: speedup of this branch over upstream HMMER 3.4 for eleven scenarios, from 1.2× (--max, domain stage only) to 12× (all of Pfam-A against 5 sequences); with --dbcache 5.4×, 23× and 7.2× on the three many-model scenarios. Values are in the table below." src="results/speedup-light.svg">
</picture>

`--dbcache` keeps the parsed target database in memory across queries (Phase 2). It is opt-in, so its bars apply only where you pass it; it pays off when one run searches many queries.

## Wall time at the end of every phase

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/phases-dark.svg">
  <img alt="Small multiples: wall time of each scenario at upstream 3.4, at the end of phases 1 to 6, and now. The large drops: P3 for many models against small databases, P5 for runs that parse big FASTA databases, P6 for filter-heavy runs, and 'now' (FASTA residue table) for uncached threaded runs." src="results/phases-light.svg">
</picture>

`3.4` is upstream, `P1`-`P6` the end of each phase, `now` the current `dev` (Phase 6 plus step 5.6, a faster FASTA residue path done after Phase 6).

| Phase | What changed | Where it shows |
|---|---|---|
| 1. Bias filter | Streaming 2-state forward instead of a per-call matrix | 4-5% where the bias filter runs often (1 model × Swiss-Prot, the realistic run) |
| 2. Target database once | Opt-in `--dbcache`: the first query's pass keeps the parsed targets in memory | Multi-query runs with the option: 3.0× on the realistic run |
| 3. Per-query costs | Persistent worker threads, reused pipelines, queries prepared ahead in parallel, query-parallel mode for small cached databases | Many models × small databases: 12× on Pfam-A × 5 sequences |
| 4. Domain stage | Exact rewrites of the stochastic traceback, null2 and clustering hot spots | Hit-rich and `--max` runs: 5-14% |
| 5. FASTA parser (Easel) | 16 KB reads, one-pass validate + digitize, header scans with `memchr` | Runs that parse a large FASTA: realistic run −46%, Pfam 100 × Swiss-Prot 50k −51% |
| 6. Wider SIMD | AVX2 SSV, MSV and Viterbi filters (AVX-512 for models longer than 448), chosen at run time | Filter-heavy runs: 15-27% |
| 5.6 (after 6) | FASTA residues digitized through one table, 8 at a time | Uncached threaded runs: 14-18% |

## All numbers

Median wall time of 3 rounds (each round ran every binary once, in rotated order, after one untimed round that checked outputs).

| Scenario | upstream 3.4 | P1 | P2 | P3 | P4 | P5 | P6 | now | speedup |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pfam 1k × human proteome | 59.0 s | 56.7 s | 56.5 s | 52.7 s | 50.6 s | 27.3 s | 25.5 s | 21.4 s | 2.8× |
| Pfam 1k × human proteome (--dbcache) | — | — | 19.1 s | 15.2 s | 13.3 s | 13.3 s | 11.2 s | 11.0 s | 5.4× |
| All of Pfam-A × 5 sequences | 35.0 s | 34.9 s | 34.9 s | 2.9 s | 2.9 s | 2.9 s | 2.9 s | 2.9 s | 12× |
| All of Pfam-A × 5 sequences (--dbcache) | — | — | 33.4 s | 1.5 s | 1.5 s | 1.5 s | 1.5 s | 1.5 s | 23× |
| Pfam 1k × Swiss-Prot 1k | 6.7 s | 6.6 s | 6.4 s | 4.2 s | 4.2 s | 3.5 s | 2.9 s | 2.7 s | 2.5× |
| Pfam 1k × Swiss-Prot 1k (--dbcache) | — | — | 5.1 s | 0.98 s | 0.99 s | 0.97 s | 0.94 s | 0.93 s | 7.2× |
| Pfam 100 × Swiss-Prot 50k | 9.5 s | 9.2 s | 9.1 s | 8.7 s | 8.6 s | 4.2 s | 4.0 s | 3.3 s | 2.9× |
| 1 model × Swiss-Prot, 16 threads | 1.1 s | 0.96 s | 0.95 s | 0.94 s | 0.91 s | 0.54 s | 0.52 s | 0.44 s | 2.4× |
| 40 models × Swiss-Prot 10k | 3.5 s | 3.5 s | 3.3 s | 3.4 s | 3.3 s | 2.9 s | 2.1 s | 2.1 s | 1.7× |
| 1 model × Swiss-Prot, 1 thread | 4.0 s | 3.8 s | 3.9 s | 3.7 s | 3.5 s | 3.0 s | 2.5 s | 2.5 s | 1.6× |
| 1 short model × Swiss-Prot | 1.2 s | 1.1 s | 1.2 s | 1.1 s | 1.2 s | 0.71 s | 0.71 s | 0.63 s | 2.0× |
| 10 models × long sequences | 0.36 s | 0.36 s | 0.36 s | 0.35 s | 0.34 s | 0.27 s | 0.22 s | 0.20 s | 1.7× |
| 10 hit-rich models × 10k | 0.87 s | 0.88 s | 0.85 s | 0.85 s | 0.81 s | 0.72 s | 0.69 s | 0.69 s | 1.3× |
| 10 models × 1k, --max | 3.6 s | 3.6 s | 3.6 s | 3.6 s | 3.1 s | 3.0 s | 3.1 s | 3.0 s | 1.2× |

Peak memory: the realistic run takes 331 MB upstream and 386 MB now (per-thread state kept across queries), 337 MB with `--dbcache`.

## How this was measured

- Machine: 16-vCPU Intel Xeon Platinum 8468 (Sapphire Rapids) KVM guest, 8 cores × 2 hyperthreads.
- Binaries: default `-O3` builds of upstream (the HMMER 3.4 source, f8bebd8) and of each phase's last commit (4088a0d, 5ad2fde, fa18adc, 4f62b07, f9dcb6b, aadf224, 7d637cb), each with the Easel commit it pins. Upstream differs from `master` only in `dev/` and `.gitignore`.
- Workloads: [`dev/perf/workloads.tsv`](../perf/workloads.tsv) (Pfam 37.0, Swiss-Prot, the human proteome; see [`dev/bench`](../bench)). Threaded runs pin `-Z`, so their outputs are deterministic.
- Wall times on this VM vary: the spread (max-min)/median of a binary's 3 rounds was up to 12% on the shortest or 16-thread runs, typically 3-6%. That doesn't affect the conclusions here; per-step decisions in PROGRESS.md used CPU counters and interleaved A/Bs.
- Reproduce (about one hour):

```
for p in p0:f8bebd8 p1:4088a0d p2:5ad2fde p3:fa18adc p4:4f62b07 p5:f9dcb6b p6:aadf224 now:7d637cb; do
  dev/perf/snap.sh --rev ${p#*:} --name series-${p%%:*}
done
dev/perf/series.py series-p0 series-p1 series-p2 series-p3 series-p4 series-p5 series-p6 series-now \
    --labels upstream,P1,P2,P3,P4,P5,P6,final --reps 3 --out phases.json \
    --workloads '^(W[1-5]_|R1_|X1_p100_50k_cN$|X2_|X3_p1k_s1k_cN|X4_pfamA_h3m_cN|X5_)'
dev/perf/plot_series.py phases.json --outdir dev/profiling/results --xlabels "3.4,P1,P2,P3,P4,P5,P6,now"
```

The raw samples (wall, user, sys, cycles, instructions, peak RSS per run) are in [`results/phases.json`](results/phases.json).
