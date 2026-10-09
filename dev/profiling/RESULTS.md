# hmmsearch speedups: results

What the work in [`PLAN.md`](PLAN.md) achieved. Upstream HMMER 3.4 and the end of every phase ran on the same workloads in one session (2026-10-09). **Every speedup comes with byte-identical output**: in each scenario, every binary's output files were compared with upstream's (dates and timings masked), and all matched. Per-step details are in [`PROGRESS.md`](PROGRESS.md), lessons in [`LEARNINGS.md`](LEARNINGS.md).

**Headline:**

- Realistic annotation, Pfam (1,000 families) against the human proteome with `--cut_ga`, 16 threads: **56 s → 7.9 s (7.1×)** by default, and about the same with the opt-in `--dbcache` (7.7 s, 7.3×).
- All of Pfam-A (22k HMMs) against a few sequences, 16 threads: **24 s → 1.5 s (16×)** from the pressed file, 1.1 s (21×) with `--dbcache`; **41 s → 7.3 s (5.6×)** from the text file.
- Single-threaded runs: **1.3-2.1×**, depending on where the time goes (filters and parsing gained most, the domain stage least).
- Annotating 1,000 Swiss-Prot proteins per taxon (bacteria, phages, fungi, viruses, eukaryotes, human) with all of Pfam-A, 16 threads: **5.2 → 71 sequences per second (14×)** from the pressed file, **3.8 → 28 (7.4×)** from the text file ([below](#annotating-proteins-with-all-of-pfam-a)).

## Speedup over upstream, per scenario

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/speedup-dark.svg">
  <img alt="Horizontal bar chart: speedup of this branch over upstream HMMER 3.4 for twelve scenarios, from 1.3× (--max, domain stage only) to 16× (all of Pfam-A, pressed, against 5 sequences); with --dbcache 7.3×, 21× and 11× on the three many-model scenarios. Values are in the table below." src="results/speedup-light.svg">
</picture>

`--dbcache` keeps the parsed target database in memory across queries (Phase 2). It is opt-in, so its bars apply only where you pass it. Since Phase 7 a run with many queries reads each block of targets once for a whole batch of them, so the option gains little at these database sizes; it still saves reading where the targets fit in memory and the run has thousands of queries.

## Wall time at the end of every phase

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/phases-dark.svg">
  <img alt="Small multiples: wall time of each scenario at upstream 3.4, at the end of phases 1 to 6, after step 5.6 and at the end of phase 7. The large drops: P3 for many models against small databases, P5 for runs that parse big FASTA databases, P6 for filter-heavy runs, 5.6 for uncached threaded runs, and P7 for every threaded run with many queries." src="results/phases-light.svg">
</picture>

`3.4` is upstream, `P1`-`P7` the end of each phase, and `5.6` Phase 6 plus step 5.6 (a faster FASTA residue path, done after Phase 6).

| Phase | What changed | Where it shows |
|---|---|---|
| 1. Bias filter | Streaming 2-state forward instead of a per-call matrix | 4% where the bias filter runs often (1 model × Swiss-Prot, 1 thread) |
| 2. Target database once | Opt-in `--dbcache`: the first query's pass keeps the parsed targets in memory | Multi-query runs with the option: 2.9× on the realistic run |
| 3. Per-query costs | Persistent worker threads, reused pipelines, queries prepared ahead in parallel, query-parallel mode for small cached databases | Many models × small databases: 9.0× on Pfam-A × 5 sequences (19× with `--dbcache`) |
| 4. Domain stage | Exact rewrites of the stochastic traceback, null2 and clustering hot spots | Hit-rich and `--max` runs: 9-14% |
| 5. FASTA parser (Easel) | 16 KB reads, one-pass validate + digitize, header scans with `memchr` | Runs that parse a large FASTA: realistic run −45%, Pfam 100 × Swiss-Prot 50k −50% |
| 6. Wider SIMD | AVX2 SSV, MSV and Viterbi filters (AVX-512 for models longer than 448), chosen at run time | Filter-heavy runs: 16-29% |
| 5.6 (after 6) | FASTA residues digitized through one table, 8 at a time | Uncached threaded runs: 13-19% |
| 7. Query batches, paired parsers | Up to 16 queries per pass over the targets, with consecutive batches overlapped; an exact fast `atof()` for text HMM files; the Forward/Backward parsers on two targets at once; a shorter bias-filter chain | Threaded runs with many queries: realistic 2.7×, Pfam 1k × Swiss-Prot 1k 4.5×, Pfam-A from text 2.4×; single-threaded runs 5-15%, except the reading-bound one (1%) |

## All numbers

Median wall time of 3 rounds (each round ran every binary once, in rotated order, after one untimed round that checked outputs).

| Scenario | upstream 3.4 | P1 | P2 | P3 | P4 | P5 | P6 | 5.6 | P7 | speedup |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pfam 1k × human proteome | 56.4 s | 56.4 s | 55.4 s | 52.2 s | 49.8 s | 27.3 s | 25.8 s | 21.0 s | 7.9 s | 7.1× |
| Pfam 1k × human proteome (--dbcache) | — | — | 19.5 s | 15.6 s | 13.7 s | 13.7 s | 11.5 s | 11.5 s | 7.7 s | 7.3× |
| All of Pfam-A × 5 sequences | 23.6 s | 23.0 s | 23.2 s | 2.6 s | 2.6 s | 2.2 s | 2.2 s | 2.2 s | 1.5 s | 16× |
| All of Pfam-A × 5 sequences (--dbcache) | — | — | 21.8 s | 1.2 s | 1.2 s | 1.2 s | 1.1 s | 1.1 s | 1.1 s | 21× |
| All of Pfam-A (text) × 5 sequences | 40.6 s | 40.9 s | 40.4 s | 17.9 s | 18.0 s | 17.7 s | 17.6 s | 17.6 s | 7.3 s | 5.6× |
| Pfam 1k × Swiss-Prot 1k | 6.1 s | 6.0 s | 6.0 s | 4.1 s | 4.1 s | 3.4 s | 2.8 s | 2.6 s | 0.58 s | 11× |
| Pfam 1k × Swiss-Prot 1k (--dbcache) | — | — | 4.6 s | 0.96 s | 0.96 s | 0.98 s | 0.93 s | 0.93 s | 0.54 s | 11× |
| Pfam 100 × Swiss-Prot 50k | 9.5 s | 9.2 s | 9.0 s | 8.6 s | 8.4 s | 4.2 s | 3.9 s | 3.2 s | 1.2 s | 7.8× |
| 1 model × Swiss-Prot, 16 threads | 1.0 s | 0.93 s | 0.98 s | 0.95 s | 0.92 s | 0.51 s | 0.50 s | 0.43 s | 0.44 s | 2.4× |
| 40 models × Swiss-Prot 10k | 3.5 s | 3.5 s | 3.5 s | 3.4 s | 3.4 s | 3.0 s | 2.1 s | 2.1 s | 1.9 s | 1.8× |
| 1 model × Swiss-Prot, 1 thread | 4.1 s | 3.9 s | 4.0 s | 3.9 s | 3.6 s | 3.1 s | 2.6 s | 2.5 s | 2.4 s | 1.7× |
| 1 short model × Swiss-Prot | 1.2 s | 1.1 s | 1.1 s | 1.2 s | 1.2 s | 0.71 s | 0.71 s | 0.62 s | 0.61 s | 2.0× |
| 10 models × long sequences | 0.38 s | 0.38 s | 0.37 s | 0.36 s | 0.36 s | 0.28 s | 0.23 s | 0.21 s | 0.18 s | 2.1× |
| 10 hit-rich models × 10k | 0.92 s | 0.91 s | 0.91 s | 0.90 s | 0.82 s | 0.73 s | 0.70 s | 0.69 s | 0.61 s | 1.5× |
| 10 models × 1k, --max | 3.6 s | 3.7 s | 3.7 s | 3.7 s | 3.2 s | 3.2 s | 3.1 s | 3.1 s | 2.8 s | 1.3× |

Peak memory: the realistic run takes 341 MB upstream and 350 MB now (per-thread state kept across queries, and the hit lists of a batch of up to 16 queries held until they're printed), 352 MB with `--dbcache`. Many models against a tiny database take 35 → 62 MB.

## Annotating proteins with all of Pfam-A

Six sets of 1,000 UniProtKB/Swiss-Prot proteins: bacteria, phages, fungi, viruses of eukaryotes, eukaryotes other than human, and human. Each set is split by length into four groups of 250, from L1 (shortest) to L4 (longest), at the set's own quartiles. Every group was searched with all 21,979 Pfam-A models (`--cut_ga`, 16 threads), by upstream 3.4 and by this branch at the end of Phase 7. Speed is sequences per second: 250 divided by the median wall time of 3 rounds. In all 48 runs (24 groups, from the pressed and from the text Pfam-A) the two binaries' outputs were identical. How the sets were drawn: [`dev/bench/README.md`](../bench/README.md), "UniProt protein sets".

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/annot-pressed-dark.svg">
  <img alt="Six grouped bar charts, one per taxon, with length groups L1 to L4: sequences per second searched against the pressed Pfam-A, upstream HMMER 3.4 against this branch. Upstream does 2.2-9.7 sequences per second, this branch 18-155; 14-17× in most groups, 8.2× for the longest human proteins. Values are in the table below." src="results/annot-pressed-light.svg">
</picture>

- **All 6,000 proteins: 5.2 → 71 sequences per second (14×).** 21 of the 24 groups are 14-17× faster.
- Upstream is slow on these runs because 250 targets fit in one of its read blocks (1,000 sequences): every model is searched by a single worker, and the run keeps 1.4-2.6 cores busy. This branch spreads each batch of up to 16 models over all the workers (Phase 7) and keeps 5.7-9.9 cores busy.
- The longest eukaryotic proteins gain least: human L4 (696-13,477 residues) 8.2×, the other eukaryotes' L4 11×, human L3 12×. Human L4 keeps the fewest cores busy (5.7): each pass over a batch of models probably waits for the worker that has the longest proteins.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="results/annot-text-dark.svg">
  <img alt="The same six grouped bar charts for the text Pfam-A file: upstream 1.8-5.9 sequences per second, this branch 13-34, 5.7-10× faster. Values are in the table below." src="results/annot-text-light.svg">
</picture>

- From the text file: **3.8 → 28 sequences per second (7.4×)**. Reading the 21,979 models from text adds 4.3-6.1 s to each of this branch's runs (median 5.6 s; upstream 14-24 s). That caps it near 30 sequences per second whatever the protein length, with 2.2-6.4 cores busy. Pressing Pfam-A once with `hmmpress` removes it: `hmmsearch` reads `Pfam-A.hmm.h3m` instead of the text when it sits next to `Pfam-A.hmm`.

Sequences per second (median of 3 rounds):

| Set | Group | Lengths | pressed: upstream | this branch | speedup | text: upstream | this branch | speedup |
|---|---|---|---:|---:|---:|---:|---:|---:|
| Bacteria | L1 | 29-165 | 9.0 | 144 | 16× | 5.7 | 34 | 5.9× |
|  | L2 | 166-277 | 7.4 | 118 | 16× | 5.0 | 33 | 6.7× |
|  | L3 | 278-407 | 6.4 | 103 | 16× | 4.7 | 33 | 7.1× |
|  | L4 | 411-1,457 | 4.6 | 73 | 16× | 3.5 | 30 | 8.7× |
| | **all 1,000** | | **6.5** | **103** | **16×** | **4.6** | **32** | **7.1×** |
| Phages | L1 | 23-93 | 9.7 | 155 | 16× | 5.9 | 33 | 5.7× |
|  | L2 | 93-164 | 8.6 | 141 | 16× | 5.5 | 33 | 6.0× |
|  | L3 | 164-336 | 7.2 | 119 | 17× | 4.9 | 33 | 6.7× |
|  | L4 | 336-3,500 | 4.7 | 70 | 15× | 3.6 | 29 | 8.0× |
| | **all 1,000** | | **7.0** | **110** | **16×** | **4.8** | **32** | **6.6×** |
| Fungi | L1 | 15-256 | 8.1 | 125 | 15× | 5.0 | 31 | 6.2× |
|  | L2 | 256-400 | 6.3 | 103 | 16× | 4.2 | 30 | 7.1× |
|  | L3 | 400-593 | 5.2 | 84 | 16× | 3.7 | 29 | 7.9× |
|  | L4 | 595-2,904 | 3.0 | 41 | 14× | 2.5 | 23 | 9.5× |
| | **all 1,000** | | **5.0** | **74** | **15×** | **3.6** | **28** | **7.8×** |
| Eukaryotic viruses | L1 | 31-179 | 8.6 | 137 | 16× | 5.3 | 32 | 6.0× |
|  | L2 | 179-318 | 6.8 | 112 | 16× | 4.7 | 31 | 6.6× |
|  | L3 | 318-570 | 5.5 | 90 | 16× | 3.9 | 29 | 7.5× |
|  | L4 | 570-7,182 | 2.7 | 40 | 15× | 2.2 | 23 | 10× |
| | **all 1,000** | | **4.9** | **77** | **16×** | **3.6** | **28** | **7.9×** |
| Non-human eukaryotes | L1 | 6-200 | 8.7 | 136 | 16× | 5.3 | 31 | 6.0× |
|  | L2 | 202-355 | 6.5 | 89 | 14× | 4.3 | 30 | 7.0× |
|  | L3 | 355-537 | 5.5 | 85 | 16× | 3.9 | 29 | 7.6× |
|  | L4 | 537-5,801 | 2.7 | 29 | 11× | 2.2 | 19 | 8.7× |
| | **all 1,000** | | **4.9** | **62** | **13×** | **3.6** | **27** | **7.5×** |
| Human | L1 | 16-246 | 8.0 | 116 | 15× | 5.1 | 32 | 6.3× |
|  | L2 | 248-395 | 6.2 | 91 | 15× | 4.3 | 30 | 7.1× |
|  | L3 | 396-695 | 4.5 | 55 | 12× | 3.3 | 26 | 8.0× |
|  | L4 | 696-13,477 | 2.2 | 18 | 8.2× | 1.8 | 13 | 7.4× |
| | **all 1,000** | | **4.1** | **42** | **10×** | **3.1** | **22** | **7.3×** |
| **All six** | **6,000** | | **5.2** | **71** | **14×** | **3.8** | **28** | **7.4×** |

Peak memory grows with protein length here: 34-192 MB upstream and 65-754 MB with this branch, the top of both for human L4 (pressed or text alike). Since Phase 3 every worker keeps its DP matrices across models, sized for the largest domain it has met; upstream frees them after each model and, on these runs, has one worker searching at a time.

## How this was measured

- Machine: 16-vCPU Intel Xeon Platinum 8468 (Sapphire Rapids) KVM guest, 8 cores × 2 hyperthreads.
- Binaries: default `-O3` builds of upstream (the HMMER 3.4 source, f8bebd8) and of each phase's last commit (4088a0d, 5ad2fde, fa18adc, 4f62b07, f9dcb6b, aadf224, 7d637cb for `5.6`, 50882a6 for `P7`), each with the Easel commit it pins. Upstream differs from `master` only in `dev/` and `.gitignore`.
- Workloads: [`dev/perf/workloads.tsv`](../perf/workloads.tsv) (Pfam 37.0, Swiss-Prot, the human proteome; see [`dev/bench`](../bench)). Threaded runs pin `-Z`, so their outputs are deterministic.
- Wall times on this VM vary: the spread (max-min)/median of a binary's 3 rounds was 2.7% at the median and up to 12% (the 2 s pressed Pfam-A run), and upstream's realistic run measured 54-59 s on other days. That doesn't affect the conclusions here; per-step decisions in PROGRESS.md used CPU counters and interleaved A/Bs.
- Threaded scenarios are timed without `perf stat`; the single-threaded ones (`W1`-`W5`, `R1`) also collect counters with it. Until 2026-10-07 every timed round ran under `perf stat`, which slows threaded runs that synchronize often on this VM, upstream's most ([`LEARNINGS.md`](LEARNINGS.md), measuring item 20). The previous version of this page overstated two speedups that way: all of Pfam-A × 5 sequences took 34.4 s upstream under `perf stat` and takes 23.6 s without (20× → 16×), and from the text file 50.8 s → 40.6 s (6.6× → 5.6×). The other scenarios moved by less than the day-to-day noise.
- Annotation runs: the six protein sets come from UniProt release 2026_03 (`dev/bench/uniprot.py fetch` downloads exactly those proteins; see [`dev/bench/README.md`](../bench/README.md)), and the workloads are [`dev/perf/annot.tsv`](../perf/annot.tsv). Upstream and `P7` ran round-robin, 1 untimed round that checked outputs and 3 timed rounds, without `perf stat` (it slowed threaded runs, upstream's most; see [`LEARNINGS.md`](LEARNINGS.md), measuring item 20). The pressed runs and the bacteria and phage text runs were measured on 2026-10-07, the other text runs on 2026-10-09 (a VM reboot stopped the first run). Within a group, both binaries ran in the same session. The spread (max-min)/median of a binary's 3 rounds was 1.7% at the median, at most 6.1%. Raw samples: [`results/annot.json`](results/annot.json).
- Reproduce (about one and a half hours; the second `series.py` run adds the text Pfam-A scenario, and the two JSON files are merged by concatenating their `workloads` lists):

```
for p in p0:f8bebd8 p1:4088a0d p2:5ad2fde p3:fa18adc p4:4f62b07 p5:f9dcb6b p6:aadf224 p56:7d637cb p7:50882a6; do
  dev/perf/snap.sh --rev ${p#*:} --name series-${p%%:*}
done
S="series-p0 series-p1 series-p2 series-p3 series-p4 series-p5 series-p6 series-p56 series-p7"
L=upstream,P1,P2,P3,P4,P5,P6,5.6,P7
dev/perf/series.py $S --labels $L --reps 3 --out phases-a.json \
    --workloads '^(W[1-5]_|R1_|X1_p100_50k_cN$|X2_|X3_p1k_s1k_cN|X4_pfamA_h3m_cN|X5_)'
dev/perf/series.py $S --labels $L --reps 3 --out phases-b.json --workloads '^X4_pfamA_txt_cN$'
dev/perf/plot_series.py phases.json --outdir dev/profiling/results --xlabels "3.4,P1,P2,P3,P4,P5,P6,5.6,P7"

# annotation runs (about 4 hours)
dev/bench/uniprot.py fetch
dev/perf/series.py series-p0 series-p7 --labels upstream,P7 --file annot.tsv --reps 3 --out annot.json
dev/perf/plot_annot.py annot.json --outdir dev/profiling/results
```

The raw samples (wall, user, sys, cycles, instructions, peak RSS per run) are in [`results/phases.json`](results/phases.json).
