# Upstream HMMER vs pyhmmer vs this fork: plan and findings (2026-10-20)

[pyhmmer](https://github.com/althonos/pyhmmer) (M. Larralde) is a Cython binding to HMMER 3.4. It says that it is often faster than the `hmmsearch` binary with several CPUs. This page plans a comparison of the three on the same jobs and records what came out. The charts and tables are in [`RESULTS.md`](RESULTS.md).

## What pyhmmer does differently

Read from the 0.12.3 sources (`pyhmmer/hmmer/_hmmsearch.py`, `_base.py`, `plan7.pyx`):

- **The same pipeline.** `Pipeline.search_hmm()` runs HMMER's own `p7_pli_NewModel`, `p7_pli_NewSeq`, `p7_bg_SetLength`, `p7_oprofile_ReconfigLength` and `p7_Pipeline` loop on each target, then `p7_tophits_SortBySortkey` and `p7_tophits_Threshold`. The vendored HMMER 3.4 and Easel 0.49 carry only two kinds of patches (`patches/`). Some are portability fixes (`uintptr_t` casts, `ctime_s`, header guards). The others drop `static` from readers and helpers, so that the Cython layer can call them. The scores, E-values and alignments should therefore equal upstream's.
- **The same kernels as upstream.** The PyPI wheel's `liblibhmmer.so` is HMMER's SSE implementation, built with SSE4.1 and GCC 14.2. It has no AVX2 or AVX-512 code (no `ymm`/`zmm` registers), so its filters match upstream's speed, not the fork's (Phase 6).
- **Targets in memory.** The recommended use reads the whole sequence file into a `DigitalSequenceBlock` once, and the threads share it. Upstream re-reads and re-parses the file for every query. The fork's opt-in `--dbcache` (Phase 2) does the same as pyhmmer, and since Phase 7 the fork reads each block once for a batch of up to 16 queries.
- **One query per thread** (`parallel="queries"`, the default). Each worker searches one HMM against all the targets, with no barrier between queries. A queue of 2 × threads chores sits in between, and results come back in query order. Upstream's threads share each query's targets instead, in blocks of 1,000 sequences, and wait for each other at the end of every query. The fork's Phase 3 query-parallel mode and Phase 7 batches remove that barrier too.
- **Few queries: split the targets** (`parallel="targets"`). This mode is used only when pyhmmer knows that there are fewer queries than threads. The targets are cut into one residue-balanced chunk per thread, and the threads search each query together before moving to the next one.
- **Threads: physical cores by default.** `cpus=0` means `psutil.cpu_count(logical=False)`, which is 8 on this 8-core × 2-HT VM. The README suggests no more than the physical cores for `hmmsearch`.
- **No main output.** `TopHits.write()` writes `--tblout`, `--domtblout` and Pfam tables through HMMER's own `p7_tophits_TabularTargets` / `TabularDomains`. The human-readable report (`-o`) and `-A` alignments have no writer. Tables come without the `# Program: ...` block that `hmmsearch` appends (`p7_tophits_TabularTail`).
- **Subnormal numbers are kept.** `hmmsearch` calls `impl_Init()`, which sets flush-to-zero and denormals-are-zero in the MXCSR. The worker threads inherit them. pyhmmer leaves that call out on purpose (`plan7.pyx`, module init): it would change the float arithmetic of the whole Python process. So pyhmmer computes with subnormal numbers. Found in the first output check: on X5 it prints a best-domain E-value of `3.2e-315` (and `1.6e-319`, `3.2e-315` in the domain table) where `hmmsearch` prints `0`. Every other field of the 1,599 + 1,984 rows was identical.

## Plan

**Question.** On the jobs that `RESULTS.md` already measures, how fast is pyhmmer compared with upstream 3.4 and with the fork? Do all three produce the same hits?

**Contenders.**

| Label | What | Build |
|---|---|---|
| upstream | HMMER 3.4 `hmmsearch` (f8bebd8, snapshot `series-p0`) | `-O3`, GCC 13.3 |
| pyhmmer | pyhmmer 0.12.3 from PyPI, via `dev/perf/pyhmmsearch.py` (snapshot `pyhmmer-0.12.3`, made by `dev/perf/pyhmmer.sh`) | wheel: SSE4.1, GCC 14.2 |
| fork | this branch at 5935bd9 (Phase 7 + step 7E.1; no code changed after it; snapshot `series-p7e1`) | `-O3`, GCC 13.3 |

**The pyhmmer side** is a command-line driver, `dev/perf/pyhmmsearch.py`. It takes the `hmmsearch` options that the workloads use and runs pyhmmer the way its documentation recommends:

- It reads the targets once into a `DigitalSequenceBlock`, in the HMMs' alphabet (as `hmmsearch` does).
- It streams the HMMs from the file. When there are fewer HMMs than threads, it passes them as a list, so that pyhmmer picks its `targets` mode by itself.
- It writes both tables with `TopHits.write()`, with the column header on the first query only.

Options map one to one:

| `hmmsearch` option | pyhmmer equivalent |
|---|---|
| `--cut_ga` | `bit_cutoffs="gathering"` |
| `-Z`, `-E` | `Z=`, `E=` |
| `--max` | `bias_filter=False, F1=F2=F3=1` (`pli->do_max` is read nowhere else) |
| `--cpu N` | `cpus=N` |
| `--cpu 0` | `cpus=1`: the main thread, as `hmmsearch --cpu 0` |

Any other option is an error, not ignored. The driver writes where its time went (imports, reading the targets, search + writing) to `.pyhmmsearch.json`, and `series.py` stores that with the sample. The wrapper runs `python -I`, and the wall time includes the interpreter's start and the imports (about 0.05 s).

Two variants of the same snapshot, made by `pyhmmer.sh` options:

- **`pyhmmer-0.12.3-c8`** (`--cpus 8`) runs 8 threads wherever a workload asks for more than 0. 8 is the number of physical cores here, and pyhmmer's own default.
- **`pyhmmer-0.12.3-ftz`** (`--ftz`) sets flush-to-zero and denormals-are-zero through `fesetenv()` before anything else, as `hmmsearch` does. This variant is an experiment: it shows what the subnormals cost pyhmmer, and whether they are the only difference from `hmmsearch`. It is not pyhmmer as people run it, so it is not in the charts.

**Same job for all three.**

- The workloads are [`dev/perf/vs_pyhmmer.tsv`](../perf/vs_pyhmmer.tsv): the 12 scenarios of `RESULTS.md` (and the `--dbcache` variants, run by the fork only) and the 48 annotation runs of `annot.tsv`.
- The only change from the original workloads: every run writes `--tblout` and `--domtblout`, and `hmmsearch` sends its main output to `-o /dev/null`. `-A` is dropped from W1. `hmmsearch` still formats the report, alignments included, into `/dev/null`. That is a small cost that pyhmmer doesn't pay. Every contender computes the alignments, since the domain stage needs them.
- Threaded runs pin `-Z` or use `--cut_ga`, as before, so their outputs are deterministic.

**Threads.** `hmmsearch` runs with `--cpu 16`, as in `RESULTS.md`. pyhmmer's default (8) differs from that, so a calibration run came first. It ran pyhmmer with 8 and with 16 threads, interleaved, 3 rounds (`ab.sh --b-env PYHMMSEARCH_CPUS=8`):

| Workload | 16 threads | 8 threads | 8 / 16 |
|---|---:|---:|---:|
| X5 Pfam 1k × human | 14.5 s | 14.7 s | 1.01 (~) |
| X4 Pfam-A pressed × 5 | 2.23 s | 2.04 s | 0.91 (faster) |
| X4 Pfam-A text × 5 | 16.9 s | 16.8 s | 0.99 (~) |
| X3 Pfam 1k × Swiss-Prot 1k | 0.94 s | 0.89 s | 0.94 (~) |
| X1 Pfam 100 × Swiss-Prot 50k | 2.75 s | 3.11 s | 1.13 (slower) |
| X2 1 model × Swiss-Prot | 2.54 s | 2.83 s | 1.12 (slower) |
| bacteria L1 (pressed Pfam-A) | 2.32 s | 1.91 s | 0.82 (faster) |
| human L4 (pressed Pfam-A) | 9.43 s | 12.8 s | 1.36 (slower) |
| bacteria L1 (text Pfam-A) | 16.6 s | 16.8 s | 1.01 (~) |

Neither setting is best everywhere. So the main runs give pyhmmer both settings, as two contenders. Its charts and tables show, per workload, the faster of the two: `plot_*.py --best-of pyhmmer,pyhmmer-c8`, marked `*` in a table when 8 won. Picking the lower of two medians favors pyhmmer slightly. The text annotation runs keep 16 only: they are bound by reading the text HMMs, and the two settings were within 1%.

**Output check.**

- In the untimed first round, `series.py` compares each contender's files with upstream's after the usual normalization.
- `--mask-tail` drops the `# Program: ...` block from the tables, since pyhmmer doesn't write it.
- `--flush-subnormals` compares numbers with an exponent of -309 or below as 0. All such numbers are below DBL_MIN, 2.2e-308. The run's notes count them for each contender, so the scale of the tolerance stays visible. The `-ftz` variant must have none.
- A negative control checked that the masks don't hide real differences: an X5 table with one bit score changed by 0.1, an E-value of `3.2e-305` (normal), or one row removed still compares as different.
- Whitespace-only differences count as `ws` (column padding), as in `ab.sh`. Any other difference fails the run.
- The fork's `--dbcache` runs are checked against the fork's own default runs here; `RESULTS.md` already checked them against upstream.

**Timing.**

- `dev/perf/series.py`: 1 untimed round, then 3 timed rounds. Each round runs every contender once per workload, in an order rotated each round. The metric is the median wall time.
- No `perf stat` (it slows threaded runs on this VM; LEARNINGS measuring 20). GNU `time` records user, sys and peak RSS.
- All three contenders run in the same session, interleaved, so VM drift hits them equally. Old numbers from `RESULTS.md` are not mixed in.
- Three sessions (`bench-results/pyhmmer/run-all.sh`), run detached, because background tasks here stop after 2 hours:
  1. the 15 scenario workloads: upstream, pyhmmer, pyhmmer-c8, pyhmmer-ftz, fork;
  2. the 24 pressed annotation groups, with the same five;
  3. the 24 text annotation groups: upstream, pyhmmer, fork.

**Reporting** (`RESULTS.md`).

- The per-scenario speedup chart becomes speedup over upstream for pyhmmer and the fork (fork with `--dbcache` where it applies). A table of wall times for the three goes with it.
- The annotation charts get a third bar per length group. They stay one grouped bar chart per taxon, in sequences per second.
- The headline gets the pyhmmer numbers.
- The per-phase charts stay as they are: they measure the fork's history, not pyhmmer.

**Not in scope.**

- pyhmmer built from source with `-march=native`, or a pyhmmer patched with wider SIMD.
- pyhmmer's `hmmscan`, `phmmer` and `nhmmer`.
- Pressed optimized profiles (`HMMPressedFile`) as `hmmsearch` queries. That is a pyhmmer-only shortcut, and it skips the profile setup that all three otherwise do per query.
- Databases larger than memory: pyhmmer's recommended use needs the whole target set in RAM.

## Findings

### Outputs

- **The fork's tables matched upstream's in every run. pyhmmer's matched in all 12 scenarios and in 44 of the 48 annotation runs** (22 of 24 pressed, 22 of 24 text), for every pyhmmer contender.
- **The only differences come from subnormal numbers**, in two ways:
  - **Printed subnormals.** Stock pyhmmer printed numbers below 2.2e-308 where `hmmsearch` prints `0`: 3 on X5 (1 in `t.tbl`, 2 in `t.dom`), and 1-10 in each of 5 annotation groups, from the pressed and the text file alike (bacteria L4, phages L4, viruses L3 and L4, human L3). `--flush-subnormals` accepts these.
  - **E-values just above them.** In 2 annotation groups, from the pressed and the text file alike, one hit's E-values differ, and nothing else in their tables. Bacteria L4, NhaB in NHAB_ECOK1: `hmmsearch` prints `0`, pyhmmer `1.5e-307`. Human L4, Myosin_head in MYH9_HUMAN: `0` against `3.7e-306`. The cause is the same, one step earlier: the hit's P-value is below 2.2e-308, so `hmmsearch` flushes it to 0 before multiplying it by Z = 250 into the E-value. pyhmmer keeps the subnormal P-value, and the product is a normal number. A rerun of the two models alone (`hmmsearch` from this branch, same 250 targets) printed the same rows as upstream, `0` in those fields. After the sessions, upstream's tables for the two groups were made again and compared field by field with pyhmmer's 4 runs (pressed and text): in each table, that one row differs, and only in its E-value fields (2 in `t.tbl`, 3 in `t.dom`). Bit scores, bias, coordinates and every other row are the same.
  - **With flush-to-zero on** (`pyhmmer-ftz`), pyhmmer printed neither kind. Its tables were byte-identical to `hmmsearch`'s in all 36 runs it took part in, the `# Program:` block aside.
- **The column widths matched too.** pyhmmer keeps the same hits as `hmmsearch` even when it splits the targets among threads (X2).

### Scenarios (2026-10-20, 16 threads, median of 3)

| Scenario | upstream | pyhmmer | this branch | pyhmmer vs upstream | this branch vs upstream | this branch vs pyhmmer |
|---|---:|---:|---:|---:|---:|---:|
| X5 Pfam 1k × human proteome | 56.0 s | 14.5 s | 7.8 s | 3.9× | 7.2× | 1.9× |
| X4 All of Pfam-A × 5 sequences | 22.8 s | 2.1 s* | 1.4 s (1.1 s with `--dbcache`) | 11× | 16× (20×) | 1.5× (1.9×) |
| X4 the same from the text file | 38.4 s | 15.8 s | 6.8 s | 2.4× | 5.6× | 2.3× |
| X3 Pfam 1k × Swiss-Prot 1k | 5.9 s | 0.84 s* | 0.55 s (0.49 s) | 7.0× | 11× (12×) | 1.5× (1.7×) |
| X1 Pfam 100 × Swiss-Prot 50k | 8.5 s | 2.6 s | 1.2 s | 3.2× | 7.2× | 2.2× |
| X2 1 model × Swiss-Prot | 0.99 s | 2.5 s | 0.46 s | 0.39× | 2.1× | 5.4× |
| W2 40 models × Swiss-Prot 10k, 1 thread | 3.5 s | 3.3 s | 1.9 s | 1.0× | 1.8× | 1.7× |
| W3 1 model × Swiss-Prot, 1 thread | 3.9 s | 6.0 s | 2.3 s | 0.65× | 1.7× | 2.6× |
| R1 1 short model × Swiss-Prot, 1 thread | 1.2 s | 2.1 s | 0.57 s | 0.54× | 2.0× | 3.7× |
| W4 10 models × long sequences, 1 thread | 0.32 s | 0.30 s | 0.15 s | 1.1× | 2.1× | 2.0× |
| W1 10 hit-rich models × 10k, 1 thread | 0.83 s | 1.5 s | 0.57 s | 0.54× | 1.5× | 2.7× |
| W5 10 models × 1k, `--max`, 1 thread | 3.5 s | 4.1 s | 2.7 s | 0.86× | 1.3× | 1.5× |

\* pyhmmer at 8 threads, which was faster there. Everywhere else, 16 threads were faster or the run was single-threaded.

1. **pyhmmer's claim holds against upstream where there are many queries.** It is 2.4-11× faster than upstream on every threaded run with many models. The reason is its design: the targets stay in memory, and a thread takes a whole query, so threads never wait for each other between queries. Upstream re-reads the database for every query and stops all threads at the end of each one. Those are the same two costs that this branch's Phases 2, 3 and 7 removed.
2. **pyhmmer is slower than upstream in the other half of the scenarios:**
   - **One model against a big database:** X2 0.39×, W3 0.65×, R1 0.54×. pyhmmer first reads all 570k Swiss-Prot sequences into memory, serially. That takes 1.35-1.4 s, more than upstream's whole X2 run, and 670-713 MB of memory (this branch: 4-125 MB). Its `targets` mode then searches X2 in 0.86 s, while upstream needs 0.99 s for reading and searching together.
   - **Hit-rich and domain-heavy single-threaded runs:** W1 0.54×, W5 0.86×. These are subnormal numbers (next item).
   - **Otherwise even.** On single-threaded runs that are bound by the filters it matches upstream (W2 1.0×, W4 1.1×), as expected from the same SSE kernels.
3. **Subnormals cost pyhmmer up to 1.8×.**

   | Workload | pyhmmer | with flush-to-zero | gain |
   |---|---:|---:|---:|
   | W1 10 hit-rich × 10k | 1.55 s | 0.85 s | 1.82× |
   | X5 Pfam 1k × human | 14.5 s | 11.1 s | 1.30× |
   | W3 1 model × Swiss-Prot | 5.99 s | 4.82 s | 1.24× |
   | W2 40 models × 10k | 3.32 s | 2.90 s | 1.15× |
   | W5 10 × 1k `--max` | 4.09 s | 3.61 s | 1.13× |
   | X1, X2, X3 | | | 1.03-1.07× |
   | X4, R1, W4 (filter- or parse-bound) | | | ≈1.0× |

   The CPU counters confirm the mechanism (W1, `perf stat -e assists.fp,cycles:u,instructions:u`, after the sessions):

   | W1, 1 thread | FP microcode assists | cycles | instructions |
   |---|---:|---:|---:|
   | upstream | 0 | 2.85 G | 9.45 G |
   | pyhmmer | 21.1 M | 5.36 G | 8.43 G |
   | pyhmmer with flush-to-zero | 0 | 2.71 G | 8.43 G |
   | this branch | 0 | 1.90 G | 5.56 G |

   The same instructions take twice the cycles, about 125 cycles per assist: each operation on a subnormal operand or result goes through microcode. (pyhmmer runs fewer instructions than upstream because it reads the 10,000 targets once. Upstream parses them again for each of the 10 models, about 15% of its cycles here (`addbuf`, `header_fasta`, `seebuf` in `perf record`). Formatting the alignments into `/dev/null` costs it almost nothing: `--noali` saves 0.2% of the instructions.)

   With FTZ/DAZ on, pyhmmer's single-threaded runs match upstream's (W1 0.85 vs 0.83 s, W5 3.61 vs 3.53 s). pyhmmer could get this without changing Python's arithmetic. MXCSR is per thread, so it could set the two flags in its worker threads only, or around each `nogil` search call. The outputs would then also equal `hmmsearch`'s.
4. **This branch is faster than pyhmmer in every scenario: 1.5-5.4×, and 1.7-1.9× with `--dbcache`.** Most of that comes from changes that pyhmmer's design doesn't give it:
   - the AVX2/AVX-512 filters (Phase 6);
   - the paired Forward/Backward parsers (7C) and the domain-stage rewrites (Phase 4);
   - the faster FASTA reader (Phase 5, 5.6), which streams the targets while the workers search;
   - the exact `atof()` fast path for text HMMs (7B). Reading the 21,979 text models is pyhmmer's floor on X4 text: 15.8 s against 6.8 s.
   - FTZ itself.
5. **Thread count.** pyhmmer at 8 threads, its default here, was 4% faster than at 16 on X4 pressed and 11% faster on X3, the two shortest many-model runs. It was slower on the rest: X1 13%, X2 9%, X5 8%. On X5, 8 threads used 119 instead of 174 CPU-seconds for 8% more wall time: the hyperthreads add little to pyhmmer.
6. **CPU time on X5:** upstream 169 CPU-s, pyhmmer 174 (16 threads), this branch 114.

### Annotation benchmark (2026-10-20, 16 threads, median of 3)

All 21,979 Pfam-A models against each group of 250 proteins, `--cut_ga`. Sequences per second (the full tables per group are in `RESULTS.md`).

#### Pressed Pfam-A

| Set | upstream | pyhmmer | this branch | pyhmmer vs upstream | this branch vs upstream | this branch vs pyhmmer |
|---|---:|---:|---:|---:|---:|---:|
| Bacteria | 6.4 | 97 | 103 | 15× | 16× | 1.1× |
| Phages | 6.8 | 105 | 108 | 15× | 16× | 1.0× |
| Fungi | 5.0 | 76 | 74 | 15× | 15× | 0.98× |
| Eukaryotic viruses | 4.9 | 71 | 76 | 15× | 15× | 1.1× |
| Non-human eukaryotes | 4.8 | 69 | 61 | 14× | 13× | 0.88× |
| Human | 4.0 | 55 | 41 | 14× | 10× | 0.74× |
| **All 6,000, pressed Pfam-A** | **5.2** | **75** | **69** | **15×** | **13×** | **0.92×** |

pyhmmer is again the faster of its 16- and 8-thread runs per group. 8 threads won 14 of the 24 groups: every L1 and L2 group and 2 of the L3 ones.

1. **Pressed Pfam-A: pyhmmer is 8% faster than this branch overall** (75 against 69 sequences per second), and both are 13-15× faster than upstream.
   - This branch is ahead by 5-22% in 9 of the 24 groups, within ±5% in 8, and behind in 7.
   - The 7 are all eukaryotic, and the longest proteins lose most: human L4 (696-13,477 residues) 0.64×, eukaryotes L4 0.81×, human L3 0.83×, eukaryotes L2 0.88×, human L1 0.91×, fungi L4 0.92×, eukaryotes L1 0.95×.
2. **On long proteins, pyhmmer keeps more cores busy.** CPU time divided by wall time, on the L4 groups: pyhmmer 13-15 cores at 16 threads, this branch 5.7-9.6. On the L4 groups this branch needs no more CPU time than pyhmmer at 8 threads, and less on the eukaryotic ones. Human L4: 84 CPU-s, against pyhmmer's 101 at 8 threads and 145 at 16 (where a hyperthread's time counts like a core's). It loses that run on occupancy alone.
   - This is the limit found on 2026-10-09 (Phase 7 batches): results print in query order, and at most 2 batches of 16 queries are in flight. A few (query, target) pairs take 0.1-0.7 s in the domain stage, e.g. the 13,477-residue MUC3B_HUMAN against low-complexity families. While one runs, the window can't advance, and the other workers wait.
   - pyhmmer bounds only the queue of queries that no worker has taken yet (2 × threads, `_base.py`). It keeps finished results in a list until all earlier ones are done, so no worker waits for a slow query. The cost is memory: peak RSS for human L4 is 1,137 MB against this branch's 265 MB (upstream 192 MB). Over all 24 groups pyhmmer used 48-1,137 MB, this branch 71-340 MB, upstream 33-192 MB.
3. **On short proteins the two are close** (L1: 0.91-1.22×). There the run is mostly per-model work: reading and configuring 21,979 profiles for 250 short targets. Both spread it over the threads. They use about the same CPU time there (10-15 CPU-s per L1 group, pyhmmer at 8 threads), so this branch's faster filters don't show.
4. **pyhmmer's thread count matters here.** At 8 threads it caps at 8 cores, which is faster on the short groups (bacteria L1 1.95 s against 2.46 s at 16) and slower on the long ones (human L4 12.6 s against 9.5 s). At its default here, 8, it would lose to this branch on 5 of the 6 L4 groups, all but human.
5. **Flush-to-zero gains pyhmmer up to 1.19×** (pressed, 16 threads), the most on the longest groups: eukaryotes L4 1.19×, human L4 1.15×, viruses L4 1.13×. On L1 and L2 it is 0.96-1.05×.

#### Text Pfam-A

| Set | upstream | pyhmmer | this branch | pyhmmer vs upstream | this branch vs upstream | this branch vs pyhmmer |
|---|---:|---:|---:|---:|---:|---:|
| Bacteria | 4.5 | 15 | 32 | 3.3× | 7.0× | 2.1× |
| Phages | 4.7 | 15 | 32 | 3.2× | 6.7× | 2.1× |
| Fungi | 3.7 | 15 | 29 | 4.0× | 7.8× | 2.0× |
| Eukaryotic viruses | 3.7 | 15 | 29 | 4.0× | 7.9× | 2.0× |
| Non-human eukaryotes | 3.7 | 15 | 28 | 4.0× | 7.5× | 1.9× |
| Human | 3.2 | 14 | 23 | 4.5× | 7.2× | 1.6× |
| **All 6,000, text Pfam-A** | **3.9** | **15** | **28** | **3.8×** | **7.4×** | **1.9×** |

pyhmmer ran at 16 threads only here.

6. **From the text file, this branch is 1.9× faster than pyhmmer overall**, and 1.4-2.2× in every group but human L4 (1.05×).
7. **pyhmmer is flat at 13-15 sequences per second.** It reads the 21,979 text models with HMMER 3.4's parser in its main thread, while the workers search. Its runs take 16.6-19.3 s whatever the protein length, about the 15.8 s of the X4 text scenario, which has only 5 targets. Against the pressed runs, the text file adds 10-14 s to pyhmmer's runs, 15-18 s to upstream's and 3.8-5.8 s to this branch's (Phase 7's exact fast `atof()`).
8. **The two meet on human L4** (13 and 14 sequences per second): the search there takes about as long as pyhmmer's reading, and this branch's occupancy limit (item 2) applies again.
9. **Memory, text runs:** upstream 33-191 MB, pyhmmer 59-856 MB, this branch 73-347 MB.
