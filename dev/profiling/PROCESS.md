# How the hmmsearch speedup was done: process and lessons

This page describes the process behind the hmmsearch speedup work on this fork. It covers the constraints, the safety net, the profiling, the plan, the per-step loop, the reporting and how the work was kept re-runnable on other machines. It also says which practices decided the outcome and which didn't. It is meant as a playbook for performance work on a scientific tool whose outputs must not change.

Implementation details are deliberately left out. They are in the other pages:

| Page | Role |
|---|---|
| [`README.md`](README.md) | First profile (2026-10-04): where the time goes, what bounds wall time |
| [`PLAN.md`](PLAN.md) | The step plan: every step with its byte-identity argument, metric and outcome; the user's decisions |
| [`PROGRESS.md`](PROGRESS.md) | One measured row per step, including dropped steps |
| [`REPROFILE.md`](REPROFILE.md) | Second profile (2026-10-06), after the bottleneck moved |
| [`LEARNINGS.md`](LEARNINGS.md) | Numbered lessons, each with the number that taught it |
| [`RESULTS.md`](RESULTS.md) | The report: upstream vs every phase, charts, method, reproduce commands |
| [`GPU.md`](GPU.md), [`GPU-AGENT-PROMPT.md`](GPU-AGENT-PROMPT.md) | Ceiling test for a GPU path and the task prompt that encodes this process for a follow-on agent |
| [`../README.md`](../README.md), [`../bench/README.md`](../bench/README.md), [`../perf/README.md`](../perf/README.md) | The test suite, the benchmark, the A/B tooling |

---

## Outcome in brief

- **Duration and size.** The work ran from 2026-10-03 to 2026-10-09. It took 92 commits on `dev` and about 45 kept optimization steps in 7 phases. About 6,000 lines were added in `src/`, plus parser changes in a pinned fork of the Easel library.
- **Correctness.** Every kept step produces output byte-identical to upstream HMMER 3.4.
- **Speed at 16 threads** (from [`RESULTS.md`](RESULTS.md)):
  - realistic annotation (Pfam 1k × human proteome): 56 → 7.9 s (7.1×);
  - all of Pfam-A × 5 sequences: 24 → 1.5 s (16×);
  - Pfam-A annotation of 6,000 UniProt proteins: 5.2 → 71 sequences/s (14×).
- **Speed single-threaded:** 1.3-2.1×.
- **Memory.** Peak memory on long proteins grew in Phase 3 and was found only four days later. After a fix it is 220 MB on the realistic run (upstream 341 MB) and 272 MB on the longest human proteins (upstream 192 MB).
- **GPU.** A GPU path was evaluated by measuring its ceiling before writing any GPU code. The ceiling was 1.1-3.3×, so the user chose to test larger CPU machines instead and the path wasn't built.

## The process at a glance

| # | Stage | What it produced | Where | Verdict |
|---|---|---|---|---|
| 1 | Constraints and decisions | Byte-identical rule, opt-in rule, dependency policy, dated decision log | PLAN.md §6 | Crucial |
| 2 | Safety net | One command: build, unit tests, golden outputs of every program | `dev/check.sh`, `dev/regress/` | Crucial |
| 3 | Large-data benchmark | Pinned real datasets, 3 speed tiers, output hashes per case | `dev/bench/` | Hashes crucial; stored-baseline speed comparisons rarely decisive |
| 4 | Profile | Stage costs, the wall-time bound, ranked directions | README.md | Crucial |
| 5 | Plan | Small steps, each with an exactness argument, a metric and an estimate; accept/drop rules | PLAN.md | Crucial |
| 6 | Measurement protocol | Snapshots, interleaved A/B with output identity, noise calibration, a metric per bound | `dev/perf/` | Crucial |
| 7 | Step loop | ~45 steps kept, 6 dropped, evidence per step, phase milestones | PROGRESS.md | Crucial |
| 8 | Re-profile | The bound had moved; this led to the largest default-path gain | REPROFILE.md | Crucial |
| 9 | Use-case benchmark | Pfam-A annotation, 6 taxa × 4 length groups | `dev/bench/uniprot.py`, `dev/perf/annot.tsv` | Important; should have come first |
| 10 | Report | Upstream vs every phase in one session; charts; method; reproduce | RESULTS.md | Needed to communicate; not used for decisions |
| 11 | Portability | Pinned inputs and dependencies, machine info in results, run-time CPU dispatch, host-derived limits | across `dev/` | Important |
| 12 | Ceiling test before a big bet | GPU record/replay; stopped at the gate | GPU.md | Saved a large effort |
| 13 | Knowledge capture | Lessons, progress log, decision log, handoff prompt | LEARNINGS.md, GPU-AGENT-PROMPT.md | Important |

---

## 1. Constraints first

These were settled before any code changed.

- **Goal.** Make hmmsearch faster for annotation at scale: many profile HMMs (Pfam) against many proteins, eventually whole sequence databases.
- **Correctness means byte-identical.** Outputs must match upstream byte for byte after masking dates and timings. That covers scores, E-values, alignments, column widths and binary files, not just the same hits. "Is it correct?" then becomes a mechanical check with no judgement calls.
- **Compatibility.** Default command lines behave exactly as before. New behaviour, such as the in-memory target cache, is opt-in.
- **Dependencies.** Changes to the Easel library go to a fork the user owns, pinned by commit.
- **Decision log.** Some questions are the user's to answer. Each was asked once, answered and recorded with its date (PLAN.md §6). Examples:
  - the cache is opt-in, with a memory cap derived from the host;
  - whitespace-only output differences are acceptable (an upstream flake, §5.5);
  - no non-byte-identical speedups for now;
  - a proposed scheduling fix was declined, and that is recorded so it isn't proposed again;
  - the GPU work stopped at its first gate.

**Why it mattered.** The byte-identical rule ruled out a class of tempting gains: fewer stochastic samples, wider float SIMD, fast-math, and `-march=native` (which changed last digits, presumably through FMA contraction). In exchange, about 45 aggressive changes could be checked by machine instead of reviewed by eye. They included threading, a parser rewrite and new SIMD kernels.

---

## 2. Safety net before touching code

### 2.1 Golden-output regression

- **One command, two speeds.** `dev/check.sh` builds, runs the unit and integration tests (394 exercises) and compares the outputs of about 100 commands with stored goldens. The commands cover every user-facing program. The full run takes about 45 s and `--quick` about 5 s, so it runs after every edit.
- **Normalization.** Only volatile lines are masked: dates, CPU time and run time. Everything else must match byte for byte, binary files included.
- **Stable goldens.** Before the goldens were used, they were checked to be deterministic across runs, thread counts, and `-O0` vs `-O3` builds.
- **Rules for changing them:**
  - a golden update needs its own explanation in the commit, and an unexplained golden diff counts as a regression;
  - goldens for new coverage are generated with the unmodified binary and committed separately as "coverage only".

### 2.2 Benchmark on large real data, with output hashes

- **Pinned public inputs.** Pfam 37.0, Swiss-Prot, Dfam, Rfam, hg38 chr21 and E. coli are pinned by URL and md5.
- **Derived inputs made without the program under test.** Subsets, length groups, decoys and random sets are made by a script that uses only awk, gzip and the pinned Easel tools, never HMMER. They are checked by sha256 on every run.
- **Tiers.** 133 cases in three tiers: S (~9 min), M (~35 min) and L (~2.5 h).
- **Correctness on every case.** Each case's outputs are hash-checked, and repeated runs must be identical, which catches nondeterminism.
- **Speed.** Speed is compared with a per-machine baseline. This half turned out to be the least reliable part (§5.5). The hash half was used at every step.
- **A free first look.** The first full run already showed the main hotspots: thread scaling stopped at about 4 threads, hmmscan used fewer than 2 cores, and memory spiked on long queries.

### 2.3 Coverage aimed at what will change

Before each risky phase, goldens were added for the paths the plan would touch, generated by the old binary.

- **Before the threading changes:** 11 cases. They cover multi-query runs on several threads with several blocks per worker, `-Z` unset, pressed query files, the `.gz`/stdin error paths, and restricted database subsets.
- **After a scheduling change:** moving *when* output is printed changed what is printed before two fatal errors. No golden covered either. Both were fixed. An error-path comparator was then added and made a gate for every scheduling change. It runs two binaries on HMM-file errors, target-file errors and batching cases at several thread counts.

### 2.4 Diversity of inputs

Almost every optimization won on some inputs and lost on others. The workload set therefore spans these axes:

| Axis | Range covered | What it caught |
|---|---|---|
| Model length | under 60 to over 1,000 match states | Wide SIMD lost on short models; 512-bit vectors paid only on long ones |
| Target length | 60 aa or less to 35k aa | Per-sequence overhead; memory growth on long proteins |
| Hit density | random and decoy databases to hit-rich families, plus `--max` | Domain-stage cost; where pairing targets pays |
| Database size | 5 to 570k sequences (206 M residues) | A mode that wins on small databases ran 1.5× slower on big ones |
| Query count | 1 to 21,979 models | Per-query fixed costs: thread creation, model parsing, per-model setup |
| Input formats | FASTA, EMBL, gzip, stdin; text vs pressed HMMs | Text model parsing was 54% `atof()` |
| Threads | `--cpu 0`, 1 to 16 | The reader bound, the per-query barrier, hyperthreading effects |
| Options | E-value and bit-score thresholds, `-Z`, `--cut_ga`, `-A`, all tables | Determinism classes, error paths |
| Real use case | All of Pfam-A × 6 taxa × 4 length groups | A memory regression; thread occupancy on long proteins |

---

## 3. Profiling: find the bound, not just the hot function

### 3.1 First profile (analysis only, no source changes)

Method:

- **A profiling build.** Same optimization flags plus frame pointers. Its outputs were checked byte-identical and its speed within noise of the normal build.
- **Stage attribution.** Each sample was assigned to a pipeline stage by its call chain: the reader, each filter, each domain-stage step and output. This was done on about 16 workloads across the axes above.
- **Per-stage throughput.** Exact work counts (uprobe counts of cells per kernel call) give cells per second for each stage. Compared with the kernel microbenchmarks, they show whether each stage runs at its design speed.
- **Scaling.** A thread-scaling table from 0 to 16 threads, with CPU utilization per run.
- **Serial bottlenecks.** Busy time per thread shows which thread sets the pace.
- **One flags experiment.** `-march=native` gave no speedup and changed outputs.

Output:

- **A verdict per stage** (expected, inefficient or bottleneck), with its evidence.
- **A rule of thumb** for when the serial reader bounds the run: beyond about 3 cores for Pfam-like models. The scaling table confirmed it: 16 threads were slower than 8.
- **A projection for the target workload.** All of Pfam-A × Swiss-Prot had a reader floor of at least 6 h, whatever the core count.
- **Ranked directions**, each with its risk to byte-identity.

The key finding was structural. The hottest kernels already ran at their SIMD peak. Wall time at scale was set by serial work: the target database was re-parsed for every query, and every query paid for thread setup. That set the order of the plan.

### 3.2 Re-profile when the bound moves

After Phases 0-6, the same tools ran again, plus three additions:

- **Per-pass timestamps** from a throwaway patch. Workers were idle 39% of each query, waiting in its tail, and the reader was idle 40%, waiting for those tails. The new #1 bound was the per-query barrier. Per-thread CPU sums couldn't show it.
- **Cheap experiments** that tested each hypothesis before it entered the plan:
  - a stand-in for a hot library call gained 1%, not the 6.5% its profile share suggested (the samples were stalls);
  - a prototype fast path for number parsing in the text HMM reader was verified on 230 M real tokens;
  - LTO;
  - a forced execution mode on big inputs, which ran 1.5× slower;
  - two copies on hyperthread siblings, as a probe for latency-bound code.
- **Estimates that held.** This pass predicted ~7-8 s for the realistic run, and 7.85 s was measured.

Lessons:

- Profile shares mislead: a stall is charged to whatever instruction is waiting.
- Measure a ceiling with a stand-in or a replay instead of inferring it.
- Measure barriers with timestamps, not CPU sums.

---

## 4. Planning: from hypotheses to small steps

Each step in PLAN.md states four things, written before coding:

- **The change**, as one commit.
- **Why it's byte-identical.** Typical arguments:
  - the same operations in the same order;
  - integer DP whose result doesn't depend on the lane count;
  - scheduling only, which is safe because the RNG is reseeded per region and hits are sorted deterministically.
- **The workload set and metric** that judge it.
- **An expected effect**, estimated from the profile's stage shares.

**Rules written before Phase 1:**

- A step is kept only if it is byte-identical on every check and faster beyond noise on its target metric, with no other workload slower and any RSS growth explained.
- A neutral step stays only if a later step needs it, as a declared refactor.
- At the end of each phase: run everything on the default build, re-save the baseline and re-profile.

**Order by the bound.** Serial bounds go first, because CPU savings don't show at full thread count while a serial stage sets the pace. Measured in Phase 1: 2.7% fewer instructions showed as 1.6% wall time at 16 threads, within noise.

**Invariants first.** The invariants that make reordering safe were checked in the code and written down before any threading change.

**The plan was a living document:**

- About a dozen steps were added mid-phase, when per-step measurement showed the real limiter. After persistent threads came a kernel wake-up storm, then a single preparation helper on the critical path, then latency with no thread saturated.
- Six planned steps were measured and dropped.
- Estimates were replaced with milestone measurements as they came in.

**Hindsight.** The first profile offered two fixes for its #1 finding: keep the parsed database in memory, or search a batch of queries per database pass. The opt-in rule fitted the cache, so it went first. But it helped only runs that pass the flag: the default realistic run stayed at 50-56 s through Phase 4. Batching, rediscovered by the re-profile, gave the largest default-path gain (Phase 7, 2.7×). Rank fixes by what the default path gets.

---

## 5. Measurement protocol

### 5.1 Snapshots

- `snap.sh` freezes the binaries of any commit. It builds the commit in a temporary worktree, with the dependency commit that revision pins, and records the sha, the flags and the binary hash.
- Static linking makes a copied binary a complete "before".

### 5.2 Interleaved A/B with output identity

- **Interleaving.** A and B alternate (AB, BA, AB, ...), so machine drift hits both. The first pair is untimed and serves as the output reference.
- **Output identity.** Every output file of every run is normalized and hashed.
  - A must equal B.
  - Every rep must equal its first; otherwise the run is reported as NONDETERMINISTIC.
  - Whitespace-only differences are classified separately (the accepted upstream flake).
- **Verdict.** "Faster" or "SLOWER" only beyond a floor per metric and beyond both sides' spread; otherwise "~".
- **Records.** Every sample goes to a JSON file.
- **One binary, two paths.** A run-time switch selects the old or the new kernel, so A and B can be the same binary with no code-layout difference.

### 5.3 Calibrate the noise first (A/A)

The same binary on both sides gives each metric's noise on this machine:

| Metric | Spread | Use it for |
|---|---|---|
| Instructions, `--cpu 0` | 0.00% (deterministic) | Proving work was added or removed |
| Cycles, `--cpu 0` | 0.2-0.7% | Single-thread CPU steps |
| User + kernel cycles | 0.9% | Steps that move work into system calls; the VM's clock rate varied 8-14% between runs |
| Wall time, threaded, interleaved | 1.4-12% | Threading and scheduling steps; needs big effects or more reps |
| Wall time vs a stored baseline | drifts up to ~10% a day | Nothing |

### 5.4 Pick the metric that can resolve the effect

- **Single-thread CPU steps:** instructions and cycles at `--cpu 0`.
- **Threading and scheduling steps:** interleaved wall time.
- **Never:** wall time against a results file from another day.
- **Know the bound before predicting a wall-time gain.** Sum busy time per thread, or timestamp each pass.

### 5.5 What faked results, and how each was neutralized

| Artifact | Size | Fix |
|---|---|---|
| Code layout: editing one file shifted the alignment of every SIMD kernel | ±1-3% cycles with identical instructions | Build both sides with aligned functions. Build a control (B's code with the new behaviour switched off) when in doubt. Read cycles against instruction counts |
| Stored baselines | up to 10% a day; untouched kernels looked "12% faster" | Interleaved runs only |
| The harness's own timer polled the child process | wall times rounded up to the next 50 ms | Blocking wait |
| `perf stat` on threaded runs | +18% to +41% wall time, more on the side that synchronizes more | Count only the workloads judged by counters. Published speedups were corrected (one went from 20× to 16×) |
| VM speed drift | ~12% slower for hours | Compare only within one interleaved session |
| Short runs (under 0.5 s) | +13.7% by wall time where counters showed +4.6% | Counters, or many reps |
| Stale incremental builds (Makefiles didn't track a header) | a snapshot silently ran the old code | Negative controls; check the binary for a new symbol |
| Upstream nondeterminism without `-Z` | column widths varied in about 1 of 6 threaded runs | Pin `-Z` in threaded workloads; classify the difference as `ws`; count its frequency before and after scheduling changes |

### 5.6 Workload sets keyed to what a step targets

| Set | What it stresses | Metric |
|---|---|---|
| `cpu` | Filters and domain stage, single thread | cycles |
| `reader` | Sequence parsing | user + kernel cycles |
| `scale` | Thread scaling, reader bound | wall |
| `query` | Many models × small databases (per-query costs) | wall |
| `milestone` | Realistic annotation run | wall |
| `annot`, `annot_txt` | Use case: all of Pfam-A × 6 taxa × 4 length groups | wall, reported as sequences/s |

Each step runs the sets it targets. A milestone runs all of them.

---

## 6. The step loop

For each step:

1. Branch from `dev` and snapshot A.
2. Implement the step.
3. Run `dev/check.sh`: build, tests and goldens.
4. Check the real-data hashes (benchmark tier S, hmmsearch cases).
5. Add the gates for the kind of change:
   - kernel rewrites: a crosscheck build, run over the goldens and all workloads;
   - opt-in features: all goldens rerun with the option forced on, and identical instruction counts on the default path;
   - threading changes: ThreadSanitizer, AddressSanitizer, a stress loop of 100+ runs, and the error-path comparator.
6. A/B on the target sets, with both sides built with the same flags.
7. Record a row in PROGRESS.md, also for a dropped step. The commit message carries the exactness argument and the numbers.
8. At the end of a phase: all sets on the default build, benchmark tier M hashes, a re-profile, then fast-forward `dev`.

### 6.1 Layers of evidence for correctness

| Layer | What it proves | Used for |
|---|---|---|
| Golden regression | Every program, small inputs, byte for byte | Every edit |
| Benchmark output hashes | Large real inputs | Every step |
| A-vs-B output identity | Big workloads, every rep | Every step |
| Forced-option golden rerun | An opt-in feature gives the same outputs on every case | Opt-in features |
| Golden rerun with each code path forced | Every SIMD level gives the same outputs | Run-time-dispatched kernels |
| Bitwise unit test against the **original** routine | Kernel equivalence on random and edge inputs | Kernel rewrites |
| Mutation tests | The unit test can fail on a near-miss (a 6-ULP difference, an off-by-one overflow test) | Kernel rewrites |
| Crosscheck build | Old and new code run side by side on real data and abort on any bit difference. It caught boundary mutants that random tests missed | Kernel rewrites |
| Differential test against a verbatim copy of the old implementation, with read sizes down to 1 byte | Parser equivalence at every buffer boundary. It found an upstream buffer overflow and several boundary-dependent upstream behaviours | The parser rewrite |
| Sanitizers and stress loops | Races, leaks, hangs | Threading changes |
| Error-path comparator | Output before fatal errors | Scheduling changes |
| Negative control for every checker | The checker can fail | Every new checker |

Negative controls paid off twice:

- a crosscheck build turned out to have the check compiled out of the kernels' directory;
- a snapshot turned out to contain the old code.

### 6.2 Keep or drop

- About 45 steps were kept.
- Six planned steps were measured and dropped or deferred. Each gained under ~1% or lost time, and the numbers are recorded.
- More variants were rejected in experiments: 512-bit vectors for every model (slower), 64 KB-1 MB read buffers (slower than 16 KB), a first version of finer work splitting (more futex churn, and the upstream width flake in 9 of 10 runs), an exact inline `log()` (1%), and LTO (0.5-1%).
- Count how often a shortcut fires before adding it. One exact early-out almost never fired on real data and cost 2%.

### 6.3 Milestones catch what step A/Bs miss

- The Phase 6 milestone ran all 21 workloads and found two regressions that the step A/Bs, which ran only each step's target sets, had missed: short models, and per-query setup on many-model runs.
- The use-case benchmark found the Phase 3 memory growth four days after it was introduced: 754 MB peak against upstream's 192 MB.

---

## 7. Reporting

- **One results page** ([`RESULTS.md`](RESULTS.md)). The user wanted a single report, with new benchmarks added as sections.
- **One session for all binaries.** Upstream and the end of every phase ran on the same workloads in one session, round-robin in rotated order. An untimed first round checked every binary's outputs against upstream's.
- **Content:**
  - headline numbers;
  - a speedup-per-scenario chart and per-phase small multiples, light and dark SVG, drawn by a stdlib-only script (no matplotlib on the machine);
  - a "what changed / where it shows" table per phase;
  - the full table of medians and the peak memory;
  - the use-case section;
  - the method: machine, commit SHAs, rounds, spread;
  - reproduce commands and the raw JSON.
- **Metrics that fit the reader.** Annotation results are throughput (sequences/s), shown as one grouped bar chart per taxon with length groups as the groups. The user chose this format.
- **Corrections in the open.** When `perf stat` was found to inflate threaded speedups, every phase was re-measured. The page states the old numbers and why they changed.
- **Honest labels.** "Outputs identical" is claimed only for what was actually compared, and estimates are labelled as estimates.

---

## 8. Portability: re-runnable on another machine

Everything was measured on one 16-vCPU KVM guest: a Xeon Platinum 8468 (Sapphire Rapids) with 8 cores × 2 hyperthreads, 196 GB RAM and gcc 13.3. The aim was that the suite, the benchmarks and the report can be re-run elsewhere, a larger VM in particular.

| Concern | How it's handled |
|---|---|
| Inputs | Public sources pinned by md5. Derived sets are regenerated by a script and checked by sha256. The UniProt sets come with a manifest (accession, entry and sequence versions); entries that changed since the draw are fetched from UniSave at the recorded version, and the files are checked by sha256 |
| Dependencies | Easel is pinned by commit in `dev/easel.ref`, and `check.sh` clones it. Snapshots of old revisions build the Easel commit that revision pins |
| Binaries | HMMER and Easel are linked statically (only libc and libm are shared). Each snapshot has an `info.txt`, so binaries can be copied between machines |
| Results files | Record host, CPU model, core count, compiler and flags. Benchmark comparisons print a note when these differ. Baselines are per machine and not committed |
| Core counts | Workloads use `$NCPU`, and output hashes don't depend on the core count |
| CPU features | The SIMD level is chosen at run time from CPU feature detection, with an environment override. The goldens are rerun with each level forced, so a machine without AVX-512 or AVX2 is covered for correctness |
| Memory | The cache's default cap is derived from physical RAM or a lower cgroup limit |
| Goldens on other platforms | Deterministic across thread counts and `-O` levels on x86-64 with gcc 13. Another architecture or compiler may legitimately differ in a last digit; regenerate the goldens there (documented in `dev/README.md`) |
| Machine-tuned choices | The 512-bit threshold by model length, the read buffer size, the batch size and the work-slice count were tuned on Sapphire Rapids. Re-check them on other CPUs, such as AMD, where 512-bit execution differs |
| Scale-out | Only 16 vCPU were measured. Query sharding across processes was measured before Phase 7: 4 processes × 4 threads beat 1 × 16 (9.6 vs 22.1 s uncached), and the concatenated outputs were byte-identical. Phase 7 then brought 1 × 16 to 7.9 s. Scaling past 16 threads is unmeasured; GPU.md has a suggested test for 32-128 vCPU and a labelled GenBank-scale estimate |

**Re-running on a new machine:**

1. `dev/check.sh` clones the pinned Easel, builds, tests and runs the goldens. If goldens differ on another architecture or compiler, find out why before regenerating them.
2. Fetch the data: `dev/bench.sh prepare` (~2 GB) and `dev/bench/uniprot.py fetch`.
3. Record the machine: cores vs hyperthreads, AVX2/AVX-512 support, the cgroup memory limit, and which SIMD level the programs choose.
4. Calibrate noise with an A/A run (the same snapshot on both sides of `dev/perf/ab.sh`). The floors in §5.3 belong to the old VM.
5. Save a benchmark baseline for this machine: `dev/check.sh --bench M --bench-args "--save-baseline"`.
6. Snapshot upstream (`f8bebd8`) and current `dev` with `dev/perf/snap.sh --rev`. Run the series commands at the end of RESULTS.md, without `perf stat` on threaded runs.
7. On a larger VM, also measure thread scaling (1 × N) and process sharding (k processes × N/k threads) on the case suggested in GPU.md, and compare with the 16-vCPU numbers.
8. Re-check the machine-tuned choices above.

---

## 9. Bigger bets: measure the ceiling first

The GPU question followed the same pattern as the small experiments, at a larger scale.

1. **Ceiling before code.** A throwaway patch recorded every filter outcome in a single-threaded run, then replayed the run at 16 threads with the filter stages skipped. The outputs stayed identical, so the replay's time is the floor for a GPU that does the filters for free.
2. **Result.**
   - Realistic run: 1.76-2.03×.
   - Small annotation sets: 1.1-1.25×.
   - Pfam-A × 20-50k proteins: 1.6-3.3×.
   - The profile-share estimate had been 4.8×.
   - The domain stage, per-model overhead and the reader stay on the CPU.
3. **Gate.** The task prompt said to stop and ask the user if the ceilings came out low. The user decided to test larger CPU machines instead.

The cost was a few hours and no product code. The same idea, a stand-in or a replay before the real implementation, decided several smaller steps too (§3.2).

---

## 10. Knowledge capture and working model

- **The documents** (table at the top) were written as the work went:
  - PROGRESS.md gets a row per step, including dropped steps;
  - LEARNINGS.md gets each lesson with the number that taught it;
  - PLAN.md records every decision with its date.
- **Roles.**
  - The user set the goals and constraints, made the calls at decision points, reviewed, and did all pushes (to this repository and the Easel fork).
  - The agent (Claude Code) profiled, planned, implemented, measured and documented. It committed small steps on branches and fast-forwarded `dev` after each milestone.
- **State across sessions.** Short persistent notes carried decisions, gotchas and current status from one session to the next, alongside the repository's own documents.
- **Handoff.** [`GPU-AGENT-PROMPT.md`](GPU-AGENT-PROMPT.md) packs the whole process into one task prompt for the next piece of work: ground rules, reading order, the known state, decision gates, quality gates, the performance protocol and reporting rules. Its G0 gate stopped the GPU work as intended.

---

## 11. What was crucial, what helped, what didn't

### Crucial

| Practice | Why it mattered | Evidence |
|---|---|---|
| Byte-identical rule plus a one-command golden gate, built before any change | Turned correctness into a mechanical check; aggressive changes became safe | ~45 kept steps, including threading, a parser rewrite and new SIMD kernels, with byte-identical output |
| Profile for the bound, and order the work by it | CPU savings are invisible while a serial stage sets the pace | Phase 1: 2.7% fewer instructions, 1.6% wall time at 16 threads |
| Interleaved A/B with output identity, calibrated noise and a metric per bound | Without it, layout effects and drift produce false verdicts | Layout ±1-3%; stored baseline showed untouched code "12% faster" |
| Re-profile when the bound moves | The bound moved from the reader to the per-query barrier | Phase 7: 2.7× on the default realistic run, the largest default-path gain |
| Small steps with a written exactness argument and a keep/drop rule; negatives recorded | Each gain is attributable; nothing below noise accumulates | 6 planned steps dropped with numbers; about a dozen steps added from measurements |
| Layered evidence for rewrites, with a negative control for every checker | Random tests miss boundaries; checkers can be silently broken | Mutants caught only by crosscheck or scans; a crosscheck build that checked nothing |
| Full-set milestones and a real use-case benchmark | Step A/Bs see only their target sets | 2 regressions caught at the Phase 6 milestone; a 4× memory growth found by the annotation benchmark |
| Decision gates with the user, decisions logged | Scope questions answered once; nothing re-litigated | Opt-in cache, whitespace rule, declined fix, GPU stop |
| Ceiling test before a big investment | Profile-share estimates overstate gains | GPU ceiling 1.8-2× measured vs 4.8× estimated |

### Helpful, not decisive

- **Benchmark speed tiers against stored baselines.** The hash half was used at almost every step: 36 runs in 7 days, tier S per step, tier M at milestones. The speed half rarely decided anything because of drift. Tier L (~2.5 h) ran once, on day 1, to record its hashes.
- **Kernel microbenchmarks and exact cell counts.** They confirmed which stages already ran at their design speed. This was a one-off, but it shaped the plan.
- **Charts.** Needed to communicate the results; no decision depended on them.
- **The process-sharding measurement.** It informs how to use bigger machines, but isn't product work.

### Low value or harmful

- **Micro-steps worth under ~1%.** Five of the six dropped steps were this kind. Under layout noise they can't be resolved. Count how often the code fires first.
- **Compiler and build flags.** `-march=native` changed outputs with no gain. LTO gave 0.5-1%. PGO wasn't tried.
- **"Wider or bigger is faster" assumptions.** 512-bit vectors for every model were slower (fewer ports, lower clock at 16 threads), and 64 KB-1 MB read buffers were slower than 16 KB.
- **Gains estimated from profile shares.** A 6.5% libm share was worth 1%. A 4.8× GPU estimate measured 1.8-2×.
- **`perf stat` on wall-time runs, and comparisons against stored baselines.** Both produced wrong numbers that had to be corrected.

---

## 12. What to do differently next time

1. **Build the use-case benchmark on day 1.** The Pfam-A annotation runs came at the end. They have a different bound (per-model overhead, thread occupancy on long proteins), and they found the memory growth four days after it was introduced.
2. **Make peak memory an acceptance metric** at every milestone, on long-target workloads in particular, not just "reported".
3. **Validate the harness itself in Phase 0.** Check the timer's resolution and the profiler's overhead with runs of known length. Both problems were found late, and one forced a correction of published numbers.
4. **Rank fixes for the top bound by what the default path gets.** The batching fix was on the table from the first profile.
5. **Add error-path goldens before any scheduling change**, not after one breaks them.
6. **Use more than one machine.** All numbers come from one 16-vCPU Sapphire Rapids VM. Thresholds tuned there, and scaling past 16 threads, are unverified.

**Open decisions:**

- **Non-byte-identical speedups.** Wider float SIMD for Forward/Backward is the largest lever left: 2-4× on 16-40% of CPU. It needs a reviewed golden update, which is the user's call.
- **In-process gzip decompression.** It adds a dependency and was not started.

---

## Appendix A: Timeline

| Date | Work |
|---|---|
| 10-03 | Test suite (golden regression, one-command gate); benchmark with output hashes |
| 10-04 | First profile; coverage goldens for the planned changes; A/B tooling; the plan; Phase 1 |
| 10-05 | Phase 1 milestone; Phase 2 (target cache); LEARNINGS.md; Phases 3 and 4; Easel fork; start of Phase 5 |
| 10-06 | Phase 5 milestone; Phase 6; a late parser step (5.6); first RESULTS.md; re-profile; start of Phase 7 |
| 10-07 | Phase 7 and its milestone; RESULTS.md through Phase 7 |
| 10-09 | Annotation benchmark; RESULTS.md re-measured without `perf stat`; memory fix; GPU ceiling test, stopped at the gate |

## Appendix B: The phases in one line each

| Phase | Hypothesis from the profile | Where it showed (16 threads unless noted) |
|---|---|---|
| 0 | Tooling first: snapshots, A/B, coverage, workload sets | — |
| 1 | The bias filter is inefficient (per-call allocation, matrix) | ~4% where it runs often; protocol warm-up |
| 2 | The target database is re-parsed per query (opt-in cache) | Realistic run with the option 2.9× |
| 3 | Per-query fixed costs (threads, setup, model preparation) | Many models × small database 9× (19× with the cache) |
| 4 | Domain-stage hot spots, under an exact-RNG-stream constraint | Hit-rich and `--max` runs 9-14% |
| 5 | The FASTA parser is slow (Easel) | Parsing-heavy runs −45 to −50% wall |
| 6 | Integer filters are at the 128-bit peak; use wider SIMD | Filter-heavy runs 16-29% |
| 7 | (Re-profile) The per-query barrier, text model parsing, latency-bound Forward/Backward | Many-query threaded runs 2.4-4.5×; single thread 5-15% |
| 7E.1 | (Use-case benchmark) Memory kept by reused workers | Peak memory 754 → 272 MB on the longest proteins |
