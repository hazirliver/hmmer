# HMMER speed benchmark

The benchmark times a HMMER build on large real-world and synthetic datasets. Every case also checks correctness: its outputs must hash-match `golden.json`, and repeated runs must be identical.

```
dev/check.sh --bench          # full dev suite + benchmark tier S (~9 min), compared to your baseline
dev/bench.sh run S            # benchmark only (build-dev must already be built)
```

## First-time setup

```
dev/bench.sh prepare          # ~2 GB download (~5 min) + ~3.5 min to derive 3.1 GB of datasets
dev/check.sh --bench --bench-args "--save-baseline"   # record this machine's baseline
```

- **Downloads** go to `~/.cache/hmmer-bench/` (override with `HMMER_BENCH_DATA`). They're pinned by md5 in `sources.tsv`.
- **Derived datasets** are made by `prepare.sh` using only awk, gzip and the pinned Easel tools, never HMMER. HMMER changes therefore can't alter the inputs. They're bit-reproducible and checked against `data.sha256` on every run.

## Tiers

| Tier | Approx. time | Use |
|---|---|---|
| **S** | ~8.5 min, 92 cases | Every feature branch, before merging to dev |
| **M** | ~35 min, 121 cases (includes S) | Performance-oriented changes; adds 50k/full-Pfam/full-Swiss-Prot I/O scale and thread-scaling series |
| **L** | ~2.5 h, 133 cases (includes M) | Before releases; adds full Swiss-Prot searches, all of Pfam-A.seed and human Dfam vs chr21 |

The times above were measured on a 16-vCPU Xeon 8468 VM. Each case runs 3 reps, except that cases over 60 s run once.

## Datasets

**Protein sequences**

| Dataset | Source | Contents |
|---|---|---|
| `sprot*.fa` | Swiss-Prot (Pfam 37.0 snapshot) | Full set: 570k sequences, 206M aa. Deterministic 50k, 10k and 1k subsets |
| `sprot_10k.dat` | Swiss-Prot | The same 10k sequences in UniProt/EMBL text format |
| `sprot_10k.fa.gz`, `sprot_50k.fa.gz` | Swiss-Prot | gzip-compressed copies |
| `sprot_long.fa` | Swiss-Prot | Sequences of 3,000 aa or more, up to 35k |
| `sprot_short.fa` | Swiss-Prot | Sequences of 60 aa or fewer |
| `human.fa`, `q_human_{10,100,1000}.fa` | Swiss-Prot | Human proteome (20k sequences) and query sets drawn from it |
| `q_*.fa` | Swiss-Prot | Hand-picked queries: INS (tiny), HBA, P53, EGFR (kinase, many hits), BRCA1 and TITIN (34k aa) |
| `random_100k.fa` | Generated | 100k iid random sequences, where no hit survives the filters |
| `sprot_10k_dishuffled.fa` | Generated | Composition-preserving decoys |

**Profile HMMs**

| Dataset | Source | Contents |
|---|---|---|
| `Pfam-A.hmm` | Pfam 37.0 | All 21,979 models, plus deterministic 1k, 100 and 10 subsets |
| `pfam_len_{short,med,long,xlong,mix}` | Pfam 37.0 | Models grouped by length, from under 60 to over 1,000 match states |
| `pfam_hitters` | Pfam 37.0 | Pkinase, 7tm_1, WD40, zf-C2H2 and other families with very many hits |
| `kern_M{50,200,400,1000}` | Pfam 37.0 | Single models of fixed length, for the kernel microbenchmarks |
| `dfam_human.hmm` | Dfam 3.7 curated | 1,365 human repeat families (median length 572, max 16k), plus 100 and 15 subsets and AluY and L1HS single models |

**Alignments**

| Dataset | Source | Contents |
|---|---|---|
| `Pfam-A.seed.sto` | Pfam 37.0 | All seed alignments, plus 1k, 100 and largest-20 subsets |
| `msa_LRR_1.*` | Pfam 37.0 | One large alignment in each of 7 MSA formats: sto, afa, a2m, psiblast, clustal, phylip and selex |
| `rfam_*.sto` | Rfam 14.10 | RNA seed alignments: 100 by stride, plus tRNA, 5S, SSU, RNaseP, tmRNA and 6S |

**DNA**

| Dataset | Source | Contents |
|---|---|---|
| `chr21.fa` | hg38 | Human chromosome 21 (46.7 Mb), plus 5 Mb and 1 Mb slices and a gzip copy |
| `q_chr21_*` | hg38 | 100 kb chunks for nhmmscan, plus 2 kb and 6 kb sequence queries |
| `ecoli.fa` | E. coli K-12 genome | 4.6 Mb, a different base composition |
| `random_dna_10M.fa` | Generated | Random DNA |

**UniProt protein sets (Pfam-A annotation benchmark)**

Six sets of 1,000 UniProtKB/Swiss-Prot proteins (reviewed, no fragments), each split by length into four groups of 250, `L1` (shortest) to `L4` (longest). The bounds are each set's own quartiles. They're the targets of the annotation workloads in [`../perf/annot.tsv`](../perf/annot.tsv).

```
dev/bench/uniprot.py fetch       # ~10 s: 6,000 proteins by accession -> ~/.cache/hmmer-bench/uniprot-<hash>/
```

| Set | UniProtKB query | L1 | L2 | L3 | L4 |
|---|---|---|---|---|---|
| `bacteria` | `taxonomy_id:2` | 29-165 | 166-277 | 278-407 | 411-1,457 |
| `phage` | `taxonomy_id:10239 AND virus_host_id:2` (viruses of bacteria) | 23-93 | 93-164 | 164-336 | 336-3,500 |
| `fungi` | `taxonomy_id:4751` | 15-256 | 256-400 | 400-593 | 595-2,904 |
| `virus` | `taxonomy_id:10239 AND virus_host_id:2759 AND NOT virus_host_id:2` (viruses of eukaryotes) | 31-179 | 179-318 | 318-570 | 570-7,182 |
| `eukaryota` | `taxonomy_id:2759 AND NOT organism_id:9606` (fungi included) | 6-200 | 202-355 | 355-537 | 537-5,801 |
| `human` | `organism_id:9606` | 16-246 | 248-395 | 396-695 | 696-13,477 |

Files: `<set>_L<1-4>.fa` (250 proteins each) and `<set>.fa` (all 1,000), with headers `>sp|ACCESSION|ENTRY_NAME OX=taxid`.

- **Drawing the sets** (`uniprot.py select`, done once, UniProt release 2026_03). Every query also gets `AND reviewed:true AND fragment:false`. The script lists every matching accession, ranks them by sha256(`<set>:<accession>`) and keeps the first 1,000 not drawn by an earlier set (in the table's order), so no protein is in two sets. The 1,000 are then sorted by length (ties by accession) and cut into four groups. Within a group, proteins keep the draw order, which is effectively random.
- **Why Swiss-Prot.** Every reviewed entry can be listed, so the draw is uniform over the whole query. Host annotations, which define the phage and virus sets, are near complete there: 99% of reviewed viral entries have one, against 18% of unreviewed ones. Fragments are excluded by the query.
- **Reproducible in later releases.** [`uniprot/manifest.tsv`](uniprot/manifest.tsv) records every protein's set, group, accession, entry name, entry and sequence versions, length and taxon. `fetch` downloads them from the UniProtKB REST API in batches of 250. A protein whose sequence version has changed since 2026_03, or that was deleted or merged, comes from UniSave at the manifest's entry version instead. Headers are built from the manifest, so annotation updates don't change the files. The files are checked against [`uniprot/data.sha256`](uniprot/data.sha256).

## Cases

`cases.tsv` covers every stage of the toolchain:

- **build:** configure, `make -j`, `make dev`.
- **kernel:** each pipeline stage is timed in isolation with HMMER's own benchmark drivers. The stages are SSV/MSV, Viterbi filter, Forward/Backward, posterior decoding, optimal-accuracy alignment, stochastic traceback and null2. Each runs at model lengths of 50 to 1,000 and is timed with both the SIMD and generic implementations. Times come from the driver's internal timer, which excludes setup, and Mc/s is recorded.
- **hmmbuild:** multi-MSA, big MSAs, RNA, every input format, and weighting/effective-number variants.
- **hmmsearch:** model-length and sequence-length extremes, hit-heavy families, null and decoy databases, and input formats (FASTA, EMBL, gzip, stdin). Also `--max`, `--nobias`/`--nonull2` and `--cut_ga`, plus thread scaling at 0/1/2/4/8/N threads.
- **hmmscan:** against the Pfam 1k subset and full Pfam, including titin.
- **phmmer and jackhmmer:** iterative search with checkpoints.
- **hmmalign:** large sequence sets and output formats.
- **nhmmer:** repeat families on chr21, a random genome, RNA models on E. coli, the FM index built by makehmmerdb, sequence and MSA queries, and strand and gzip variants.
- **nhmmscan, hmmpress, hmmstat, hmmconvert, hmmfetch, hmmemit, hmmsim, hmmlogo and alimask.**

## What gets measured and checked

For each case the harness records:

- wall time (median of `--reps`, default 3; a case slower than 60 s on its first run is run once)
- min wall time
- user and sys CPU
- peak RSS
- Mc/s, for kernel cases

Results go to `bench-results/<time>_<sha>_<tier>[_label].json`, which is gitignored.

**Correctness:** a case fails when any of these happens:

- its normalized outputs differ from `golden.json`
- two reps produce different outputs (`NONDETERMINISTIC`)
- it exits nonzero

Normalization is the same as `dev/regress`: only dates and timings are masked.

**Speed:** results are compared with `bench-results/baseline-<tier>.json`, or with whatever `--baseline` names.

- **SLOWER:** the best new run is more than `--tolerance` (default 10%) slower than the baseline median.
- **faster:** the new median is that much faster.
- **Ignored:** cases under `--min-time` (0.5 s). They're still checked for correctness.
- **Summary:** the comparison also prints the geometric-mean ratio per group and overall.
- **Exit code:** a slowdown doesn't make the run exit nonzero unless you pass `--fail-on-slowdown`. Timing noise on shared VMs is real, so re-run before you trust a single SLOWER.

## Typical workflows

```
# Before starting a feature: record the baseline on dev
git checkout dev && dev/check.sh --bench M --bench-args "--save-baseline --label dev"

# On the feature branch: full correctness + speed vs that baseline
dev/check.sh --bench M

# A/B two builds explicitly, for example dev vs a feature built elsewhere
HMMER_BUILD_DIR=$PWD/build-feat dev/check.sh --quick           # builds into build-feat
dev/bench.sh run M --build-dir build-feat --label feat --baseline baseline-M.json

# Compare any two result files later
dev/bench.sh compare bench-results/A.json bench-results/B.json

# Only some cases or groups (regex on name or group)
dev/bench.sh run M --cases 'hmmsearch|kernel'

# Intentional output change: refresh the hashes and commit golden.json
dev/bench.sh run L --reps 1 --update-hashes --baseline none
```

**Tips for low-noise timing:**

- Keep the machine idle. The harness warns if the load average is above 1.
- Use the same build flags as the baseline. The comparison prints a note when host, CPU, compiler or flags differ.
- Use `--reps 5` for small effects.

## Hotspots in the 3.4 baseline

These are good optimization targets. Numbers come from `baseline-M` and the L run. The machine was a 16-vCPU Xeon Platinum 8468 VM with gcc 13 at -O3, using SSE.

- **hmmsearch thread scaling stops at about 4 threads.** `hmmsearch_p100_50k_cpu{1,2,4,8,N}` took 17.4, 10.6, 8.6, 8.4 and 9.6 s. CPU utilization peaks at 2.9×. The single reader/master thread is the bottleneck, and 16 threads are slower than 8.
- **hmmscan barely uses its threads.** `hmmscan_h1000_pfam` ran for 267 s at 1.9× CPU, and `hmmscan_h100_pfam` at 1.8×.
- **The nhmmer FM-index path is single-threaded and about 8× slower than plain FASTA.**
  - `nhmmer_alu_5M_fm` took 8.3 s against 1.0 s for `nhmmer_alu_5M`.
  - `nhmmer_dfam100_5M_fm` took 877 s against 108 s, at 1.0× CPU.
  - Its hit sets also differ slightly, because the FM seed heuristics aren't the same as the plain path's.
- **Memory spikes on long or repeat-rich queries.**
  - `phmmer_titin_50k`: 37 GB RSS for a 34k-aa query.
  - `nhmmer_seqquery6k_5M`: 13 GB for a 6 kb repeat-rich DNA query.
  - `nhmmer_l1_5M`: about 1 GB for a single L1 model.
- **Input-parsing overhead.**
  - gzip input is about 1.9× slower than plain FASTA on full Swiss-Prot, because of single-threaded `gzip -dc`.
  - EMBL/UniProt text is about 2.6× slower than FASTA (50k sequences).
- **Kernel throughput, M=400 and L=400.**
  - SIMD: MSV 46.7 Gc/s, Viterbi filter 2.8 Gc/s, Forward/Backward 121 Mc/s.
  - Generic: MSV 345 Mc/s, Viterbi 164 Mc/s, Forward/Backward 16 Mc/s.

