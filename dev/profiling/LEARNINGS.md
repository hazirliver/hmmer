# Lessons from the hmmsearch speedup work (Phases 0-7)

These are practical lessons from doing [`PLAN.md`](PLAN.md) step by step. The measured numbers are in [`PROGRESS.md`](PROGRESS.md). Read this before measuring or changing anything.

## Measuring

1. **Pick the metric that can actually resolve the effect.** Measured on this 16-vCPU VM:

   | Metric | Typical spread | Use it for |
   |---|---|---|
   | `instructions:u` at `--cpu 0` | 0.00% | Proving work was added or removed. It's deterministic. |
   | `cycles:u` at `--cpu 0` | 0.2-0.7% | Single-thread CPU steps |
   | Wall time, threaded, interleaved A/B | 1.4-12% | Scaling steps; it needs big effects or many reps |
   | Wall time vs. an old baseline file | Drifts up to ~10% day to day | Nothing. The Phase 1 bench M run showed untouched kernels "12% faster" than the stored baseline |

2. **Always interleave A and B in the same session** (`dev/perf/ab.sh`). Never judge a step against a results file from another day.

3. **Code layout moves cycles by 1-3% without any real change.** Editing `p7_bg.c` shifted every later function by 96 bytes. That changed the 64-byte alignment of all SSE kernels. A workload that never calls the edited code (W5, `--max`) got 1.1-1.6% "slower" with an identical instruction count.
   - Step A/Bs build both sides with `-O3 -falign-functions=64` (`snap.sh --perf`), which removed the artifact.
   - Default `-O3` builds still show ±1.4% in both directions at phase milestones. Read those against the instruction counts.

4. **Identical instructions plus different cycles means re-run before concluding.** The aligned W3 was once +1.3% cycles with an identical instruction count. On the next step it read 1.0000.

5. **Profile shares can mislead.** After the bias-filter rewrite, libm's share of the profile rose from 3.9% to 7.1% while total cycles fell. Stalls just land on different instructions. Judge by totals.

6. **Know which resource bounds the workload before predicting a wall-time gain.** While the serial reader set the pace at 16 threads, a 2.7% instruction cut showed up as about 1.6% wall time on the realistic run (Phase 1, within noise). The same run got 2.8× faster once the reader was fixed (Phase 2).

7. **`-falign-functions=64` doesn't remove every layout effect.** Step 3.1 edited `p7_pipeline.c`, and the serial X3 read +1.2% cycles with *fewer* instructions, cache misses, TLB misses and page faults. To separate layout from cost, build a **control**: B's exact library code with the new behaviour switched off at the call site (here: destroy and recreate the pipeline anyway). The control was +1.5% vs A, so the change itself was −0.5%. Compare B against the control, not only against A.

8. **Short runs (< 0.5 s) have noisy wall times.** W4 showed +13.7% wall time in `ab.sh` where a direct `perf stat -r 15` measured +4.6%. For short workloads use counters, or `perf stat -r N`.

9. **Find the real limit before optimizing the planned one.** Per-thread `perf script -F tid,period` sums showed, in turn: thread creation (sys), the wake-up storm (kernel spinlocks), the single helper (64% of samples), then latency with no thread saturated (4.7 of 16 cores busy). Each pointed to a different fix than the plan had.

10. **Count how often a shortcut fires before adding it.** An exact early-out that skips a divide (4.8, `nov <= 0` in the clustering link test) almost never fired on real data. It added +0.15% instructions and cost 2% cycles on W5. Likewise, removing an integer divide that was off the critical path (4.3) gained nothing measurable. A quick counter in a throwaway build answers this before any A/B.

11. **Check the measuring tool's clock.** `ab.py` timed runs around `Popen.wait(timeout=...)`, which polls with sleeps growing to 50 ms: every wall time was rounded up to the next poll (X2 samples sat on a 0.865/0.915/0.965/1.016 s ladder, ~5% of a 0.9 s run). A blocking `wait()` plus a timer that kills on timeout gives the exact exit time. `dev/bench` had the same flaw.
12. **On this VM the clock rate varies between runs.** On 2026-10-05 R1's task-clock spread was 8-14% and X2's wall 16-25% while cycles stayed within 1-2%. For single-threaded steps that move work between user and kernel (system calls), judge by user + kernel `cycles` (ab.sh `reader` set, `cycles_uk`), not task-clock.
13. **Loop alignment can hide a small gain completely.** 5.1 removed 80k `lseek`s (kernel cycles −9%), but `seebuf()`/`header_fasta()` lost 2% user cycles because an inlined callee changed and moved their loops; a "control" build landed 2% *below* A, so ±3% was layout. `-falign-functions=64` doesn't align loops. A real −0.5% can't be resolved that way; it was dropped (5.2 cut the calls 4× anyway).
14. **Bigger buffers aren't automatically faster.** 16 KB reads beat 64 KB-1 MB by 1-7% cycles because the old parser made two passes over each buffer and needs it in L1; the plan had guessed 64-256 KB.
15. **A large share in a library call can be stall time, not the call's cost.** After Phase 6, libm `log()` held 6.5% of the realistic run's samples, called from the bias filter. An exact inline log (verified on all 2^32 inputs) saved only ~1%: the filter's critical path is its `mul → add → max → div` chain, and `log` runs in its shadow while the samples land on it. Before rewriting a hot call, test the hypothesis cheaply: put a stand-in into a throwaway build and judge it by cycles.
16. **Measure the barrier, not only per-thread CPU.** Per-thread cycle sums made uncached X5 look reader-bound. Timestamps per pass (start, each worker's end, time inside `esl_sqio_ReadBlock`) showed workers idle 39% of each query in its tail and the reader idle 40%, waiting for those tails. A ~10-line throwaway patch answers this.
17. **Hyperthread siblings make a quick latency-vs-throughput probe.** Pin two copies to sibling CPUs. Each running ~1.0× slower means latency-bound code, where interleaving independent work in one thread pays. ~2× slower means port-bound. Measured: Forward/Backward and traceback 1.38×, SSV 1.58×.
18. **Check float functions exhaustively when the input is a float.** All 2^32 bit patterns take about a minute. That turns "within 1 ulp, so the rounding should agree" into a proof for the libm in use. Keep a guard band with a fallback for other libms.
19. **A mode that wins on small inputs can lose badly on big ones.** The query-parallel mode (3.5b), forced on X5, was 1.53× slower: output is in query order, so one long model blocks the job window. Test a mode on the input sizes it would newly cover before widening its threshold.

## Byte-identical rewrites of floating-point code

1. **Same operations, same order, same types.** IEEE `+` and `*` are commutative but not associative, so keep the evaluation tree. Keep `(float) log((double) x)` and float accumulators as they were.
   - `0.0f + a` equals `a` for `a ≥ +0`. A term that is always `+0.0` can be dropped exactly only when the accumulator can't be `-0`.
   - Kahan summation (`esl_vec_FSum`) over two elements is a plain sum. Over more elements it must be replicated exactly.
   - Inlining a library routine is exact if its statements are copied verbatim (4.1-4.2: `esl_vec_FNorm` + `esl_rnd_FChoose` as one `static inline` with a constant `n`). The compiler may only apply value-preserving rewrites, and without `-mfma` it can't contract.
   - **A dependent add chain can be broken exactly by reordering loops, not sums** (4.6): run q outer and residues inner, so the per-residue sums become independent chains while each still adds in ascending q. Keep the block width a compile-time constant so the accumulators stay in registers, and check the disassembly for spills.
2. **Check what the compiler emitted.** It could narrow `(float)log((double)x)` to `logf`, which isn't bit-identical. Verify with `objdump -d | grep call`. FMA contraction only appears with `-mfma`/`-march=native`, which the SSE goldens don't use. Mirror the original statement structure anyway (`n = 0; n += a*b; n += c*d;`).
3. **Three layers of evidence for a kernel rewrite:**
   1. a permanent bitwise unit test against the old routine (random inputs and edge lengths 0, 1, 2; all residue codes);
   2. a mutation test showing the unit test fails on a slightly non-identical rewrite (`x/max` → `x*(1/max)` was caught at a 6-ULP difference);
   3. a `-Dp7_CROSSCHECK` build that runs old and new code side by side and aborts on any bit difference, run over the goldens and all real-data workloads.
4. **Random inputs rarely hit the boundary of an RNG-driven choice.** To test a sampler bitwise, peek the next roll with a third RNG of the same seed and build the probabilities so a cumulative sum lands on it (4.1). Then rounding decides the outcome, and a `<` vs `<=` or float-vs-double mutant is caught within a few trials.

## Rewriting a parser byte-identically (Phase 5)

1. **Keep the old implementation as a verbatim, renamed copy and drive both through the public API.** `testsuite/sqio_ascii_ref.c` is upstream `esl_sqio_ascii.c` byte for byte; the test driver `#include`s it with every function renamed `ref_*` by macros. `ref_esl_sqascii_Open()` then installs the old routines in an `ESL_SQFILE`'s dispatch table, so every `esl_sqio_*` call (Read, ReadBlock, ReadWindow, Position, Echo, SSI fetches, guessing, Parse) runs on either parser with identical arguments. Compare everything that leaks out: return codes, error messages, line numbers, bookkeeping fields, every `ESL_SQ` field.
2. **Make the buffer size a variable in the test build and run down to 1 byte.** A boundary-handling mistake then shows on almost every file instead of once per 16 KB. The unchanged parser at sizes 1-4095 against the 4 KB reference proved the old code boundary-independent where it matters, and found where it isn't (next item).
3. **Old code can depend on where reads end; find those cases before changing the read size.** Found here: on pipes `ftello()` fails, so record offsets are -1 plus the position in the current read; daemon-format files read a byte past the buffer at record ends; with an alphabet that disagrees with the format's input map, `seebuf()` and `addbuf()` fall out of step and which residues survive depends on the chunking; `skip_whitespace()` indexes the input map with a negative `char`. Keep such inputs at the old read size where that's cheap (pipes, daemon format), and leave the undefined ones out of the comparison with the reason written next to the exclusion.
4. **ASan the differential test, not only the product.** It found an upstream line-buffer overflow that only bigger reads trigger (`if` instead of `while` when growing; fixed in 5.2) and the negative-index read. A test that compares garbage-dependent results flakes; the sweep over seeds shows it as rare, irreproducible failures.
5. **Make the generator reach each path's preconditions.** Two mutants survived the random files: advancing an offset by the read size instead of the byte count (needs a header at EOF after a partial read) and not restoring the line number on hand-over (needs the disagreeing byte after the first lines). Hand-written edge-case files run at every read size, and generators that place features deliberately, caught them.
6. **A fast path needs an exact hand-over.** `seeaddbuf()` returns `eslENORESULT` with all state restored when it can't stand in for the old pair (line-based formats, another record terminator, a disagreeing alphabet), and the caller runs the old code. The hand-over itself needs a test that reaches it, or it is dead weight that can be wrong.
7. **Re-profile after each step, down to the instruction.** After 5.4 the header parser still took 31% of the reader; `perf annotate` put 85% of it on one byte loop (the description scan, two taken branches per byte, ~3 cycles per character on ~210-byte Swiss-Prot headers). Replacing it with `memchr()` (5.4b) gave another −11% on R1 and X2, more than 5.4 itself.
8. **A per-byte lookup can be cached if the cache checks what it was made from.** 5.6 replaced two map lookups and three range checks per residue with one table lookup, 8 residues per test. Maps can be edited in place by callers, so the cache keeps copies of both 128-byte maps and compares them on every call (one `memcmp()` per record: noise); a test edits a map between reads, and a mutant without the comparison fails it. Write 8-byte stores with the compiler's byte order (`__BYTE_ORDER__`), not a byte array: GCC didn't merge the array's stores, which gave back a third of the gain.

## Wider SIMD for exact integer DP (Phase 6)

1. **Prove first that the lane count can't reach the result.** Striped integer DP is layout-independent only where every cell is the same function of the same predecessors. List what *does* depend on the layout and show it can't leak:
   - padding nodes (their scores must keep them at the floor: SSV 127 saturates to -128, MSV 255 gives 0, VF -32768 stays -32768);
   - work done after an overflow (SSV ORs every shifted lane with 0x80; harmless only because the overflow is already decided by then);
   - shortcuts whose trigger is computed per layout (the Viterbi filter's DD bound uses the row's max M->D value, which is the same in every layout; the D values it leaves are either exact or exactly the M->D values).
   Write the argument in the kernel's header; it tells you which mutants and test targets matter.
2. **Fewer instructions can be more cycles in a row recurrence.** The first AVX2 full MSV had 3% fewer instructions on W3 and was 25% slower than SSE in that kernel: each row's max went to a scalar and xB came back by broadcast, and that round trip is on the critical path from row to row (row i needs xB(i-1)). Keeping the reduction in vector registers (halves swap, then in-lane shuffles, as the SSE code does) made it 22% faster than SSE. For per-row reductions, count latency, not instructions.
3. **Per-query preparation shows up on many-models × tiny-database runs.** Restriping the profile for the wide kernels is O(Kp × M) per query: invisible on X5, but 3.6% of X4 serial in its first version (wrapping counters: a branch miss every few nodes), 1.9% after copying the word scores too. Fixes, in order: nested loops / precomputed offset tables (no data-dependent branches), then making the rarely used part lazy (the Viterbi copies, needed by ~2% of targets). Always run the query set (X4) for a step that adds per-model work, and measure it in cycles: by wall the 1.9% was inside X4's noise.
4. **Lazy shared state without blocking.** Clones share the copies, so the "built" flag lives in the shared block, not in the (per-clone) struct. The first caller wins a compare-and-swap, builds, and publishes with a release store; readers use an acquire load. Callers that find it being built use the old (SSE) path, which gives identical results, so nobody waits. Tests: a re-conversion test (a profile reused for another model must rebuild; a mutant without the reset failed), a threaded test with clones racing to build (a mutant that used half-built copies failed it), TSan on both (clean; it flagged the racy mutant), and threaded hmmsearch under TSan.
5. **Random tests miss exact boundaries; scans hit them.** Mutants found two holes the random targets (and the goldens) didn't cover: an off-by-one in the MSV overflow test (`>` for `>=`: needs a row max landing exactly on 255 - bias) and dropped M->D values at segment ends (needs a deletion at exactly that node). Both were caught by the crosscheck build on real data, and then by scans built for them: a long random flank plus every prefix of the consensus (the score climbs through the boundary one residue at a time) and the consensus with 1-8 nodes deleted at every position.
6. **Random test models can be unrealistic in ways that hide paths.** `p7_hmm_Sample()` gives weak M->M transitions, so the Viterbi filter never overflowed on DNA; `p7_hmm_SampleUniform()` (Pfam-like transitions, random emissions) fixed it. Make the test assert that every outcome (ok, overflow, J state, possible overflow) was reached, so a generator change can't silently remove coverage.
7. **One binary, two paths.** With a run-time switch (`HMMER_SIMD`) and `ab.sh --a-env/--b-env`, the old and new kernels run from the same binary: no code-layout difference at all between A and B. The crosscheck can use the same switch: rerun the old path on a struct copy of the profile with the new feature off, compare status and score bits, abort on a difference. No refactoring of the old code is needed.
8. **The build must notice what changed.** `src/*.o` didn't depend on `impl_sse.h` (a struct change left stale objects in incremental builds), neither `dev/check.sh` nor `snap.sh` reran configure when a `Makefile.in` changed, and `snap.sh`'s dirty hash ignored new untracked files. A snapshot built that way silently ran the old code (it passed a negative control it should have failed). After any build-system or header change, check that the binary has the new code: a new symbol in `nm`, and a negative control (`HMMER_SIMD=bogus` must fail).
9. **Wider isn't faster by default: judge a new vector width by wall time at full thread count.** The AVX-512 kernels ran 9-25% fewer instructions than AVX2 and were still slower on most workloads: on Sapphire Rapids 512-bit integer ops issue on 2 ports (AVX2 on 3), fixed per-residue work (the SSV's cross-lane shift, the MSV's row reduction) is spread over half as many vectors, and at 16 threads the clock ran 1.0-1.6% lower (cycles per task-clock in `ab.sh`'s JSON). Cycles alone hide that last part: a per-model choice (AVX-512 only where AVX2 needs several SSV passes, M > 448) was 1% fewer cycles on X5 and exactly as fast in wall time. Break kernels down by model length before choosing (`perf report` share × total cycles, per level, in one binary).
10. **Check a bottleneck hypothesis with a cheap experiment before building on it.** The AVX2 SSV at M=200 ran ~7 cycles per residue for a band of 7 vectors, which looked like the single max accumulator's dependency chain; splitting it in two changed nothing (+0.6% instructions). The vector ports were the limit.
11. **Run every set at the milestone, even when each step's A/B looked clean.** The step A/Bs ran each step's target sets; the first Phase 6 milestone (all 21 workloads) found two regressions none of them showed: R1 (one short model, M = 50: the AVX2 SSV with only two vectors is latency-bound on its cross-lane shift; 6.4c) and per-query setup on X4 (6.4b, below). Model length changes which kernel wins; keep a short-model workload in the sets.
12. **A constant that becomes a run-time value can change the generated code completely.** When 6.4 made the vector width a variable, 16 small `memcpy()`s per profile row compiled to `rep movsq` (large start-up cost) instead of two vector moves: +1.7% on X4. `perf annotate` pointed at one instruction. Copy only what is needed, in one piece, or keep sizes compile-time constants per kernel.


## Scheduling queries and pairing targets (Phase 7)

1. **The barrier was the bound, not the kernels.** Phase 6 left X5 at 16 threads with workers busy 61% of each query; the per-query tail (one slow hit-bearing target) and, uncached, the reader waiting out each tail cost more than any kernel. Searching a batch of queries per pass over the targets (7A.2, 2.3× on X5 uncached) and then overlapping consecutive batches (7A.4, another 11%) fixed it. Measure idle time per pass before tuning kernels for a threaded run.
2. **Batch size is a memory decision too.** A batch's hit lists are held until it's printed: 10 hit-rich queries × Swiss-Prot with `-A` went 285 → 366 MB at 16 queries per pass. The sweep (4-64) flattened after 16-32; 16 keeps the worst case (huge families × big databases) at half of 32's.
3. **Moving work changes what's printed before a fatal error.** Printing a query's header after its pass instead of before it (7A.2) changed two error paths no golden covered: `--cut_ga` failing on the very first model ("Unexpected error 0 reading sequence file") and a parse error in the *target* file (the failing query's header was missing). `dev/perf/errpaths.sh A B` now compares two binaries on HMM-file errors, target-file errors and batching cases at several thread counts; run it for every scheduling change, and add a golden for each new case, checked against upstream.
4. **Shared scratch must not be reset by setup.** Setting up the next batch's slots with `p7_pipeline_Reinit()` reset the DP matrices the worker was using for the current batch. Configure slots from a template copy instead, and let the worker itself do per-batch resets (`worker_Fresh()`). Check first what results really depend on: the RNG is reseeded for every region, so where its state comes from doesn't matter.
5. **TSan also reports races that write the same value.** The master rewrote each worker's `reuse_seqs` for every batch while streaming workers read it. The value never changed; still a race. Set such fields once, while the workers are parked.
6. **Deferring a target's later stages means restoring everything they read.** To pair a waiting target with the next one, its own length configuration (`om->xf`, `xw`, `tjb_b`, `L`; `bg->p1`) and the pipeline's running Z have to be put back before its domain stage. List what each later stage reads; for the parsers it's only `xmx`, `L`, `has_own_scales` and `totscale` of their matrices.
7. **Two targets per 256-bit register is exact where a wider stripe isn't.** A wider stripe would change the float summation order; one target per 128-bit half keeps each target's SSE operation sequence: in-lane shifts and shuffles, the same operand nesting, per-target specials, a blend instead of `x * 1.0` for the half that isn't rescaled, and per-target early exits (zero the stopped half's carry). Gains are bounded by how often two targets need the parsers: 1.4-1.9× per pair, but only ~0.2% of targets reach them on the realistic run.
8. **Read granularity is a cache decision.** Reading 1000 targets before searching them (to pair within a block) made a reading-heavy serial run 5% slower in cycles with fewer instructions, and allocating that block per pass cost 4% on 22k-query runs. Reading one target at a time and keeping only the waiting ones alive (a small pool) gets the pairs without either cost.
9. **Some mutants survive because the case is rare, not because the code is wrong.** Dropping the per-target stop of the paired Forward's D->D passes survived thousands of random and deletion-heavy pairs: an extra pass changes a stopped target only for a D value near the top of its binade receiving a carry just below the next value's half-ulp. Keep such guards when they are the original code's behaviour, write the argument down, and don't count the test as covering them.

## Determinism

1. **Upstream threaded hmmsearch isn't byte-deterministic without `-Z`.** Each worker uses its own running sequence count as Z (`p7_pipeline.c:580`), so the set of *unreported* kept hits depends on scheduling. That set feeds column widths (`p7_tophits.c:1134`). Reported results never change.
   - Decision: whitespace-only differences are acceptable (`ab.sh` reports them as `ws`).
   - Strict threaded workloads pin `-Z`.
   - The small golden cases were stable 30/30.
2. **Domain definition reseeds its RNG per region**, and hits are sorted by score, then name, then position. Reordering or redistributing targets across threads can't change results. This made the cache and threading work safe.

3. **Scheduling-only changes still change *how often* nondeterminism shows.** Splitting blocks smaller (3.3c, first version) made the upstream width variation appear in 9 of 10 threaded runs instead of 0 of 10: each worker's running Z was smaller, so more borderline hits were kept. Check a threaded no-`-Z` case 10 times before and after.

## Design patterns that worked

1. **Opt-in features keep the default path's code literally unchanged.** A per-worker flag (`reuse_seqs`) was enough. Then prove it twice: instruction counts identical on the default path, and `dev/check.sh --force-opt <opt>` reruns every golden case with the option forced on.
2. **Change the payload, not the protocol.** The `--dbcache` reader keeps the work queue's fixed pool of block tokens. A pool block becomes a *carrier* whose `list` points at cached sequences while in flight, and its own list is restored after the pass. There's no queue surgery and no ownership confusion at cleanup.
3. **Make the first pass identical to the uncached path, then copy.** Reading into the usual reused block and appending a compact copy had several payoffs:
   - outputs are trivially the same;
   - workers can reset blocks as before;
   - exceeding the memory limit mid-run just means dropping the copies.
4. **The allocator was the overhead, not the copying.** Fresh 1000-sequence blocks with Easel's default growable buffers (32/32/128/256 B, plus doubling) cost +39% wall time and 653 MB.
   - An arena with exact sizes brought it to +26% and 547 MB.
   - 2 MB-aligned chunks with `madvise(MADV_HUGEPAGE)` cut first-touch page faults on the serial reader and brought it to +15%.
   - What's left is the data itself: residues, 120 MB of Swiss-Prot headers and about 230 B of `ESL_SQ` per sequence.
5. **Reused objects must be reset to *exactly* the created state.** Share one configure routine between `Create` and the reset (`p7_pipeline_Reinit()`), so they can't drift apart. Watch for state that `Create` gets for free:
   - Grow-only DP matrices keep the row stride of the largest model seen. Smaller later queries then walk a sparser matrix (+0.7% cycles on X3) until the row pointers are invalidated (`validR = 0`).
   - `p7_ProfileConfig` copies RF/MM/CS only when the HMM has them. A reused `P7_PROFILE` must go through `p7_profile_Reuse` first, or the previous query's RF line leaks into the next query's alignments.
6. **Gotcha:** on the protein path, `esl_sqio_ReadBlock` does **not** reset block slots. `sqascii_Read` expects them already `esl_sq_Reuse`'d (only the DNA branch resets them itself). Whoever owns a block between reads must reset it.
7. **Don't switch execution modes per item when switching has a barrier cost.** Per-query mode choice put a drain before every switch: 1.52× slower on X5. Decide once per run from something fixed (db size).


## Verification habits

- **Run a negative control for every new checker.**
  - `ab.sh` with a B that adds `--seed 7` must report DIFF.
  - `--force-opt --no-such-opt` must fail every hmmsearch case.
  - The bitwise unit test must fail on a mutation.
  - A crosscheck must fire on a mutant. In 4.3 it didn't: `build-xchk` had `-Dp7_CROSSCHECK` hand-added to `src/Makefile` only, so `src/impl_sse/` compiled without it. Configure crosscheck builds with `CFLAGS="-O3 -Dp7_CROSSCHECK"`, never by editing generated Makefiles (`config.status` also silently reverts such edits).
- **The reference in a crosscheck or unit test must be the original code**, not a wrapper around the new code. In 4.7 the Easel-callback reference first wrapped the new predicate, so it couldn't have caught a mistake in it (noticed while trying 4.8). It now holds the upstream body verbatim.
- **Generate goldens for new coverage with the unmodified binary**, and commit them separately ("coverage only").
- **New feature goldens get cross-checked** against the equivalent existing goldens, minus the option's echo lines, before they're committed.
- **New synchronization gets a ThreadSanitizer run and a stress loop.** TSan on this kernel needs ASLR off: configure, build and run under `setarch -R` (otherwise "unexpected memory mapping"). Check that it flags a deliberately racy toy program first. Then loop the threaded binary under `timeout` at several thread counts to catch hangs.
- **Manual memory management gets a sanitizer run.** Valgrind isn't installed here, so use `-fsanitize=address` in a separate build dir and confirm with `ldd` that libasan is linked.

## Tooling gotchas (this machine)

- `pkill -f <pattern>` matches the shell running it when the pattern appears in its own command line, which kills the command. Use `pgrep -x <name>`.
- `argparse` options taking a value that starts with `--` need the `=` form: `--b-args=--dbcache`.
- New unit-test drivers go into `UTESTS` in `src/Makefile.in` and `testsuite/testsuite.sqc`. Existing build dirs need `./config.status src/Makefile` to pick them up.
- Man pages are generated by `configure`. After editing `*.man.in`, run `./config.status --file=<out>:<in>` in the build dir.
- `git worktree` plus a symlink to the pinned `easel/` builds any commit without touching the working tree (`snap.sh --rev`).
- A background shell task is capped at 2 h. Long benchmark tiers have to be split.
- `grep` in a UTF-8 locale treats a file with invalid UTF-8 (some Swiss-Prot descriptions, so most hmmsearch outputs) as binary and drops its lines, silently with `-h` across several files. That made identical outputs look different. Compare outputs with `LC_ALL=C grep -a`, or with `cmp`/`md5sum` on whole files.
- libm has no frame pointers: its samples are charged to the caller's *caller* (bias filter `log()` showed up as "pipeline glue" in `p7_Pipeline`). Find the real caller with `perf record --call-graph dwarf` on a short run.
- `gdb` isn't installed. To prove a code path runs, look for one of its non-inlined functions in a `perf record` profile.
- `perf record -g` makes thread creation look far more expensive than it is (counter setup per new thread). Profile thread-creation-heavy runs without `-g`, or after the threads are made persistent.
- HMMER's `hmmpgmd_ga`/`hmmpgmd_shard_ga` sqc tests bind fixed ports 51373/51374, inside Linux's ephemeral range (32768-60999). A long-lived outbound connection (here the Claude Code session's own HTTPS connection) can hold one, and the master dies with "bind error 98"; the tests then fail and leave a worker running. Rerun later, or run a copy on ports 61373/61374 (hmmpgmd requires > 49151).
- The Claude Code safety check refuses `rm -rf` of directories and `xargs sh -c '...'`; use fresh `mktemp -d` directories and a small script file for parallel runs.
- `snap.sh --rev` must build the Easel commit that revision pins (it used to symlink the current `easel/`, which silently gave an old A the new Easel).
- HMMER's Makefiles don't regenerate themselves, and `src/*.o` didn't depend on `impl_sse.h` (fixed in 6.0): after editing `configure.ac`, a `Makefile.in` or a header, make sure the build dir was reconfigured (`dev/check.sh` and `snap.sh` now do it) and prove the binary has the change (`nm` for a new symbol, or a negative control).
- The Bash tool's foreground commands stop after 10 minutes whatever `timeout` says; wait for long A/Bs through their background task's completion notice, not a polling loop.
