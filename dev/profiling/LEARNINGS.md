# Lessons from the hmmsearch speedup work (Phases 0-3)

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

## Byte-identical rewrites of floating-point code

1. **Same operations, same order, same types.** IEEE `+` and `*` are commutative but not associative, so keep the evaluation tree. Keep `(float) log((double) x)` and float accumulators as they were.
   - `0.0f + a` equals `a` for `a ≥ +0`. A term that is always `+0.0` can be dropped exactly only when the accumulator can't be `-0`.
   - Kahan summation (`esl_vec_FSum`) over two elements is a plain sum. Over more elements it must be replicated exactly.
2. **Check what the compiler emitted.** It could narrow `(float)log((double)x)` to `logf`, which isn't bit-identical. Verify with `objdump -d | grep call`. FMA contraction only appears with `-mfma`/`-march=native`, which the SSE goldens don't use. Mirror the original statement structure anyway (`n = 0; n += a*b; n += c*d;`).
3. **Three layers of evidence for a kernel rewrite:**
   1. a permanent bitwise unit test against the old routine (random inputs and edge lengths 0, 1, 2; all residue codes);
   2. a mutation test showing the unit test fails on a slightly non-identical rewrite (`x/max` → `x*(1/max)` was caught at a 6-ULP difference);
   3. a `-Dp7_CROSSCHECK` build that runs old and new code side by side and aborts on any bit difference, run over the goldens and all real-data workloads.

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
- **Generate goldens for new coverage with the unmodified binary**, and commit them separately ("coverage only").
- **New feature goldens get cross-checked** against the equivalent existing goldens, minus the option's echo lines, before they're committed.
- **New synchronization gets a ThreadSanitizer run and a stress loop.** TSan on this kernel needs ASLR off: configure, build and run under `setarch -R` (otherwise "unexpected memory mapping"). Check that it flags a deliberately racy toy program first. Then loop the threaded binary under `timeout` at several thread counts to catch hangs.
- **Manual memory management gets a sanitizer run.** Valgrind isn't installed here, so use `-fsanitize=address` in a separate build dir and confirm with `ldd` that libasan is linked.

## Tooling gotchas (this machine)

- `pkill -f <pattern>` matches the shell running it when the pattern appears in its own command line, which kills the command. Use `pgrep -x <name>`.
- `argparse` options taking a value that starts with `--` need the `=` form: `--b-args=--dbcache`.
- Man pages are generated by `configure`. After editing `*.man.in`, run `./config.status --file=<out>:<in>` in the build dir.
- `git worktree` plus a symlink to the pinned `easel/` builds any commit without touching the working tree (`snap.sh --rev`).
- A background shell task is capped at 2 h. Long benchmark tiers have to be split.
- `gdb` isn't installed. To prove a code path runs, look for one of its non-inlined functions in a `perf record` profile.
- `perf record -g` makes thread creation look far more expensive than it is (counter setup per new thread). Profile thread-creation-heavy runs without `-g`, or after the threads are made persistent.
