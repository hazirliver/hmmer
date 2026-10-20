#!/usr/bin/env python3
"""Run several hmmsearch snapshots (e.g. one per phase) on the same workloads, in one session.

  dev/perf/series.py SNAP1 SNAP2 ... [--labels L1,L2,...] [--workloads REGEX] [--reps 3] [--out FILE]

SNAPs are anything ab.sh accepts (binary, dir, or snapshot name under bench-results/bins/).
Workloads come from workloads.tsv, or --file (all of them by default; --workloads filters by name).

Each round runs every binary once per workload, in an order rotated by one each round, so
machine drift spreads over all binaries; the first round is an untimed warm-up whose outputs
are the reference: every binary's files must equal the first binary's after the dev/bench
normalization (whitespace-only differences accepted, as in ab.sh). A workload that uses an
option a binary doesn't have (--dbcache before Phase 2) is skipped for that binary.

For pyhmmer (dev/perf/pyhmmsearch.py, dev/profiling/PYHMMER.md): --mask-tail drops the block that
hmmsearch appends to its tables, and --flush-subnormals compares numbers below 2.2e-308 as 0 (hmmsearch
flushes them to zero, pyhmmer doesn't), noting how many each binary had.

Writes JSON (default bench-results/series/<time>.json): per workload, per binary label, the
samples of wall, user, sys and peak RSS, plus cycles:u, cycles (user+kernel) and instructions:u
for workloads judged by counters (the cpu and reader sets; --perf-all for all of them: perf stat
inflates threaded wall times, see ab.COUNTER_SETS).
dev/perf/plot_series.py draws it.
"""
import argparse, datetime, json, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import ab      # noqa: E402  (runs, hashing, workloads)
bench = ab.bench   # dev/bench/bench.py: data dir, results dir, normalization
# the lines of p7_tophits_TabularTail(), which ends every --tblout / --domtblout (a lone '#' occurs only there)
TAIL = re.compile(r"^#(?: (?:Program|Version|Pipeline mode|Query file|Target file|Option settings|Current dir|Date):.*| \[ok\])?\n?$")


def supports(binary, option):
    r = subprocess.run([binary, "-h"], capture_output=True, text=True, timeout=60)
    return option in r.stdout


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("snaps", nargs="+")
    ap.add_argument("--labels", help="comma-separated labels, one per snapshot (default: snapshot names)")
    ap.add_argument("--workloads", help="regex on workload names")
    ap.add_argument("--file", default="workloads.tsv", help="workloads file, in dev/perf or a path [workloads.tsv]")
    ap.add_argument("--reps", type=int, default=3, help="timed rounds [3]")
    ap.add_argument("--warmup", type=int, default=1, help="untimed rounds first [1]")
    ap.add_argument("--timeout", type=float, default=3600)
    ap.add_argument("--out")
    ap.add_argument("--perf-all", action="store_true",
                    help="perf stat on every workload, not only those judged by counters (inflates threaded wall times)")
    ap.add_argument("--mask-tail", action="store_true",
                    help="drop the '# Program: ... # [ok]' block hmmsearch appends to its tables before comparing "
                         "(pyhmmer's tables don't have it; dev/perf/vs_pyhmmer.tsv)")
    ap.add_argument("--flush-subnormals", action="store_true",
                    help="compare numbers below 2.2e-308 as 0, as hmmsearch prints them (pyhmmer doesn't flush "
                         "subnormals to zero); each binary's count is noted")
    a = ap.parse_args()
    if a.mask_tail:
        ab.EXTRA_RULES.append((TAIL, ""))
    ab.FLUSH_SUBNORMALS = a.flush_subnormals

    bins   = [ab.resolve_bin(s) for s in a.snaps]
    labels = a.labels.split(",") if a.labels else [ab.describe_bin(b).get("name", s) for b, s in zip(bins, a.snaps)]
    if len(labels) != len(bins): ab.die("--labels: one label per snapshot")
    ncpu = os.cpu_count()
    fixtures, wls = ab.load_workloads(ab.workloads_file(a.file), ncpu)
    wls = [w for w in wls if not a.workloads or re.search(a.workloads, w["name"])]
    if not wls: ab.die("no workloads selected")
    data = bench.derived_dir()
    ab.check_data(wls)
    ab.make_fixtures(fixtures, [w["args"] for w in wls], data, os.path.dirname(bins[0]))
    has_dbcache = [supports(b, "--dbcache") for b in bins]
    use_perf = ab.perf_ok()

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    outdir = os.path.join(bench.RESULTS, "series")
    work = os.path.join(outdir, "work", stamp)
    os.makedirs(work, exist_ok=True)
    out = a.out or os.path.join(outdir, f"{stamp}.json")
    res = dict(date=stamp, machine=bench.sh("lscpu | sed -n 's/^Model name: *//p'"), ncpu=ncpu,
               reps=a.reps, warmup=a.warmup, perf=use_perf,
               binaries=[dict(label=l, **ab.describe_bin(b)) for l, b in zip(labels, bins)], workloads=[])
    failed = False

    for w in wls:
        args = ab.expand(w["args"], data)
        needs_cache = "--dbcache" in w["args"].split()
        runs = {l: [] for l in labels}
        ref, ref_label, notes, problems = None, None, [], []
        print(f"\n{ab.BLD(w['name'])}  {ab.DIM(w['args'])}", flush=True)
        for r in range(a.warmup + a.reps):
            order = list(range(len(bins)))
            order = order[r % len(order):] + order[:r % len(order)]
            line = []
            for i in order:
                if needs_cache and not has_dbcache[i]: continue
                cwd = os.path.join(work, w["name"], f"{labels[i]}_{r}")
                os.makedirs(cwd)
                s = ab.run_once(bins[i], args, cwd, use_perf and ab.wants_perf(w, a.perf_all) and r >= a.warmup, a.timeout)
                if s.get("rc") != 0:
                    problems.append(f"{labels[i]} round {r}: exit {s.get('rc')}"); failed = True
                if r < a.warmup:
                    h = ab.hash_run(cwd)
                    if a.flush_subnormals:
                        n = {f: ab.count_subnormals(os.path.join(cwd, f)) for f in h}
                        if sum(n.values()):
                            notes.append(f"{labels[i]}: {sum(n.values())} numbers below 2.2e-308 compared as 0 ("
                                         + ", ".join(f"{f} {k}" for f, k in sorted(n.items()) if k) + ")")
                    if ref is None: ref, ref_label = h, labels[i]
                    else:
                        kind, bad = ab.compare(ref, h)
                        if kind == "ws":   notes.append(f"{labels[i]} vs {ref_label}: whitespace only ({', '.join(bad)})")
                        if kind == "diff":
                            problems.append(f"{labels[i]} vs {ref_label}: {', '.join(bad)} DIFFER"); failed = True
                            continue                                  # keep the run dir
                else:
                    runs[labels[i]].append(s)
                shutil.rmtree(cwd)
                line.append(f"{labels[i]} {s['wall']:.2f}s")
            print("  " + "  ".join(line), flush=True)
        for n in notes:    print(ab.YEL("  note: " + n))
        for p in problems: print(ab.RED("  " + p))
        res["workloads"].append(dict(name=w["name"], args=w["args"], sets=w["sets"], runs=runs,
                                     identical=not any("DIFFER" in p for p in problems), notes=notes, problems=problems))
        with open(out, "w") as f: json.dump(res, f, indent=1)   # after every workload: survives an interruption

    print(f"\nresults: {out}")
    print(ab.RED("FAILED: a run failed or outputs differ") if failed else "outputs identical across all binaries")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
