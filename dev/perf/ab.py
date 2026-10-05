#!/usr/bin/env python3
"""Interleaved A/B comparison of two hmmsearch binaries: output identity + speed.

  dev/perf/ab.sh A [B] [--set cpu|scale|query|milestone|all] [--workloads REGEX] [--reps N]

A and B are a binary, a directory holding one, or a snapshot name under
bench-results/bins/ (made by dev/perf/snap.sh; a unique prefix is enough).
B defaults to build-dev/src/hmmsearch.

For each workload (dev/perf/workloads.tsv) both binaries run alternately
(A B, B A, A B, ...) so machine drift hits both equally. Each run is wrapped in
`perf stat` (task-clock, instructions:u, cycles:u) and GNU time (wall, user,
sys, peak RSS). Every file a run writes, plus stdout/stderr, is hashed after the
dev/bench normalization; A and B must match exactly, and every rep of a binary
must match its first rep.

Verdict per workload: the metric is cycles:u for the `cpu` set, wall time
otherwise. A ratio B/A counts as faster/slower only when it differs from 1 by
more than both the metric's floor and the measured spread (max-min)/median of
either binary's reps.

Exit status: 0 ok; 1 outputs differ, nondeterministic or a run failed;
3 a slowdown with --fail-on-slowdown. See dev/profiling/PLAN.md §3.
"""
import argparse, datetime, hashlib, json, os, re, shutil, statistics, subprocess, sys, time

HERE = os.path.dirname(os.path.realpath(__file__))
ROOT = os.path.realpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "dev", "bench"))
import bench  # noqa: E402  (normalization, data dir, colors)

RESULTS = os.path.join(bench.RESULTS, "ab")
BINS = os.path.join(bench.RESULTS, "bins")
FIXTURES = os.path.join(bench.DATA, "ab-fixtures")
EVENTS = ["task-clock", "instructions:u", "cycles:u"]
# metric -> (label, smallest relative change that can count, lower is better)
FLOOR = {"cycles": 0.005, "instructions": 0.001, "wall": 0.02, "task_clock": 0.01}
OWN_FILES = {".perf", ".rusage"}
WS = re.compile(rb"[ \t]+")
RED, GRN, YEL, BLD, DIM = bench.RED, bench.GRN, bench.YEL, bench.BLD, bench.DIM


def die(msg):
    bench.die(msg)


# ------------------------------------------------------------------ inputs
def resolve_bin(spec):
    """Binary path, dir containing hmmsearch, or snapshot name/prefix in bench-results/bins."""
    cands = [spec, os.path.join(spec, "hmmsearch"), os.path.join(spec, "src", "hmmsearch")]
    if os.path.isdir(BINS):
        names = sorted(n for n in os.listdir(BINS) if n.startswith(spec))
        if spec in names: names = [spec]
        if len(names) == 1: cands.append(os.path.join(BINS, names[0], "hmmsearch"))
        elif len(names) > 1 and not os.path.exists(spec):
            die(f"snapshot prefix '{spec}' is ambiguous: {', '.join(names)}")
    for c in cands:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.realpath(c)
    die(f"can't find an hmmsearch binary for '{spec}' (path, dir, or name under {BINS})")


def describe_bin(path):
    info = os.path.join(os.path.dirname(path), "info.txt")
    if os.path.exists(info):
        kv = dict(l.split(None, 1) for l in open(info) if l.strip() and len(l.split(None, 1)) == 2)
        return {k: v.strip() for k, v in kv.items()} | {"path": path}
    src = os.path.realpath(os.path.join(os.path.dirname(path), ".."))
    try:
        sha = bench.sh(f"git -C {src} rev-parse --short=8 HEAD 2>/dev/null || true")
        dirty = bench.sh(f"git -C {src} status --porcelain --untracked-files=no 2>/dev/null || true")
    except Exception:
        sha, dirty = "", ""
    return {"path": path, "name": (sha + ("-dirty" if dirty else "")) or os.path.basename(src)}


def load_workloads(path, ncpu):
    fixtures, wls = [], []
    for ln in open(path):
        if not ln.strip() or ln.lstrip().startswith("#"): continue
        if ln.startswith("@fixture"):
            _, made, cmd = ln.rstrip("\n").split(None, 2)
            fixtures.append(dict(made=made, cmd=cmd)); continue
        name, sets, args = ln.rstrip("\n").split(None, 2)
        wls.append(dict(name=name, sets=sets.split(","), args=args.replace("$NCPU", str(ncpu))))
    return fixtures, wls


def expand(s, data):
    s = re.sub(r"(?<![\w/.])D/", data + "/", s)
    return re.sub(r"(?<![\w/.])F/", FIXTURES + "/", s)


def make_fixtures(fixtures, needed_args, data, bindir):
    os.makedirs(FIXTURES, exist_ok=True)
    env = dict(os.environ, PATH=bindir + os.pathsep + os.environ["PATH"], LC_ALL="C", TZ="UTC")
    for fx in fixtures:
        if not any("F/" + fx["made"] in a for a in needed_args): continue
        if os.path.exists(os.path.join(FIXTURES, fx["made"])): continue
        print(DIM(f"  fixture {fx['made']}: {fx['cmd']}"), flush=True)
        r = subprocess.run(expand(fx["cmd"], data), shell=True, cwd=FIXTURES, env=env)
        if r.returncode != 0 or not os.path.exists(os.path.join(FIXTURES, fx["made"])):
            die(f"fixture {fx['made']} failed")


# ------------------------------------------------------------------ running
def perf_ok():
    try:
        r = subprocess.run(["perf", "stat", "-x,", "-e", ",".join(EVENTS), "true"],
                           capture_output=True, text=True, timeout=30)
        return r.returncode == 0 and "not supported" not in r.stderr and "<not counted>" not in r.stderr
    except Exception:
        return False


def run_once(binary, args, cwd, use_perf, timeout):
    """Run `hmmsearch args` (argv[0] = 'hmmsearch', so outputs don't embed the path)."""
    bindir = os.path.join(cwd, ".bin")
    os.makedirs(bindir); os.symlink(binary, os.path.join(bindir, "hmmsearch"))
    env = dict(os.environ, PATH=bindir + os.pathsep + os.environ["PATH"], LC_ALL="C", TZ="UTC")
    env.pop("HMMER_NCPU", None)
    cmd = ["hmmsearch"] + args.split()
    if use_perf:
        cmd = ["perf", "stat", "-x,", "-o", os.path.join(cwd, ".perf"), "-e", ",".join(EVENTS), "--"] + cmd
    cmd = ["/usr/bin/time", "-q", "-o", os.path.join(cwd, ".rusage"), "-f", "%U %S %M %x"] + cmd
    t0 = time.perf_counter()
    with open(os.path.join(cwd, "stdout"), "wb") as so, open(os.path.join(cwd, "stderr"), "wb") as se:
        try:
            rc = subprocess.run(cmd, cwd=cwd, env=env, stdout=so, stderr=se, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = "timeout"
    r = dict(wall=time.perf_counter() - t0)
    try:
        u, s, m, x = open(os.path.join(cwd, ".rusage")).read().split()[-4:]
        r.update(user=float(u), sys=float(s), maxrss_mb=int(m) / 1024, rc=int(x) if rc != "timeout" else rc)
    except Exception:
        r.update(rc=rc)
    if use_perf and os.path.exists(os.path.join(cwd, ".perf")):
        for ln in open(os.path.join(cwd, ".perf")):
            f = ln.strip().split(",")
            if len(f) < 3 or ln.startswith("#"): continue
            key = {"task-clock": "task_clock", "instructions:u": "instructions", "cycles:u": "cycles"}.get(f[2])
            try:
                if key: r[key] = float(f[0]) / (1000.0 if key == "task_clock" else 1.0)
            except ValueError:
                pass
    shutil.rmtree(bindir)
    return r


def ws_hash(path):
    """Hash after the dev/bench normalization AND collapsing runs of blanks / trailing blanks.
    Two text files with equal ws_hash differ at most in column padding (accepted, see PLAN.md §6)."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        if b"\0" in f.read(8192): return None          # binary: no whitespace tolerance
        f.seek(0)
        for line in f:
            if line[:1] in (b"D", b"#", b"r"):
                t = line.decode("latin-1")
                for rx, rep in bench.TEXT_RULES: t = rx.sub(rep, t)
                line = t.encode("latin-1")
            h.update(WS.sub(b" ", line).rstrip() + b"\n")
    return h.hexdigest()


def hash_run(cwd):
    """{relpath: (exact normalized hash, whitespace-insensitive hash)}"""
    files = {}
    for dp, dns, fs in os.walk(cwd):
        for f in fs:
            p = os.path.join(dp, f)
            rel = os.path.relpath(p, cwd)
            if rel in OWN_FILES or rel.endswith(".norm") or os.path.islink(p): continue
            files[rel] = (bench.norm_hash(p), ws_hash(p))
    return files


def compare(h1, h2):
    """-> ('same'|'ws'|'diff', [differing files])"""
    keys = sorted(set(h1) | set(h2))
    bad = [k for k in keys if h1.get(k, (None,))[0] != h2.get(k, (None,))[0]]
    if not bad: return "same", []
    if all(k in h1 and k in h2 and h1[k][1] is not None and h1[k][1] == h2[k][1] for k in bad):
        return "ws", bad
    return "diff", bad


def normalized_copy(path):
    """Copy of a file with the volatile lines masked (as hashed), so diffs show real changes."""
    out = path + ".norm"
    with open(path, "rb") as f, open(out, "wb") as o:
        for line in f:
            if line[:1] in (b"D", b"#", b"r"):
                t = line.decode("latin-1")
                for rx, rep in bench.TEXT_RULES: t = rx.sub(rep, t)
                line = t.encode("latin-1")
            o.write(line)
    return out


def show_diff(da, db, files):
    for f in files[:3]:
        pa, pb = os.path.join(da, f), os.path.join(db, f)
        if not (os.path.exists(pa) and os.path.exists(pb)):
            print(f"      {f}: only in {'A' if os.path.exists(pa) else 'B'}"); continue
        na, nb = normalized_copy(pa), normalized_copy(pb)
        r = subprocess.run(["diff", "-u", "--label", "A/" + f, "--label", "B/" + f, na, nb],
                           capture_output=True, text=True, errors="replace")
        os.remove(na); os.remove(nb)
        print("\n".join("      " + (l if len(l) <= 160 else l[:157] + "...") for l in r.stdout.splitlines()[:20]))


# ------------------------------------------------------------------ stats
def med(xs): return statistics.median(xs) if xs else None
def spread(xs):
    m = med(xs)
    return (max(xs) - min(xs)) / m if xs and len(xs) > 1 and m else None


def verdict(metric, a, b):
    ma, mb = med(a), med(b)
    if not ma or mb is None: return None, None, ""
    ratio = mb / ma
    sp = [s for s in (spread(a), spread(b)) if s is not None]
    noise = max(sp) if sp else None
    thr = max(FLOOR[metric], noise or 0.0)
    v = "faster" if ratio < 1 - thr else "SLOWER" if ratio > 1 + thr else "~"
    return ratio, noise, v


def fmt(metric, x):
    if x is None: return "-"
    if metric in ("cycles", "instructions"): return f"{x / 1e9:.3f}G"
    return f"{x:.3f}s"


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("A"); ap.add_argument("B", nargs="?", default=os.path.join(ROOT, "build-dev", "src", "hmmsearch"))
    ap.add_argument("--set", default="cpu", help="cpu|scale|query|milestone|all, or comma list [cpu]")
    ap.add_argument("--workloads", help="regex on workload name (applied after --set)")
    ap.add_argument("--reps", type=int, default=5, help="timed runs per binary [5]")
    ap.add_argument("--warmup", type=int, default=1, help="untimed runs per binary first; its outputs are the A/B reference [1]")
    ap.add_argument("--check-only", action="store_true", help="outputs only: one run each, no timing (e.g. crosscheck builds)")
    ap.add_argument("--no-perf", action="store_true", help="don't use perf stat (wall/user/sys/RSS only)")
    ap.add_argument("--timeout", type=float, default=3600)
    ap.add_argument("--label")
    ap.add_argument("--keep", action="store_true", help="keep all run dirs (default: keep only mismatching ones)")
    ap.add_argument("--fail-on-slowdown", action="store_true")
    ap.add_argument("--list", action="store_true", help="list workloads and exit")
    a = ap.parse_args()

    ncpu = os.cpu_count()
    fixtures, wls = load_workloads(os.path.join(HERE, "workloads.tsv"), ncpu)
    sets = None if a.set == "all" else set(a.set.split(","))
    wls = [w for w in wls if (sets is None or sets & set(w["sets"])) and (not a.workloads or re.search(a.workloads, w["name"]))]
    if a.list:
        for w in wls: print(f"{w['name']:20s} {','.join(w['sets']):10s} {w['args']}")
        return 0
    if not wls: die("no workloads selected")
    if a.check_only: a.reps, a.warmup = 0, 1
    if a.warmup < 1 and a.reps < 1: die("need at least one run")

    data = bench.derived_dir()
    if not os.path.isdir(data): die(f"no benchmark data in {data}; run: dev/bench.sh prepare")
    A, B = resolve_bin(a.A), resolve_bin(a.B)
    da, db = describe_bin(A), describe_bin(B)
    use_perf = not a.no_perf and not a.check_only and perf_ok()
    make_fixtures(fixtures, [w["args"] for w in wls], data, os.path.dirname(A))

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    tag = f"{stamp}_{da.get('name', 'A')}_vs_{db.get('name', 'B')}" + (f"_{a.label}" if a.label else "")
    work = os.path.join(RESULTS, "work", tag)
    os.makedirs(work, exist_ok=True)
    load = os.getloadavg()[0]
    print(BLD(f"A: {da.get('name')}  {A}"))
    print(BLD(f"B: {db.get('name')}  {B}"))
    print(f"set={a.set} workloads={len(wls)} reps={a.reps} warmup={a.warmup} perf={'yes' if use_perf else 'no'} ncpu={ncpu}")
    if load > 1.0: print(YEL(f"WARNING: load average {load:.1f} > 1; timings will be noisy"))
    if A == B or da.get("hmmsearch") and da.get("hmmsearch") == db.get("hmmsearch"):
        print(DIM("(A and B are the same binary: this is an A/A noise run)"))

    results, failed, slower = [], False, False
    for w in wls:
        args = expand(w["args"], data)
        metric = "cycles" if w["sets"][0] == "cpu" and use_perf else "wall"
        runs = {"A": [], "B": []}
        ref, ref_dir, problems, notes = {}, {}, [], []
        nrun = a.warmup + a.reps
        print(f"\n{BLD(w['name'])}  {DIM(w['args'])}", flush=True)
        for i in range(nrun):
            for v in (("A", "B") if i % 2 == 0 else ("B", "A")):
                cwd = os.path.join(work, w["name"], f"{v}{i}")
                os.makedirs(cwd)
                r = run_once(A if v == "A" else B, args, cwd, use_perf and i >= a.warmup, a.timeout)
                h = hash_run(cwd)
                if r.get("rc") != 0:
                    problems.append(f"{v} run {i}: exit {r.get('rc')}"); failed = True
                if v not in ref:
                    ref[v], ref_dir[v] = h, cwd
                else:
                    kind, bad = compare(ref[v], h)
                    if kind == "ws":
                        notes.append(f"{v} run {i} vs its first run: whitespace-only ({', '.join(bad)})")
                    elif kind == "diff":
                        problems.append(f"{v} run {i} differs from its first run: {', '.join(bad)} (NONDETERMINISTIC)")
                        print(RED(f"  {v}{i} DIFFERS"))
                        show_diff(ref_dir[v], cwd, bad); failed = True
                        continue
                if i >= a.warmup: runs[v].append(r)
                if cwd not in ref_dir.values() and not a.keep: shutil.rmtree(cwd)
                print(DIM(f"  {v}{i} {r['wall']:.2f}s") + ("" if r.get("rc") == 0 else RED(f" rc={r.get('rc')}")), end="", flush=True)
        print()
        kind, bad = compare(ref["A"], ref["B"])
        same = kind != "diff"
        if kind == "ws":
            notes.append(f"A vs B: whitespace-only differences ({', '.join(bad)})")
        if kind == "diff":
            problems.append(f"OUTPUT DIFFERS A vs B: {', '.join(bad)}"); failed = True
            show_diff(ref_dir["A"], ref_dir["B"], bad)
        elif not a.keep and not problems:
            for d in ref_dir.values(): shutil.rmtree(d, ignore_errors=True)

        row = dict(name=w["name"], sets=w["sets"], args=w["args"], metric=metric, identical=same,
                   output=kind, notes=notes, problems=problems, files=sorted(ref["A"]), runs=runs, verdicts={})
        for m in ("cycles", "instructions", "task_clock", "wall"):
            xa = [r[m] for r in runs["A"] if m in r]; xb = [r[m] for r in runs["B"] if m in r]
            if xa and xb:
                ratio, noise, v = verdict(m, xa, xb)
                row["verdicts"][m] = dict(a=med(xa), b=med(xb), ratio=ratio, spread=noise, verdict=v)
        for m in ("user", "sys", "maxrss_mb"):
            xa = [r[m] for r in runs["A"] if m in r]; xb = [r[m] for r in runs["B"] if m in r]
            if xa and xb: row[m] = dict(a=med(xa), b=med(xb))
        if row["verdicts"].get(metric, {}).get("verdict") == "SLOWER": slower = True
        results.append(row)
        for p in problems: print(RED("  " + p))
        for n in notes: print(YEL("  note: " + n))

    # ---------------------------------------------------------------- report
    print("\n" + BLD("==> A/B summary") + f"  (ratio = B/A, <1 is faster; spread = (max-min)/median over reps)")
    print(f"  {'workload':20s} {'out':5s} {'metric':12s} {'A':>9s} {'B':>9s} {'ratio':>7s} {'spread':>7s}  verdict   {'instr B/A':>9s} {'wall B/A':>8s}  rss A→B MB")
    for row in results:
        out = (RED("DIFF ") if row["problems"] or not row["identical"]
               else YEL("ws   ") if row["output"] == "ws" or row["notes"] else GRN("same "))
        vm = row["verdicts"].get(row["metric"])
        ins = row["verdicts"].get("instructions", {}).get("ratio")
        wal = row["verdicts"].get("wall", {}).get("ratio")
        rss = row.get("maxrss_mb")
        if vm:
            col = GRN if vm["verdict"] == "faster" else RED if vm["verdict"] == "SLOWER" else (lambda s: s)
            print(f"  {row['name']:20s} {out} {row['metric']:12s} {fmt(row['metric'], vm['a']):>9s} {fmt(row['metric'], vm['b']):>9s} "
                  f"{vm['ratio']:7.4f} {('%6.2f%%' % (vm['spread'] * 100)) if vm['spread'] is not None else '      -'}  {col(vm['verdict'].ljust(8))}  "
                  f"{(f'{ins:.4f}' if ins else '-'):>9s} {(f'{wal:.3f}' if wal else '-'):>8s}  "
                  f"{('%.0f→%.0f' % (rss['a'], rss['b']) if rss else '-')}")
        else:
            print(f"  {row['name']:20s} {out} (no timing)")
    for m in ("cycles", "wall"):
        rs = [row["verdicts"][m]["ratio"] for row in results if m in row["verdicts"] and row["metric"] == m]
        if rs: print(f"  geometric mean {m} ratio: {bench.geomean(rs):.4f}  ({len(rs)} workloads)")

    os.makedirs(RESULTS, exist_ok=True)
    out = os.path.join(RESULTS, tag + ".json")
    json.dump(dict(date=stamp, A=da, B=db, set=a.set, reps=a.reps, warmup=a.warmup, perf=use_perf,
                   load_avg=load, machine=dict(ncpu=ncpu, host=os.uname().nodename), workloads=results),
              open(out, "w"), indent=1)
    if not any(os.scandir(work)) or (not a.keep and not failed):
        shutil.rmtree(work, ignore_errors=True)
    else:
        print(f"  run dirs kept in {work}")
    print(f"  results: {out}")
    if failed:
        print(RED(BLD("FAILED: outputs differ or a run failed"))); return 1
    ws = [r["name"] for r in results if r["output"] == "ws" or r["notes"]]
    print(GRN(BLD("outputs identical")) + (YEL(f" (whitespace-only differences: {', '.join(ws)})") if ws else "")
          + (RED("; SLOWER on some workloads") if slower else ""))
    return 3 if slower and a.fail_on_slowdown else 0


if __name__ == "__main__":
    sys.exit(main())
