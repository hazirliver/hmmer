#!/usr/bin/env python3
"""HMMER speed benchmark + output-regression harness.

Subcommands:
  fetch      download pinned public inputs (dev/bench/sources.tsv), verify md5
  prepare    derive the benchmark datasets (dev/bench/prepare.sh), verify sha256
  run        time every case of a tier against a build; check output hashes;
             compare to a baseline run
  compare    compare two result JSON files
See dev/bench/README.md.
"""
import argparse, datetime, hashlib, json, os, platform, re, shutil, signal
import socket, statistics, subprocess, sys, threading, time

ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", ".."))
BENCH = os.path.join(ROOT, "dev", "bench")
DATA = os.environ.get("HMMER_BENCH_DATA", os.path.expanduser("~/.cache/hmmer-bench"))
RESULTS = os.environ.get("HMMER_BENCH_RESULTS", os.path.join(ROOT, "bench-results"))
TIERS = ["S", "M", "L"]          # cumulative: M includes S, L includes M

TTY = sys.stdout.isatty()
def c(code, s): return f"\033[{code}m{s}\033[0m" if TTY else s
RED, GRN, YEL, BLD, DIM = (lambda s: c(31, s)), (lambda s: c(32, s)), (lambda s: c(33, s)), (lambda s: c(1, s)), (lambda s: c(2, s))

def die(msg):
    print(RED("ERROR: " + msg), file=sys.stderr); sys.exit(1)

def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True, **kw).stdout.strip()

def md5(path, algo="md5"):
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def read_table(path, ncols):
    """Lines '<f1> <f2> ... <rest-of-line>' ; skip comments/blank."""
    rows = []
    for ln in open(path):
        if not ln.strip() or ln.lstrip().startswith("#"): continue
        parts = ln.rstrip("\n").split(None, ncols - 1)
        if len(parts) < ncols: die(f"{path}: malformed line: {ln!r}")
        rows.append(parts)
    return rows

# ---------------------------------------------------------------- fetch
def cmd_fetch(a):
    raw = os.path.join(DATA, "raw"); os.makedirs(raw, exist_ok=True)
    ok = True
    for name, want, url in read_table(os.path.join(BENCH, "sources.tsv"), 3):
        dst = os.path.join(raw, name + ".gz")
        if not os.path.exists(dst):
            print(f"downloading {name}  <- {url}")
            subprocess.run(["curl", "-fL", "--retry", "3", "-o", dst + ".part", url], check=True)
            os.rename(dst + ".part", dst)
        have = md5(dst)
        if want == "-":
            print(YEL(f"  {name}: no pinned md5; have {have}"))
        elif have != want:
            print(RED(f"  {name}: md5 mismatch (have {have}, want {want})")); ok = False
        else:
            print(f"  {name}: ok")
    if not ok: die("checksum mismatch; delete the bad file(s) in " + raw + " and re-run")

# ---------------------------------------------------------------- prepare
def data_version():
    h = hashlib.sha256()
    for f in ("sources.tsv", "prepare.sh"):
        h.update(open(os.path.join(BENCH, f), "rb").read())
    return h.hexdigest()[:12]

def derived_dir():
    return os.path.join(DATA, "derived-" + data_version())

def verify_manifest(d, quiet=False):
    """Compare derived data against the committed dev/bench/data.sha256."""
    golden = os.path.join(BENCH, "data.sha256")
    if not os.path.exists(golden):
        print(YEL("  no dev/bench/data.sha256 yet (run 'prepare --pin' to create it)")); return True
    want = dict(reversed(l.split()) for l in open(golden) if l.strip())
    bad = [f for f, h in want.items() if not os.path.exists(os.path.join(d, f)) or (h != md5(os.path.join(d, f), "sha256"))]
    if bad:
        print(RED(f"  {len(bad)} derived file(s) differ from dev/bench/data.sha256: " + ", ".join(bad[:8])))
        return False
    if not quiet: print(f"  derived data verified: {len(want)} files match dev/bench/data.sha256")
    return True

def cmd_prepare(a):
    cmd_fetch(a)
    d = derived_dir()
    if os.path.exists(os.path.join(d, ".complete")) and not a.force:
        print(f"derived data present: {d}")
    else:
        build = os.path.realpath(a.build_dir)
        esl = os.path.join(build, "easel", "miniapps")
        if not os.path.exists(os.path.join(esl, "esl-reformat")): die(f"no Easel miniapps in {esl}; build first")
        if os.path.exists(d): shutil.rmtree(d)
        os.makedirs(d)
        t0 = time.time()
        print(f"preparing derived data in {d} (takes a few minutes)...")
        r = subprocess.run(["bash", os.path.join(BENCH, "prepare.sh"), os.path.join(DATA, "raw"), d, esl])
        if r.returncode: die("prepare.sh failed")
        files = sorted(os.path.relpath(os.path.join(dp, f), d) for dp, _, fs in os.walk(d) for f in fs)
        with open(os.path.join(d, "MANIFEST.sha256"), "w") as m:
            for f in files: m.write(f"{md5(os.path.join(d, f), 'sha256')}  {f}\n")
        open(os.path.join(d, ".complete"), "w").close()
        print(f"prepared in {time.time()-t0:.0f}s")
    if a.pin:
        shutil.copy(os.path.join(d, "MANIFEST.sha256"), os.path.join(BENCH, "data.sha256"))
        print("pinned derived-data checksums into dev/bench/data.sha256")
    elif not verify_manifest(d):
        die("derived data does not match the pinned checksums")
    sizes = sh(f"du -sh {d} | cut -f1")
    print(f"dataset: {d} ({sizes})")

# ---------------------------------------------------------------- normalization / hashing
TEXT_RULES = [(re.compile(p, re.M), r) for p, r in [
    (r"^(DATE  ).{24}$", r"\1Thu Jan  1 00:00:00 1970"),
    (r"^# Date:.*$", "# Date:            [normalized]"),
    (r"^# Current dir:.*$", "# Current dir:     [normalized]"),
    (r"^# CPU time:.*$", "# CPU time: [normalized]"),
    (r"^# Mc/sec:.*$", "# Mc/sec: [normalized]"),
    (r"^run time: .*$", "run time: [normalized]"),
    # cases use --cpu $NCPU: keep hashes portable across machines with different core counts
    (r"^# number of worker threads:.*$", "# number of worker threads: [normalized]"),
]]
CTIME = re.compile(rb"(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun) (?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [ \d]\d \d\d:\d\d:\d\d \d{4}")

def norm_hash(path):
    """sha256 of a file after masking volatile content; streams (outputs can be GBs)."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        binary = b"\0" in f.read(8192)
        f.seek(0)
        if binary:
            carry = b""
            for chunk in iter(lambda: f.read(1 << 24), b""):
                buf = CTIME.sub(b"Thu Jan  1 00:00:00 1970", carry + chunk)
                h.update(buf[:-23]); carry = buf[-23:]
            h.update(carry)
        else:
            for line in f:
                if line[:1] in (b"D", b"#", b"r"):   # all TEXT_RULES anchor on these
                    s = line.decode("latin-1")
                    for rx, rep in TEXT_RULES: s = rx.sub(rep, s)
                    line = s.encode("latin-1")
                h.update(line)
    return h.hexdigest()

def hash_dir(d):
    files = {}
    for dp, dns, fs in os.walk(d):
        dns[:] = [x for x in dns if not os.path.islink(os.path.join(dp, x))]
        for f in fs:
            p = os.path.join(dp, f)
            if os.path.islink(p): continue
            files[os.path.relpath(p, d)] = norm_hash(p)
    total = hashlib.sha256("".join(f"{k}\t{v}\n" for k, v in sorted(files.items())).encode()).hexdigest()
    return total, files

# ---------------------------------------------------------------- run
def tier_ok(case_tier, tier):
    return TIERS.index(case_tier) <= TIERS.index(tier)

def load_cases(tier, regex, exact=False):
    fx = [dict(name=n, tier=t, cmd=cmd) for n, t, cmd in read_table(os.path.join(BENCH, "fixtures.tsv"), 3) if tier_ok(t, tier)]
    cases = []
    for n, t, g, cmd in read_table(os.path.join(BENCH, "cases.tsv"), 4):
        if not tier_ok(t, tier) or (exact and t != tier): continue
        if regex and not re.search(regex, n) and not re.search(regex, g.split(":")[0]): continue
        pre, _, timed = cmd.rpartition(";;")
        flags = g.split(":")[1:]
        cases.append(dict(name=n, tier=t, group=g.split(":")[0], hash="nohash" not in flags, kernel="kernel" in flags,
                          pre=pre.strip(), cmd=timed.strip()))
    return fx, cases

def wait_timed(p, timeout, kill):
    """Wait for Popen <p>; if it runs longer than <timeout> s, call kill(). Returns
    True if it timed out. The wait blocks, so the caller's clock stops when the
    process exits: Popen.wait(timeout=...) instead polls with sleeps that grow to
    50 ms, which rounded every measured wall time up to the next poll."""
    fired = threading.Event()
    def on_timeout():
        fired.set(); kill()
    timer = threading.Timer(timeout, on_timeout)
    timer.start()
    try:
        p.wait()
    finally:
        timer.cancel()
    return fired.is_set()

def run_timed(cmd, cwd, env, timeout):
    """Run cmd under bash; return wall, user, sys, maxrss_kb, rc.
    Rusage comes from GNU time, a small intermediate process: Linux carries a
    process's peak RSS across exec, so measuring from python would report
    python's own RSS as a floor."""
    tf = os.path.join(cwd, ".rusage")
    t0 = time.perf_counter()
    p = subprocess.Popen(["/usr/bin/time", "-q", "-o", tf, "-f", "%U %S %M %x", "bash", "-o", "pipefail", "-c", cmd], cwd=cwd, env=env,
                         start_new_session=True, stdout=open(os.path.join(cwd, "stdout"), "wb"),
                         stderr=open(os.path.join(cwd, "stderr"), "wb"))
    rc = "timeout" if wait_timed(p, timeout, lambda: os.killpg(p.pid, signal.SIGKILL)) else None
    wall = time.perf_counter() - t0
    try:
        u, s_, m, x = open(tf).read().split()[-4:]
        r = dict(wall=wall, user=float(u), sys=float(s_), maxrss_kb=int(m), rc=rc or int(x))
    except Exception:
        r = dict(wall=wall, user=0.0, sys=0.0, maxrss_kb=0, rc=rc or p.returncode)
    if os.path.exists(tf): os.remove(tf)
    return r

KTIME = re.compile(r"^# CPU time:.*Elapsed: (\d+):(\d+):([\d.]+)", re.M)
KMCS = re.compile(r"^# ([\d.]+) Mc/s", re.M)

def kernel_time(cdir, r):
    """For kernel drivers: replace process wall time by the driver's own timed-loop time."""
    out = open(os.path.join(cdir, "stdout"), errors="replace").read()
    m = KTIME.search(out)
    if not m: r["rc"] = r["rc"] or "no-kernel-time"; return r
    r["proc_wall"] = r["wall"]
    r["wall"] = int(m.group(1)) * 3600 + int(m.group(2)) * 60 + float(m.group(3))
    mc = KMCS.search(out)
    if mc: r["mcs"] = float(mc.group(1))
    return r

def src_dir(build):
    """Top source dir the build was configured from."""
    m = re.search(r"^srcdir\s*=\s*(\S+)", open(os.path.join(build, "Makefile")).read(), re.M)
    return os.path.realpath(os.path.join(build, m.group(1) if m else "."))

def fresh_dir(d, links):
    if os.path.exists(d): subprocess.run(["chmod", "-R", "u+w", d]); shutil.rmtree(d)
    os.makedirs(d)
    for k, v in links.items(): os.symlink(v, os.path.join(d, k))

def machine_info(build, src):
    def q(cmd):
        try: return sh(cmd)
        except Exception: return ""
    cfg = os.path.join(build, "config.log")
    conf = ""
    if os.path.exists(cfg):
        m = re.search(r"^\s+\$ (.*configure.*)$", open(cfg, errors="replace").read(), re.M)
        conf = m.group(1) if m else ""
    cflags = q(f"grep -E '^CFLAGS|^SIMD_CFLAGS|^CC ' {build}/Makefile {build}/src/Makefile 2>/dev/null | sort -u | tr '\\n' ' '")
    return dict(
        host=socket.gethostname(), cpu=q("lscpu | sed -n 's/^Model name:\\s*//p'"), nproc=os.cpu_count(),
        mem_gb=round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**30), kernel=platform.release(),
        compiler=q(f"$(sed -n 's/^CC *= *//p' {build}/Makefile) --version | head -1"), configure=conf, flags=cflags,
        src=src, git_sha=q(f"git -C {src} rev-parse HEAD"), git_branch=q(f"git -C {src} rev-parse --abbrev-ref HEAD"),
        git_dirty=bool(q(f"git -C {src} status --porcelain --untracked-files=no")),
        easel_sha=q(f"git -C {src}/easel rev-parse HEAD"), load_avg=os.getloadavg()[0])

def impl_dir(build):
    for i in ("impl_sse", "impl_neon", "impl_vmx"):
        if os.path.isdir(os.path.join(build, "src", i)) and os.path.exists(os.path.join(build, "src", i, "msvfilter_benchmark")):
            return os.path.join(build, "src", i)
    return os.path.join(build, "src", "impl")

def cmd_run(a):
    build = os.path.realpath(a.build_dir)
    if not os.path.exists(os.path.join(build, "src", "hmmsearch")): die(f"no HMMER build in {build}")
    d = derived_dir()
    if not os.path.exists(os.path.join(d, ".complete")): die("benchmark data not prepared: run dev/bench.sh prepare")
    if not a.skip_verify and not verify_manifest(d, quiet=True): die("derived data does not match dev/bench/data.sha256")
    fixtures, cases = load_cases(a.tier, a.cases, a.exact_tier)
    if not cases: die("no cases selected")
    work = os.path.realpath(a.work or os.path.join(build, "bench-work"))
    if os.path.exists(work): subprocess.run(["chmod", "-R", "u+w", work]); shutil.rmtree(work)
    os.makedirs(work)
    env = dict(os.environ, LC_ALL="C", TZ="UTC", ROOT=src_dir(build), BUILD=build, IMPL=impl_dir(build), NCPU=str(os.cpu_count()),
               PATH=f"{build}/src:{build}/easel/miniapps:{os.environ['PATH']}")
    env.pop("HMMER_NCPU", None)
    golden = load_golden()
    info = machine_info(build, src_dir(build))
    res = dict(meta=dict(info, tier=a.tier, reps=a.reps, label=a.label or "", build_dir=build, data=os.path.basename(d),
                         started=datetime.datetime.now().isoformat(timespec="seconds"), argv=sys.argv[1:]), fixtures={}, cases={})
    if info["load_avg"] > 1.0: print(YEL(f"note: load average is {info['load_avg']:.1f}; timings may be noisy"))
    print(BLD(f"==> bench tier {a.tier}: {len(cases)} cases, reps={a.reps}, build={build}"))

    # fixtures: untimed for comparisons but recorded
    fxd = os.path.join(work, "fx"); os.makedirs(fxd); os.symlink(d, os.path.join(fxd, "D"))
    for f in fixtures:
        print(f"  fixture {f['name']:<28}", end="", flush=True)
        r = run_timed(f["cmd"], fxd, env, a.timeout)
        res["fixtures"][f["name"]] = r
        print(f"{r['wall']:8.2f}s")
        if r["rc"] != 0:
            print(open(os.path.join(fxd, "stderr")).read()[-2000:]); die(f"fixture {f['name']} failed")
    os.remove(os.path.join(fxd, "stdout")); os.remove(os.path.join(fxd, "stderr"))
    subprocess.run(["chmod", "-R", "a-w", fxd])

    bad = 0
    for cs in cases:
        cdir = os.path.join(work, "cases", cs["name"])
        runs, hashes, files = [], [], None
        for rep in range(a.reps):
            fresh_dir(cdir, {"D": d, "fx": fxd})
            if cs["pre"]:
                p = subprocess.run(["bash", "-o", "pipefail", "-c", cs["pre"]], cwd=cdir, env=env, capture_output=True, text=True)
                if p.returncode: print(p.stderr[-2000:]); die(f"{cs['name']}: setup failed")
            r = run_timed(cs["cmd"], cdir, env, a.timeout)
            if cs["kernel"] and r["rc"] == 0: r = kernel_time(cdir, r)
            runs.append(r)
            if cs["hash"]:
                h, fl = hash_dir(cdir); hashes.append(h); files = files or fl
            if r["rc"] != 0 or (rep == 0 and r["wall"] > a.long and a.reps > 1):
                break   # failures and long cases: don't repeat
        walls = [r["wall"] for r in runs]
        e = dict(group=cs["group"], tier=cs["tier"], cmd=cs["cmd"], walls=walls, median=statistics.median(walls),
                 min=min(walls), user=statistics.median(r["user"] for r in runs), sys=statistics.median(r["sys"] for r in runs),
                 maxrss_kb=max(r["maxrss_kb"] for r in runs), rc=runs[0]["rc"])
        if "mcs" in runs[0]: e["mcs"] = statistics.median(r["mcs"] for r in runs)
        status = "nohash"
        if e["rc"] != 0:
            status = "FAILED"; tail = open(os.path.join(cdir, "stderr"), errors="replace").read()[-1500:]; e["stderr"] = tail
        elif cs["hash"]:
            e["hash"], e["files"] = hashes[0], files
            if len(set(hashes)) > 1: status = "NONDETERMINISTIC"
            elif cs["name"] not in golden: status = "new"
            elif golden[cs["name"]]["hash"] == hashes[0]: status = "ok"
            else:
                status = "MISMATCH"
                gf = golden[cs["name"]]["files"]
                e["diff_files"] = sorted(k for k in set(gf) | set(files) if gf.get(k) != files.get(k))
        e["output"] = status
        if status in ("FAILED", "MISMATCH", "NONDETERMINISTIC"): bad += 1
        res["cases"][cs["name"]] = e
        col = {"ok": GRN, "nohash": DIM, "new": YEL}.get(status, RED)
        cpu = (e["user"] + e["sys"]) / e["median"] if e["median"] > 0 and not cs["kernel"] else 1
        print(f"  {cs['group']:<11} {cs['name']:<34} {e['median']:9.3f}s  (min {e['min']:.3f}, n={len(walls)})  "
              f"cpu {cpu:4.1f}x  rss {e['maxrss_kb']/1024:7.0f}M  out:{col(status)}"
              + (f"  [{', '.join(e.get('diff_files', [])[:3])}]" if status == "MISMATCH" else "")
              + (f"\n{e['stderr']}" if status == "FAILED" else ""))
        if not a.keep: subprocess.run(["chmod", "-R", "u+w", cdir]); shutil.rmtree(cdir)
    res["meta"]["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    res["meta"]["total_wall"] = sum(e["median"] for e in res["cases"].values())

    if a.update_hashes:
        for n, e in res["cases"].items():
            if "hash" in e and e["output"] != "NONDETERMINISTIC": golden[n] = dict(hash=e["hash"], files=e["files"])
        save_golden(golden)
        print(YEL(f"updated dev/bench/golden.json for {sum('hash' in e for e in res['cases'].values())} cases"))

    os.makedirs(RESULTS, exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    out = os.path.join(RESULTS, f"{stamp}_{info['git_sha'][:8]}{'-dirty' if info['git_dirty'] else ''}_{a.tier}"
                       + (f"_{a.label}" if a.label else "") + ".json")
    json.dump(res, open(out, "w"), indent=1)
    print(f"\nresults: {out}")

    slow = 0
    base = resolve_baseline(a.baseline, a.tier, out)
    if base:
        slow = compare(json.load(open(base)), res, base, out, a.tolerance, a.min_time)
    else:
        print(DIM("no baseline to compare against (use --save-baseline to create one)"))
    if a.save_baseline:
        dst = os.path.join(RESULTS, f"baseline-{a.tier}.json"); shutil.copy(out, dst); print(f"saved as baseline: {dst}")

    print()
    if bad: print(RED(BLD(f"OUTPUT CHECK FAILED: {bad} case(s) failed / differ from golden / nondeterministic")))
    else: print(GRN(BLD("outputs: all hashed cases match golden")) if not a.update_hashes else "")
    if slow and a.fail_on_slowdown: print(RED(BLD(f"PERFORMANCE: {slow} case(s) slower than tolerance")))
    sys.exit(1 if bad or (slow and a.fail_on_slowdown) else 0)

def load_golden():
    p = os.path.join(BENCH, "golden.json")
    return json.load(open(p)) if os.path.exists(p) else {}

def save_golden(g):
    json.dump(dict(sorted(g.items())), open(os.path.join(BENCH, "golden.json"), "w"), indent=1, sort_keys=True)

def resolve_baseline(spec, tier, current):
    if spec == "none": return None
    if spec and os.path.exists(spec): return spec
    if spec and spec != "auto":
        cands = sorted(f for f in os.listdir(RESULTS) if spec in f and f.endswith(".json"))
        if not cands: die(f"no result matching baseline '{spec}' in {RESULTS}")
        return os.path.join(RESULTS, cands[-1])
    for t in TIERS[TIERS.index(tier):]:   # a bigger-tier baseline also covers this tier's cases
        p = os.path.join(RESULTS, f"baseline-{t}.json")
        if os.path.exists(p) and os.path.realpath(p) != os.path.realpath(current): return p
    return None

# ---------------------------------------------------------------- compare
def geomean(xs):
    xs = [x for x in xs if x > 0]
    return (statistics.geometric_mean(xs) if xs else float("nan"))

def compare(base, new, bpath, npath, tol, min_time):
    bm, nm = base["meta"], new["meta"]
    print(BLD(f"\n==> comparison  (ratio = new/base wall time; <1 is faster; tolerance ±{tol:.0%})"))
    print(f"  base: {os.path.basename(bpath)}  [{bm.get('git_sha','')[:8]} {bm.get('git_branch','')} {bm.get('label','')}]")
    print(f"  new:  {os.path.basename(npath)}  [{nm.get('git_sha','')[:8]} {nm.get('git_branch','')} {nm.get('label','')}]")
    for k in ("host", "cpu", "compiler", "flags"):
        if bm.get(k) != nm.get(k): print(YEL(f"  note: {k} differs: {bm.get(k)!r} vs {nm.get(k)!r}"))
    groups, slower, rows = {}, 0, []
    for n, e in new["cases"].items():
        b = base["cases"].get(n)
        if not b or b.get("rc") != 0 or e.get("rc") != 0: continue
        r = e["median"] / b["median"] if b["median"] > 0 else float("nan")
        verdict = ""
        if max(b["median"], e["median"]) >= min_time:
            # slower only if even the best new run is beyond tolerance of the base median
            if e["min"] > b["median"] * (1 + tol): verdict = "SLOWER"; slower += 1
            elif e["median"] < b["median"] * (1 - tol): verdict = "faster"
        groups.setdefault(e["group"], []).append(r)
        rows.append((e["group"], n, b["median"], e["median"], r, verdict, b["maxrss_kb"], e["maxrss_kb"]))
    print(f"  {'group':<11} {'case':<34} {'base':>9} {'new':>9} {'ratio':>7}         {'rss base→new (MB)':>18}")
    for g, n, bt, nt, r, v, br, nr in rows:
        col = RED if v == "SLOWER" else GRN if v == "faster" else (lambda s: s)
        rss = f"{br/1024:.0f}→{nr/1024:.0f}" + (YEL(" !") if nr > br * 1.2 and nr - br > 10240 else "")
        print(f"  {g:<11} {n:<34} {bt:9.3f} {nt:9.3f} {col(f'{r:7.3f}')} {col(f'{v:<7}')} {rss:>18}")
    print(BLD("  per-group geometric mean ratio:"))
    for g, rs in groups.items(): print(f"    {g:<11} {geomean(rs):6.3f}  ({len(rs)} cases)")
    allr = [r for rs in groups.values() for r in rs]
    print(BLD(f"    {'OVERALL':<11} {geomean(allr):6.3f}  ({len(allr)} cases)"))
    missing = sorted(set(base["cases"]) - set(new["cases"]))
    if missing and len(missing) < 15: print(DIM(f"  (not in new run: {', '.join(missing)})"))
    return slower

def cmd_compare(a):
    s = compare(json.load(open(a.base)), json.load(open(a.new)), a.base, a.new, a.tolerance, a.min_time)
    sys.exit(1 if s and a.fail_on_slowdown else 0)

# ---------------------------------------------------------------- main
def main():
    sys.stdout.reconfigure(line_buffering=True)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="sub", required=True)
    p = sp.add_parser("fetch"); p.set_defaults(fn=cmd_fetch)
    p = sp.add_parser("prepare"); p.set_defaults(fn=cmd_prepare)
    p.add_argument("--build-dir", default=os.path.join(ROOT, "build-dev"))
    p.add_argument("--force", action="store_true", help="re-derive even if present")
    p.add_argument("--pin", action="store_true", help="write derived checksums to dev/bench/data.sha256")
    def cmp_opts(p):
        p.add_argument("--tolerance", type=float, default=0.10, help="relative slowdown/speedup threshold [0.10]")
        p.add_argument("--min-time", type=float, default=0.5, help="ignore cases faster than this (s) for verdicts [0.5]")
        p.add_argument("--fail-on-slowdown", action="store_true")
    p = sp.add_parser("run"); p.set_defaults(fn=cmd_run)
    p.add_argument("tier", nargs="?", default="S", choices=TIERS)
    p.add_argument("--build-dir", default=os.path.join(ROOT, "build-dev"))
    p.add_argument("--reps", type=int, default=3)
    p.add_argument("--long", type=float, default=60.0, help="cases slower than this (s) run only once [60]")
    p.add_argument("--timeout", type=float, default=3 * 3600)
    p.add_argument("--cases", help="regex on case name or group")
    p.add_argument("--exact-tier", action="store_true", help="only cases of exactly this tier (not the smaller tiers)")
    p.add_argument("--label")
    p.add_argument("--baseline", default="auto", help="'auto' (baseline-<tier>.json), 'none', a path, or a substring of a results file name")
    p.add_argument("--save-baseline", action="store_true")
    p.add_argument("--update-hashes", action="store_true", help="record current output hashes as golden")
    p.add_argument("--keep", action="store_true", help="keep case work dirs")
    p.add_argument("--work")
    p.add_argument("--skip-verify", action="store_true")
    cmp_opts(p)
    p = sp.add_parser("compare"); p.set_defaults(fn=cmd_compare)
    p.add_argument("base"); p.add_argument("new"); cmp_opts(p)
    a = ap.parse_args(); a.fn(a)

if __name__ == "__main__":
    main()
