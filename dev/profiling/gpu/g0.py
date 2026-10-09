#!/usr/bin/env python3
"""G0: the measured ceiling for a perfect filter GPU (dev/profiling/GPU.md, section G0).

Throwaway tool. Needs the patched binary from g0-replay.patch (HMMER_G0=record|ctl|A|B).
Per workload:
  1. record: patched binary at --cpu 0 writes every (query, target) filter outcome;
  2. rounds, interleaved and rotated: ref (unpatched), ctl (patched, lookups only),
     A (no SSV/MSV/VF), B (also no bias filter, no Forward parser for F3 failures);
  3. every variant's outputs must equal ref's (dev/perf/ab.py normalization).
Prints one row per workload: wall medians, the record's load time, worker utilisation of B.
"""
import argparse, json, os, re, shutil, statistics, sys, tempfile

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "perf"))
import ab, bench  # noqa

def count_targets(path):
    n = 0
    with open(path, "rb") as f:
        for ln in f:
            if ln[:1] == b">": n += 1
    return n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True); ap.add_argument("--patched", required=True)
    ap.add_argument("--file", default="workloads.tsv"); ap.add_argument("--workloads", required=True)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--rec-dir", default=os.path.expanduser("~/.cache/hmmer-bench/g0-records"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--variants", default="ref,ctl,A,B")
    ap.add_argument("--cpu-override", type=int, default=None, help="run the timed variants with --cpu N (the record always runs --cpu 0)")
    a = ap.parse_args()
    ref, pat = ab.resolve_bin(a.ref), os.path.realpath(a.patched)
    ncpu = os.cpu_count()
    fixtures, wls = ab.load_workloads(ab.workloads_file(a.file), ncpu)
    wls = [w for w in wls if re.search(a.workloads, w["name"])]
    data = bench.derived_dir()
    ab.make_fixtures(fixtures, [w["args"] for w in wls], data, os.path.dirname(ref))
    os.makedirs(a.rec_dir, exist_ok=True)
    variants = a.variants.split(",")
    results = []
    for w in wls:
        args = ab.expand(w["args"], data)
        tfile = args.split()[-1]
        rec = os.path.join(a.rec_dir, w["name"] + ".g0")
        env_rec = {"HMMER_G0": "record", "HMMER_G0_FILE": rec, "HMMER_G0_NT": str(count_targets(tfile))}
        cmd0 = re.sub(r"--cpu \S+", "--cpu 0", args) if "--cpu" in args else args
        if a.cpu_override is not None: args = re.sub(r"--cpu \S+", f"--cpu {a.cpu_override}", args)
        print(f"== {w['name']}: record ({env_rec['HMMER_G0_NT']} targets)", flush=True)
        tmp = tempfile.mkdtemp(prefix="g0-", dir=os.environ.get("G0_TMP", "/tmp"))
        if not os.path.exists(rec):
            d = os.path.join(tmp, "rec"); os.makedirs(d)
            r = ab.run_once(pat, cmd0, d, False, 7200, extra_env=env_rec)
            print(f"   record wall {r['wall']:.1f} s rc {r.get('rc')}", flush=True)
            if r.get("rc") != 0: print(open(os.path.join(d, "stderr"), errors="replace").read()[-500:]); continue
            rec_hash = ab.hash_run(d)
        else:
            rec_hash = None
        bins = {"ref": (ref, {}), "ctl": (pat, {"HMMER_G0": "ctl", "HMMER_G0_FILE": rec}),
                "A": (pat, {"HMMER_G0": "A", "HMMER_G0_FILE": rec}), "B": (pat, {"HMMER_G0": "B", "HMMER_G0_FILE": rec})}
        times = {v: [] for v in variants}; info = {v: [] for v in variants}
        status = {v: "same" for v in variants}; refhash = None
        for rep in range(-1, a.reps):                                  # round -1: untimed reference
            order = variants[(rep % len(variants)):] + variants[:(rep % len(variants))] if rep >= 0 else variants
            for v in order:
                d = os.path.join(tmp, f"{v}{rep}"); os.makedirs(d)
                b, env = bins[v]
                r = ab.run_once(b, args, d, False, 7200, extra_env=env)
                if r.get("rc") != 0: status[v] = f"FAILED rc={r.get('rc')}"; continue
                se = open(os.path.join(d, "stderr"), errors="replace").read()
                os.remove(os.path.join(d, "stderr"))   # holds the patch's G0 lines; real stderr of ref was empty
                h = ab.hash_run(d)
                if v == "ref" and refhash is None: refhash = h
                elif refhash is not None:
                    c, bad = ab.compare(refhash, h)
                    if c != "same" and status[v] == "same": status[v] = c + " " + ",".join(bad[:3])
                m = re.search(r"G0: loaded .* in ([\d.]+) s", se)
                workers = [tuple(map(float, x)) for x in re.findall(r"G0W \d+ cpu ([\d.]+) search ([\d.]+) wait ([\d.]+)", se)]
                if rep >= 0:
                    times[v].append(r["wall"])
                    info[v].append(dict(user=r.get("user"), sys=r.get("sys"), rss=r.get("maxrss_mb"), load=float(m.group(1)) if m else 0.0, workers=workers))
                shutil.rmtree(d)
        row = dict(name=w["name"], args=w["args"], status=status, wall={v: times[v] for v in variants}, info=info)
        results.append(row)
        med = {v: statistics.median(times[v]) for v in variants if times[v]}
        load = statistics.median([i["load"] for i in info["B"]]) if info.get("B") else 0
        s = f"   {w['name']}: " + "  ".join(f"{v} {med[v]:.2f}s" for v in med)
        if "B" in med and "ref" in med: s += f"  | ceiling B: {med['ref']/med['B']:.2f}x ({med['ref']/(med['B']-load):.2f}x without the {load:.2f}s load)"
        if "A" in med and "ref" in med: s += f"  A: {med['ref']/med['A']:.2f}x"
        print(s); print("   outputs vs ref:", status, flush=True)
        shutil.rmtree(tmp, ignore_errors=True)
    if a.out: json.dump(results, open(a.out, "w"), indent=1)

main()
