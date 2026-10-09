#!/usr/bin/env python3
"""Charts (SVG, light + dark) and a Markdown table for the Pfam-A annotation benchmark.

  dev/perf/plot_annot.py annot.json --outdir dev/profiling/results

annot.json is a dev/perf/series.py result over dev/perf/annot.tsv (first binary = upstream, last =
the build measured against it). Group sizes and length bounds come from dev/bench/uniprot/manifest.tsv.

Speed is throughput: sequences searched per second of wall time (250 per run, against all of
Pfam-A). Writes, per mode (light, dark):
  annot-pressed-<mode>.svg  one grouped bar chart per protein set: length groups L1-L4, upstream
                            vs this branch, Pfam-A pressed with hmmpress
  annot-text-<mode>.svg     the same from the text Pfam-A file
and prints the Markdown table (sequences per second per group, per set and overall, and peak memory).
"""
import argparse, json, math, os, statistics, sys
from html import escape

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
from plot_series import MODES, text, svg_doc, legend, fmt_x, nice_step   # noqa: E402

MANIFEST = os.path.join(HERE, "..", "bench", "uniprot", "manifest.tsv")
SETS = [("bacteria", "Bacteria"), ("phage", "Phages"), ("fungi", "Fungi"),
        ("virus", "Eukaryotic viruses"), ("eukaryota", "Non-human eukaryotes"), ("human", "Human")]
FORMS = [("", "pressed", "Pfam-A pressed (hmmpress)"), ("_txt", "text", "Pfam-A text file")]
GROUPS = 4
NUM = 'font-variant-numeric="tabular-nums"'


def medians(res):
    return {w["name"]: {l: statistics.median(r["wall"] for r in runs) for l, runs in w["runs"].items() if runs}
            for w in res["workloads"]}


def peak_rss(res):
    return {w["name"]: {l: max(r["maxrss_mb"] for r in runs) for l, runs in w["runs"].items() if runs}
            for w in res["workloads"]}


def groups():
    """{(set, g): dict(n, lo, hi)} from the manifest."""
    out = {}
    for ln in open(MANIFEST):
        if ln.startswith("#") or not ln.strip(): continue
        f = ln.rstrip("\n").split("\t")
        k, L = (f[0], int(f[1])), int(f[6])
        d = out.setdefault(k, dict(n=0, lo=L, hi=L))
        d["n"] += 1; d["lo"] = min(d["lo"], L); d["hi"] = max(d["hi"], L)
    return out


def fmt_r(v):  return f"{v:.1f}" if v < 10 else f"{v:.0f}"


def vbar(x, y0, w, y1, color, tip):
    """bar from the baseline y0 up to y1, 3px rounded data end"""
    h = y0 - y1
    r = min(3, w / 2, h / 2) if h > 0 else 0
    d = (f"M{x:.1f},{y0:.1f} v{-(h - r):.1f} a{r},{r} 0 0 1 {r},{-r} h{w - 2 * r:.1f} "
         f"a{r},{r} 0 0 1 {r},{r} v{h - r:.1f} z")
    return f'<path d="{d}" fill="{color}"><title>{escape(tip)}</title></path>'


def rates(med, grp, sfx, first, last):
    """{(set, g): (upstream seq/s, now seq/s)}, and per set and "all": [n, upstream s, now s]"""
    per, tot = {}, {}
    for s, _ in SETS:
        for g in range(1, GROUPS + 1):
            m = med.get(f"A_{s}_L{g}{sfx}")
            if not m or first not in m or last not in m: continue
            n = grp[(s, g)]["n"]
            per[(s, g)] = (n / m[first], n / m[last])
            for k in (s, "all"):
                t = tot.setdefault(k, [0, 0.0, 0.0])
                t[0] += n; t[1] += m[first]; t[2] += m[last]
    return per, tot


def chart(med, grp, sfx, title, first, last, c):
    per, tot = rates(med, grp, sfx, first, last)
    up_color, now_color = c["muted"], c["s1"]
    cols, gx, gy = 3, 30, 34
    W = 900
    pw = (W - 20 - (cols - 1) * gx) / cols
    ph = 250
    top = 56
    rows = math.ceil(len(SETS) / cols)
    H = top + rows * ph + (rows - 1) * gy + 30
    vmax = max(v for pair in per.values() for v in pair)
    step = nice_step(vmax / 4)
    ymax = step * math.ceil(vmax * 1.08 / step)
    body = legend(10, 22, c, [(up_color, "upstream HMMER 3.4"), (now_color, "this branch")])
    body.append(text(W - 10, 22, f"{title}: sequences per second, 16 threads", c["ink2"], 12, "end"))
    bw, bgap = 19, 3
    for idx, (s, name) in enumerate(SETS):
        px = 10 + (idx % cols) * (pw + gx)
        py = top + (idx // cols) * (ph + gy)
        l, r, t, b = px + 34, px + pw, py + 42, py + ph - 52
        sy = lambda v: b - (b - t) * v / ymax
        n, tu, tn = tot.get(s, (0, 0, 0))
        body.append(text(px, py + 12, name, c["ink"], 13, weight=600))
        if n:
            body.append(text(px, py + 28, f"all {n:,}: {fmt_r(n / tu)} → {fmt_r(n / tn)} seq/s ({fmt_x(tu / tn)})",
                             c["ink2"], 11, extra=NUM))
        v = 0.0
        while v <= ymax + 1e-9:
            body.append(f'<line x1="{l}" y1="{sy(v):.1f}" x2="{r}" y2="{sy(v):.1f}" stroke="{c["axis"] if v == 0 else c["grid"]}" stroke-width="1"/>')
            body.append(text(l - 6, sy(v) + 4, f"{v:g}", c["muted"], 10, "end", extra=NUM))
            v += step
        slot = (r - l) / GROUPS
        for g in range(1, GROUPS + 1):
            cx = l + slot * (g - 0.5)
            k = grp[(s, g)]
            body.append(text(cx, b + 15, f"L{g}", c["ink2"], 11, "middle", 600))
            body.append(text(cx, b + 28, f"{k['lo']:,}-{k['hi']:,}", c["muted"], 10, "middle", extra=NUM))
            if (s, g) not in per: continue
            ru, rn = per[(s, g)]
            body.append(text(cx, b + 42, fmt_x(rn / ru), c["ink2"], 10, "middle", extra=NUM))
            for j, (val, color, who) in enumerate(((ru, up_color, "upstream"), (rn, now_color, "this branch"))):
                x = cx - bw - bgap / 2 + j * (bw + bgap)
                body.append(vbar(x, b, bw, sy(val), color, f"{name} L{g}, {title}, {who}: {fmt_r(val)} sequences/s"))
                body.append(text(x + bw / 2, sy(val) - 4, fmt_r(val), c["ink2"], 10, "middle", extra=NUM))
    body.append(text(10, H - 10, "Each group: 250 proteins (length range in residues) against all 21,979 Pfam-A models with --cut_ga; "
                     "under each group, the speedup.", c["muted"], 11))
    return svg_doc(W, H, c, body)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("json")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    res = json.load(open(a.json))
    labels = [b["label"] for b in res["binaries"]]
    first, last = labels[0], labels[-1]
    med, rss, grp = medians(res), peak_rss(res), groups()
    rt = {sfx: rates(med, grp, sfx, first, last) for sfx, _, _ in FORMS}
    os.makedirs(a.outdir, exist_ok=True)
    for sfx, tag, title in FORMS:
        if not rt[sfx][0]: continue
        for mode, c in MODES.items():
            open(os.path.join(a.outdir, f"annot-{tag}-{mode}.svg"), "w").write(chart(med, grp, sfx, title, first, last, c))

    def cells(sfx, key, tot=False):
        per, totals = rt[sfx]
        if tot:
            if key not in totals: return ["—"] * 3
            n, tu, tn = totals[key]
            return [f"**{fmt_r(n / tu)}**", f"**{fmt_r(n / tn)}**", f"**{fmt_x(tu / tn)}**"]
        if key not in per: return ["—"] * 3
        ru, rn = per[key]
        return [fmt_r(ru), fmt_r(rn), fmt_x(rn / ru)]

    print("Sequences per second (wall time, 16 threads):\n")
    print(f"| Set | Group | Lengths | pressed: {first} | {last} | speedup | text: {first} | {last} | speedup |")
    print("|---|---|---|---:|---:|---:|---:|---:|---:|")
    for s, name in SETS:
        for g in range(1, GROUPS + 1):
            k = grp[(s, g)]
            print(f"| {name if g == 1 else ''} | L{g} | {k['lo']:,}-{k['hi']:,} | " +
                  " | ".join(cells("", (s, g)) + cells("_txt", (s, g))) + " |")
        print("| | **all 1,000** | | " + " | ".join(cells("", s, True) + cells("_txt", s, True)) + " |")
    print("| **All six** | **6,000** | | " + " | ".join(cells("", "all", True) + cells("_txt", "all", True)) + " |")

    print("\nPeak RSS (MB), pressed Pfam-A, max over reps:\n")
    print("| Set | " + " | ".join(f"L{g} {first} | L{g} {last}" for g in range(1, GROUPS + 1)) + " |")
    print("|---|" + "---:|" * (2 * GROUPS))
    for s, name in SETS:
        row = []
        for g in range(1, GROUPS + 1):
            r = rss.get(f"A_{s}_L{g}", {})
            row += [f"{r[l]:.0f}" if l in r else "—" for l in (first, last)]
        print(f"| {name} | " + " | ".join(row) + " |")


if __name__ == "__main__":
    sys.exit(main())
