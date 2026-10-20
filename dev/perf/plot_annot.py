#!/usr/bin/env python3
"""Charts (SVG, light + dark) and a Markdown table for the Pfam-A annotation benchmark.

  dev/perf/plot_annot.py annot.json --outdir dev/profiling/results

annot.json is a dev/perf/series.py result over dev/perf/annot.tsv or the annotation part of
vs_pyhmmer.tsv (first binary = upstream, last = this branch, any in between, e.g. pyhmmer, shown
between them). Group sizes and length bounds come from dev/bench/uniprot/manifest.tsv.

Speed is throughput: sequences searched per second of wall time (250 per run, against all of
Pfam-A). Writes, per mode (light, dark):
  annot-pressed-<mode>.svg  one grouped bar chart per protein set: length groups L1-L4, one bar per
                            binary, Pfam-A pressed with hmmpress
  annot-text-<mode>.svg     the same from the text Pfam-A file
and prints Markdown tables (sequences per second per group, per set and overall, and peak memory).
--best-of / --drop as in plot_series.py (* in a table: the value of one of the --best-of alternatives).
Colors: upstream gray, this branch blue, the one in between aqua (plot_series.MODES).
"""
import argparse, json, math, os, statistics, sys
from html import escape

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
from plot_series import MODES, text, svg_doc, legend, fmt_x, nice_step, merge_best, best_args   # noqa: E402

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


def rates(med, grp, sfx, labels):
    """{(set, g): {label: seq/s}}, and per set and "all": [n, {label: total wall}] over the groups that
    every label ran"""
    per, tot = {}, {}
    for s, _ in SETS:
        for g in range(1, GROUPS + 1):
            m = med.get(f"A_{s}_L{g}{sfx}")
            if not m or any(l not in m for l in labels): continue
            n = grp[(s, g)]["n"]
            per[(s, g)] = {l: n / m[l] for l in labels}
            for k in (s, "all"):
                t = tot.setdefault(k, [0, dict.fromkeys(labels, 0.0)])
                t[0] += n
                for l in labels: t[1][l] += m[l]
    return per, tot


def names(labels, legend_text=None):
    """legend name per label"""
    legend_text = legend_text or {}
    return {l: "upstream HMMER 3.4" if i == 0 else "this branch" if i == len(labels) - 1 else legend_text.get(l, l)
            for i, l in enumerate(labels)}


def colors(labels, c):
    mid = [c["s3"]]
    return {l: c["muted"] if i == 0 else c["s1"] if i == len(labels) - 1 else mid[(i - 1) % len(mid)]
            for i, l in enumerate(labels)}


def chart(med, grp, sfx, title, labels, c, legend_text=None, won=None, note=None):
    """won: {workload: {label: the --best-of alternative its samples came from}}, named in the tooltips"""
    per, tot = rates(med, grp, sfx, labels)
    col, nm = colors(labels, c), names(labels, legend_text)
    first, last = labels[0], labels[-1]
    nb = len(labels)
    cols = 3 if nb <= 2 else 2                       # three bars per group need wider panels
    gx, gy = 30, 34
    W = 900
    pw = (W - 20 - (cols - 1) * gx) / cols
    ph = 250
    top = 56
    rows = math.ceil(len(SETS) / cols)
    H = top + rows * ph + (rows - 1) * gy + (46 if note else 30)
    vmax = max(v for d in per.values() for v in d.values())
    step = nice_step(vmax / 4)
    ymax = step * math.ceil(vmax * 1.08 / step)
    body = legend(10, 22, c, [(col[l], nm[l]) for l in labels])
    body.append(text(W - 10, 22, f"{title}: sequences per second, 16 threads", c["ink2"], 12, "end"))
    bw, bgap = (19, 3) if nb <= 2 else (22, 3)
    for idx, (s, name) in enumerate(SETS):
        px = 10 + (idx % cols) * (pw + gx)
        py = top + (idx // cols) * (ph + gy)
        l, r, t, b = px + 34, px + pw, py + 42, py + ph - 52
        sy = lambda v: b - (b - t) * v / ymax
        n, tw = tot.get(s, (0, {}))
        body.append(text(px, py + 12, name, c["ink"], 13, weight=600))
        if n and nb <= 2:
            body.append(text(px, py + 28, f"all {n:,}: {fmt_r(n / tw[first])} → {fmt_r(n / tw[last])} seq/s "
                             f"({fmt_x(tw[first] / tw[last])})", c["ink2"], 11, extra=NUM))
        elif n:
            body.append(text(px, py + 28, f"all {n:,}: " + ", ".join(f"{nm[x] if x != first else 'upstream'} "
                             f"{fmt_r(n / tw[x])}" for x in labels) + " seq/s", c["ink2"], 11, extra=NUM))
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
            rg = per[(s, g)]
            if nb <= 2:
                body.append(text(cx, b + 42, fmt_x(rg[last] / rg[first]), c["ink2"], 10, "middle", extra=NUM))
            x = cx - (nb * bw + (nb - 1) * bgap) / 2
            for lab in labels:
                val = rg[lab]
                who = "upstream" if lab == first else nm[lab]
                alt = (won or {}).get(f"A_{s}_L{g}{sfx}", {}).get(lab, lab)
                if alt != lab: who += f" ({(legend_text or {}).get(alt, alt)})"
                body.append(vbar(x, b, bw, sy(val), col[lab], f"{name} L{g}, {title}, {who}: {fmt_r(val)} sequences/s"))
                body.append(text(x + bw / 2, sy(val) - 4, fmt_r(val), c["ink2"], 10, "middle", extra=NUM))
                x += bw + bgap
    body.append(text(10, H - (26 if note else 10), "Each group: 250 proteins (length range in residues) against all "
                     "21,979 Pfam-A models with --cut_ga; " + ("under each group, the speedup." if nb <= 2 else
                     "speedups are in the tables."), c["muted"], 11))
    if note: body.append(text(10, H - 10, note, c["muted"], 11))
    return svg_doc(W, H, c, body)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("json")
    ap.add_argument("--outdir", required=True)
    best_args(ap)
    a = ap.parse_args()
    res = merge_best(json.load(open(a.json)), a.best_of, a.drop)
    labels = [b["label"] for b in res["binaries"]]
    first, last = labels[0], labels[-1]
    won = {w["name"]: w.get("best_of", {}) for w in res["workloads"]}
    med, rss, grp = medians(res), peak_rss(res), groups()
    rt = {sfx: rates(med, grp, sfx, labels) for sfx, _, _ in FORMS}
    os.makedirs(a.outdir, exist_ok=True)
    for sfx, tag, title in FORMS:
        if not rt[sfx][0]: continue
        for mode, c in MODES.items():
            open(os.path.join(a.outdir, f"annot-{tag}-{mode}.svg"), "w").write(
                chart(med, grp, sfx, title, labels, c, dict(kv.split("=", 1) for kv in a.legend), won, a.note))

    # speedups: every label over the first, and the last over each one in between
    pairs = [(x, first) for x in labels[1:]] + [(last, x) for x in labels[1:-1]]

    def cells(sfx, key, tot=False):
        per, totals = rt[sfx]
        if tot:
            if key not in totals: return ["—"] * (len(labels) + len(pairs))
            n, tw = totals[key]
            return ([f"**{fmt_r(n / tw[l])}**" for l in labels] +
                    [f"**{fmt_x(tw[y] / tw[x])}**" for x, y in pairs])
        if key not in per: return ["—"] * (len(labels) + len(pairs))
        rg, wn = per[key], won.get(f"A_{key[0]}_L{key[1]}{sfx}", {})
        return ([fmt_r(rg[l]) + ("*" if wn.get(l, l) != l else "") for l in labels] +
                [fmt_x(rg[x] / rg[y]) for x, y in pairs])

    if len(labels) == 2:
        print("Sequences per second (wall time, 16 threads):\n")
        print(f"| Set | Group | Lengths | pressed: {first} | {last} | speedup | text: {first} | {last} | speedup |")
        print("|---|---|---|---:|---:|---:|---:|---:|---:|")
        forms = [("", "")]
    else:
        forms = [(sfx, title) for sfx, _, title in FORMS if rt[sfx][0]]
    for sfx, title in forms:
        if len(labels) > 2:
            print(f"Sequences per second, {title} (wall time, 16 threads):\n")
            print("| Set | Group | Lengths | " + " | ".join(labels) + " | " +
                  " | ".join(f"{x} vs {y}" for x, y in pairs) + " |")
            print("|---|---|---|" + "---:|" * (len(labels) + len(pairs)))
        both = (lambda k, t=False: cells("", k, t) + cells("_txt", k, t)) if len(labels) == 2 else \
               (lambda k, t=False, sfx=sfx: cells(sfx, k, t))
        for s, name in SETS:
            for g in range(1, GROUPS + 1):
                k = grp[(s, g)]
                print(f"| {name if g == 1 else ''} | L{g} | {k['lo']:,}-{k['hi']:,} | " + " | ".join(both((s, g))) + " |")
            print("| | **all 1,000** | | " + " | ".join(both(s, True)) + " |")
        print("| **All six** | **6,000** | | " + " | ".join(both("all", True)) + " |\n")

    for sfx, _, title in FORMS:
        if not rt[sfx][0]: continue
        print(f"Peak RSS (MB), {title}, max over reps:\n")
        print("| Set | " + " | ".join(f"L{g} {l}" for g in range(1, GROUPS + 1) for l in labels) + " |")
        print("|---|" + "---:|" * (len(labels) * GROUPS))
        for s, name in SETS:
            row = []
            for g in range(1, GROUPS + 1):
                r = rss.get(f"A_{s}_L{g}{sfx}", {})
                row += [f"{r[l]:.0f}" if l in r else "—" for l in labels]
            print(f"| {name} | " + " | ".join(row) + " |")
        print()


if __name__ == "__main__":
    sys.exit(main())
