#!/usr/bin/env python3
"""Charts (SVG, light + dark) and a Markdown table from a dev/perf/series.py result.

  dev/perf/plot_series.py bench-results/series/phases.json --outdir dev/profiling/results

Writes, per mode (light, dark):
  speedup-<mode>.svg   speedup of the last binary over the first one, per scenario
                       (default options, and with --dbcache where the scenario has a variant)
  phases-<mode>.svg    wall time at every binary of the series, one small panel per scenario
and prints a Markdown table of the medians. No dependencies: the SVG is written by hand
(palette and mark specs from the dataviz guidelines: 2 series, validated in both modes).
"""
import argparse, json, math, os, statistics, sys
from html import escape

# Scenario per workload; a workload ending in _dbc is the --dbcache variant of its base.
SCENARIOS = [  # (workload base, short name, detail)
    ("X5_human_pfam1k_cN", "Pfam 1k × human proteome",      "realistic annotation: --cut_ga, 16 threads"),
    ("X4_pfamA_h3m_cN",    "All of Pfam-A × 5 sequences",    "many models, tiny database: pressed HMMs, 16 threads"),
    ("X4_pfamA_txt_cN",    "All of Pfam-A (text) × 5 seqs",  "the same from the text HMM file, 16 threads"),
    ("X3_p1k_s1k_cN",      "Pfam 1k × Swiss-Prot 1k",        "many models, small database, 16 threads"),
    ("X1_p100_50k_cN",     "Pfam 100 × Swiss-Prot 50k",      "16 threads"),
    ("X2_m200_sprot_cN",   "1 model × Swiss-Prot, 16 threads", "one query (M = 200), 570k targets"),
    ("W2_mix_10k",         "40 models × Swiss-Prot 10k",     "mixed model lengths, 1 thread"),
    ("W3_m200_sprot",      "1 model × Swiss-Prot, 1 thread", "M = 200"),
    ("R1_m50_sprot",       "1 short model × Swiss-Prot",     "M = 50, reading-heavy, 1 thread"),
    ("W4_med_longseqs",    "10 models × long sequences",     "targets of 3,000+ residues, 1 thread"),
    ("W1_hitters_10k",     "10 hit-rich models × 10k",       "many domains, -A alignment output, 1 thread"),
    ("W5_med_1k_max",      "10 models × 1k, --max",          "all filters off: domain stage only, 1 thread"),
]

MODES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
                  grid="#e1e0d9", axis="#c3c2b7", s1="#2a78d6", s2="#eb6834"),
    "dark":  dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781",
                  grid="#2c2c2a", axis="#383835", s1="#3987e5", s2="#d95926"),
}
FONT = 'system-ui, -apple-system, "Segoe UI", sans-serif'
S1_NAME = "this branch, default options"
S2_NAME = "this branch with --dbcache"


def medians(res):
    """{workload: {label: median wall}}"""
    out = {}
    for w in res["workloads"]:
        out[w["name"]] = {l: statistics.median(r["wall"] for r in runs) for l, runs in w["runs"].items() if runs}
    return out


def fmt_x(v):  return f"{v:.1f}×" if v < 10 else f"{v:.0f}×"
def fmt_s(v):  return f"{v:.2f} s" if v < 1 else (f"{v:.1f} s" if v < 100 else f"{v:.0f} s")


def text(x, y, s, c, size=12, anchor="start", weight=400, extra=""):
    return (f'<text x="{x:.1f}" y="{y:.1f}" fill="{c}" font-size="{size}" font-weight="{weight}" '
            f'text-anchor="{anchor}" {extra}>{escape(s)}</text>')


def hbar(x0, x1, y, h, color, tip):
    """bar from x0 (baseline) to x1, 4px rounded data end, square at the baseline"""
    r = min(4, (x1 - x0) / 2, h / 2)
    w = x1 - x0
    d = (f"M{x0:.1f},{y:.1f} h{w - r:.1f} a{r},{r} 0 0 1 {r},{r} v{h - 2 * r:.1f} "
         f"a{r},{r} 0 0 1 {-r},{r} h{-(w - r):.1f} z")
    return f'<path d="{d}" fill="{color}"><title>{escape(tip)}</title></path>'


def svg_doc(width, height, c, body):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
            f'font-family=\'{FONT}\'>\n<rect width="100%" height="100%" fill="{c["surface"]}"/>\n' + "\n".join(body) + "\n</svg>\n")


def legend(x, y, c, items):
    out, cx = [], x
    for color, name in items:
        out.append(f'<rect x="{cx:.1f}" y="{y - 9:.1f}" width="12" height="12" rx="3" fill="{color}"/>')
        out.append(text(cx + 18, y + 1, name, c["ink2"], 12))
        cx += 18 + 7.0 * len(name) + 22
    return out


# ------------------------------------------------------------------ speedup chart
def speedup_chart(med, first, last, c):
    rows = []
    for base, name, detail in SCENARIOS:
        m = med.get(base)
        if not m or first not in m or last not in m: continue
        dflt = m[first] / m[last]
        dbc  = med.get(base + "_dbc", {}).get(last)
        rows.append((name, detail, dflt, (m[first] / dbc) if dbc else None, m[first], m[last], dbc))
    W, left, right = 820, 330, 70
    x0, x1 = left, W - right
    vmax = max(max(r[2], r[3] or 0) for r in rows)
    top = max(2.0, vmax * 1.3)                       # room for the value labels past the longest bar
    ticks = [t for t in (1, 2, 5, 10, 20, 50, 100) if t <= top]
    sx = lambda v: x0 + (x1 - x0) * math.log(v) / math.log(top)
    bar, gap, pad = 14, 2, 16
    hrow = 2 * bar + gap                             # every row as tall as a two-bar row
    y = 58
    body = legend(10, 24, c, [(c["s1"], S1_NAME), (c["s2"], S2_NAME)])
    grid, marks = [], []
    for name, detail, dflt, dbc, tup, tlast, tdbc in rows:
        body.append(text(left - 14, y + hrow / 2 - 2, name, c["ink"], 13, "end", 600))
        body.append(text(left - 14, y + hrow / 2 + 13, detail, c["muted"], 11, "end"))
        ya = y if dbc else y + (hrow - bar) / 2      # a lone bar sits in the middle of its row
        marks.append(hbar(x0, sx(dflt), ya, bar, c["s1"], f"{name}: upstream {fmt_s(tup)} → {fmt_s(tlast)} ({fmt_x(dflt)})"))
        marks.append(text(sx(dflt) + 6, ya + bar - 3, fmt_x(dflt), c["ink2"], 12, extra='font-variant-numeric="tabular-nums"'))
        if dbc:
            yb = y + bar + gap
            marks.append(hbar(x0, sx(dbc), yb, bar, c["s2"], f"{name} with --dbcache: upstream {fmt_s(tup)} → {fmt_s(tdbc)} ({fmt_x(dbc)})"))
            marks.append(text(sx(dbc) + 6, yb + bar - 3, fmt_x(dbc), c["ink2"], 12, extra='font-variant-numeric="tabular-nums"'))
        y += hrow + pad
    ybot = y - pad + 8
    for t in ticks:
        grid.append(f'<line x1="{sx(t):.1f}" y1="44" x2="{sx(t):.1f}" y2="{ybot:.1f}" stroke="{c["grid"] if t != 1 else c["axis"]}" stroke-width="1"/>')
        grid.append(text(sx(t), ybot + 16, f"{t}×", c["muted"], 11, "middle"))
    grid.append(text((x0 + x1) / 2, ybot + 36, "speedup over upstream HMMER 3.4 (wall time, log scale)", c["ink2"], 12, "middle"))
    return svg_doc(W, ybot + 48, c, body + grid + marks)


# ------------------------------------------------------------------ progression chart
def phases_chart(med, labels, c, xlabels):
    panels = [(b, n, d) for b, n, d in SCENARIOS if b in med]
    cols, pw, ph = 3, 270, 190
    gx, gy, top = 18, 26, 40
    rows = math.ceil(len(panels) / cols)
    W = cols * pw + (cols - 1) * gx + 20
    H = top + rows * ph + (rows - 1) * gy + 10
    body = legend(10, 22, c, [(c["s1"], "default options"), (c["s2"], "with --dbcache")])
    for idx, (base, name, detail) in enumerate(panels):
        px = 10 + (idx % cols) * (pw + gx)
        py = top + (idx // cols) * (ph + gy)
        series = [(c["s1"], med[base], "default")]
        if base + "_dbc" in med: series.append((c["s2"], med[base + "_dbc"], "--dbcache"))
        vmax = max(v for _, m, _ in series for v in m.values())
        step = nice_step(vmax / 3)
        ymax = step * math.ceil(vmax / step)
        l, r, t, b = px + 40, px + pw - 44, py + 52, py + ph - 22
        sx = lambda i: l + (r - l) * i / (len(labels) - 1)
        sy = lambda v: b - (b - t) * v / ymax
        up, now = med[base][labels[0]], med[base][labels[-1]]
        body.append(text(px, py + 12, name, c["ink"], 13, weight=600))
        body.append(text(px, py + 27, f"{fmt_s(up)} → {fmt_s(now)} ({fmt_x(up / now)})", c["ink2"], 11))
        if len(series) > 1:
            dnow = series[1][1][labels[-1]]
            body.append(text(px, py + 41, f"with --dbcache {fmt_s(dnow)} ({fmt_x(up / dnow)})", c["ink2"], 11))
        v = 0.0
        while v <= ymax + 1e-9:
            body.append(f'<line x1="{l}" y1="{sy(v):.1f}" x2="{r}" y2="{sy(v):.1f}" stroke="{c["axis"] if v == 0 else c["grid"]}" stroke-width="1"/>')
            body.append(text(l - 6, sy(v) + 4, f"{v:g}", c["muted"], 10, "end", extra='font-variant-numeric="tabular-nums"'))
            v += step
        for i, xl in enumerate(xlabels):
            body.append(text(sx(i), b + 14, xl, c["muted"], 10, "middle"))
        for color, m, sname in series:
            pts = [(i, m[lab]) for i, lab in enumerate(labels) if lab in m]
            path = " ".join(f"{'M' if k == 0 else 'L'}{sx(i):.1f},{sy(val):.1f}" for k, (i, val) in enumerate(pts))
            body.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
            for i, val in pts:
                body.append(f'<circle cx="{sx(i):.1f}" cy="{sy(val):.1f}" r="4" fill="{color}" stroke="{c["surface"]}" stroke-width="2">'
                            f'<title>{escape(name)}, {sname}, {xlabels[i]}: {fmt_s(val)}</title></circle>')
        ends = [sy(m[labels[-1]]) for _, m, _ in series]
        if len(ends) == 1 or abs(ends[0] - ends[1]) >= 12:  # end labels only where they don't collide (header has both)
            for (color, m, sname), ye in zip(series, ends):
                body.append(text(sx(len(labels) - 1) + 8, ye + 4, fmt_s(m[labels[-1]]), c["ink2"], 10))
    return svg_doc(W, H, c, body)


def nice_step(x):
    e = 10 ** math.floor(math.log10(x))
    for m in (1, 2, 2.5, 5, 10):
        if m * e >= x: return m * e
    return 10 * e


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("json")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--xlabels", help="comma-separated x labels for the phase chart (default: the series labels)")
    a = ap.parse_args()
    res = json.load(open(a.json))
    labels = [b["label"] for b in res["binaries"]]
    xlabels = a.xlabels.split(",") if a.xlabels else labels
    med = medians(res)
    os.makedirs(a.outdir, exist_ok=True)
    for mode, c in MODES.items():
        open(os.path.join(a.outdir, f"speedup-{mode}.svg"), "w").write(speedup_chart(med, labels[0], labels[-1], c))
        open(os.path.join(a.outdir, f"phases-{mode}.svg"), "w").write(phases_chart(med, labels, c, xlabels))

    # Markdown table: wall-time medians per scenario and binary
    print("| Scenario | " + " | ".join(xlabels) + " | speedup |")
    print("|---|" + "---:|" * (len(labels) + 1))
    for base, name, detail in SCENARIOS:
        for wl, suffix in ((base, ""), (base + "_dbc", " (--dbcache)")):
            m = med.get(wl)
            if not m: continue
            cells = [fmt_s(m[l]) if l in m else "—" for l in labels]
            ref = med[base][labels[0]]
            print(f"| {name}{suffix} | " + " | ".join(cells) + f" | {fmt_x(ref / m[labels[-1]])} |")


if __name__ == "__main__":
    sys.exit(main())
