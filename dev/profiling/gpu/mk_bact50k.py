#!/usr/bin/env python3
"""50,000 random bacterial Swiss-Prot proteins (non-fragment) from raw/sprot_dat.gz, seed 1."""
import gzip, os, random, sys
src = os.path.expanduser("~/.cache/hmmer-bench/raw/sprot_dat.gz")
out = os.path.expanduser("~/.cache/hmmer-bench/g0-data/bacteria_50k.fa")
recs = []
ac = None; oc = None; frag = False; seq = []; insq = False
with gzip.open(src, "rt", errors="replace") as f:
    for ln in f:
        if ln.startswith("ID "): ac = None; oc = None; frag = False; seq = []; insq = False; name = ln.split()[1]
        elif ln.startswith("AC ") and ac is None: ac = ln.split()[1].rstrip(";")
        elif ln.startswith("OC ") and oc is None: oc = ln[5:].split(";")[0].strip()
        elif ln.startswith("DE ") and "Flags: Fragment" in ln: frag = True
        elif ln.startswith("SQ "): insq = True
        elif ln.startswith("//"):
            if oc == "Bacteria" and not frag and seq: recs.append((name, ac, "".join(seq)))
            insq = False
        elif insq: seq.append(ln.replace(" ", "").strip())
print(len(recs), "bacterial non-fragment entries", file=sys.stderr)
random.Random(1).shuffle(recs)
with open(out, "w") as o:
    for name, ac, s in recs[:50000]:
        o.write(f">{name} {ac}\n")
        for i in range(0, len(s), 60): o.write(s[i:i+60] + "\n")
