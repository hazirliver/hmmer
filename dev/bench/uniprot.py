#!/usr/bin/env python3
"""Protein sets from UniProt for the Pfam-A annotation benchmark (dev/perf/annot.tsv).

  dev/bench/uniprot.py fetch            download the proteins listed in dev/bench/uniprot/manifest.tsv,
                                        write the FASTA sets, verify them against uniprot/data.sha256
  dev/bench/uniprot.py select [--pin]   (once, done) query UniProt, draw the sets, rewrite the manifest

`select` draws 1,000 proteins per set in dev/bench/uniprot/sets.tsv (UniProtKB/Swiss-Prot,
reviewed, no fragments): it lists every matching accession, ranks them by sha256("<set>:<acc>")
and keeps the first 1,000 not drawn by an earlier set. The 1,000 are sorted by length (ties by
accession) and cut into 4 groups of 250, L1 shortest to L4 longest, so each set has its own
length bounds. The manifest records each protein's set, group, accession, entry name, entry and
sequence versions, length and taxon, plus the UniProt release.

`fetch` reproduces the sets exactly from the manifest, in any later UniProt release: proteins
come from the UniProtKB REST API in batches; one whose sequence version has changed since (or
that was deleted or merged) comes from UniSave, at the manifest's entry version. FASTA headers
are built from the manifest (">sp|ACC|NAME OX=taxid"), not taken from UniProt, so annotation
updates don't change the files. Writes <set>_L<g>.fa (250 proteins, in draw order) and <set>.fa
(all 1,000) to $HMMER_BENCH_DATA/uniprot-<manifest hash>/.
"""
import argparse, gzip, hashlib, os, sys, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import bench   # noqa: E402  (data dir, colors, die, md5)

API = "https://rest.uniprot.org"
SETS = os.path.join(HERE, "uniprot", "sets.tsv")
MANIFEST = os.path.join(HERE, "uniprot", "manifest.tsv")
PINNED = os.path.join(HERE, "uniprot", "data.sha256")
COMMON = " AND reviewed:true AND fragment:false"
PER_SET, GROUPS = 1000, 4
COLS = ["set", "group", "accession", "entry_name", "entry_version", "seq_version", "length", "taxid"]
# UniProt's field names for COLS[2:], and their TSV column headers (the API orders columns its own way)
FIELDS = "accession,id,version,sequence_version,length,organism_id"
HEADERS = {"Entry": "accession", "Entry Name": "entry_name", "Entry version": "entry_version",
           "Sequence version": "seq_version", "Length": "length", "Organism (ID)": "taxid"}


def get(url, tries=5):
    """GET -> (bytes, headers); retries with backoff on network errors and 429/5xx."""
    req = urllib.request.Request(url, headers={"User-Agent": "hmmer-bench (dev/bench/uniprot.py)",
                                               "Accept-Encoding": "gzip"})
    for t in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                body = r.read()
                while body[:2] == b"\x1f\x8b": body = gzip.decompress(body)   # /stream gzips twice
                return body, r.headers
        except Exception as e:
            code = getattr(e, "code", None)
            if t == tries - 1 or (code is not None and code < 500 and code != 429): raise
            time.sleep(2 ** t)
    raise RuntimeError("unreachable")


def read_sets():
    return bench.read_table(SETS, 2)


def read_manifest():
    rows, release = [], None
    for ln in open(MANIFEST):
        if ln.startswith("# UniProt release"): release = ln.split()[3].rstrip(";")
        if ln.startswith("#") or not ln.strip(): continue
        rows.append(dict(zip(COLS, ln.rstrip("\n").split("\t"))))
    return rows, release


def data_dir():
    return bench.uniprot_dir()


# ---------------------------------------------------------------- select
def cmd_select(a):
    taken, rows, release = set(), [], None
    for name, query in read_sets():
        q = query + COMMON
        url = f"{API}/uniprotkb/stream?format=tsv&fields={FIELDS}&query={urllib.parse.quote(q)}"
        body, hdr = get(url)
        rel = hdr.get("X-UniProt-Release")
        if release and rel != release: bench.die(f"UniProt release changed mid-run ({release} -> {rel}); re-run")
        release = release or rel
        lines = body.decode().splitlines()
        cols = [HEADERS[h] for h in lines[0].split("\t")]
        ents = [dict(zip(cols, l.split("\t"))) for l in lines[1:] if l]
        key = lambda e: hashlib.sha256(f"{name}:{e['accession']}".encode()).hexdigest()
        ents.sort(key=key)
        pick = [e for e in ents if e["accession"] not in taken][:PER_SET]
        if len(pick) < PER_SET: bench.die(f"{name}: only {len(pick)} proteins match '{q}'")
        taken.update(e["accession"] for e in pick)
        rank = {e["accession"]: i for i, e in enumerate(pick)}
        by_len = sorted(pick, key=lambda e: (int(e["length"]), e["accession"]))
        size = PER_SET // GROUPS
        for g in range(GROUPS):
            grp = sorted(by_len[g * size:(g + 1) * size], key=lambda e: rank[e["accession"]])
            for e in grp: rows.append(dict(e, set=name, group=str(g + 1)))
            print(f"  {name:10s} L{g + 1}: lengths {by_len[g * size]['length']:>5s}-{by_len[(g + 1) * size - 1]['length']:<5s}"
                  f" (of {len(ents)} matching)")
    with open(MANIFEST, "w") as f:
        f.write(f"# UniProt release {release}; made by dev/bench/uniprot.py select from uniprot/sets.tsv\n")
        f.write("# " + "\t".join(COLS) + "\n")
        for r in rows: f.write("\t".join(r[c] for c in COLS) + "\n")
    print(f"wrote {MANIFEST}: {len(rows)} proteins, UniProt release {release}")
    cmd_fetch(argparse.Namespace(force=True, pin=a.pin))


# ---------------------------------------------------------------- fetch
def parse_fasta(text):
    """{accession: (header, sequence)} from UniProt FASTA."""
    out, acc, hdr, seq = {}, None, None, []
    for ln in text.splitlines() + [">"]:
        if ln.startswith(">"):
            if acc: out[acc] = (hdr, "".join(seq))
            f = ln[1:].split("|")
            acc, hdr, seq = (f[1] if len(f) > 2 else None), ln, []
        else:
            seq.append(ln.strip())
    return out


def sv_of(header):
    for tok in header.split():
        if tok.startswith("SV="): return tok[3:]
    return None


def cmd_fetch(a):
    rows, release = read_manifest()
    d = data_dir()
    done = os.path.join(d, ".complete")
    if os.path.exists(done) and not a.force:
        print(f"UniProt sets present: {d}")
    else:
        os.makedirs(d, exist_ok=True)
        seqs, unisave = {}, []
        accs = [r["accession"] for r in rows]
        t0 = time.time()
        for i in range(0, len(accs), 250):
            batch = accs[i:i + 250]
            body, _ = get(f"{API}/uniprotkb/accessions?format=fasta&accessions={','.join(batch)}")
            got = parse_fasta(body.decode())
            for r in rows[i:i + 250]:
                h = got.get(r["accession"])
                if h and sv_of(h[0]) == r["seq_version"]: seqs[r["accession"]] = h[1]
                else: unisave.append(r)
            print(f"\r  uniprotkb: {min(i + 250, len(accs))}/{len(accs)}", end="", flush=True)
        print()
        for r in unisave:   # sequence changed, or entry gone, since the manifest's release
            body, _ = get(f"{API}/unisave/{r['accession']}?format=fasta&versions={r['entry_version']}")
            got = parse_fasta(body.decode())
            if r["accession"] not in got: bench.die(f"{r['accession']}: not in UniSave at version {r['entry_version']}")
            seqs[r["accession"]] = got[r["accession"]][1]
        if unisave: print(f"  unisave: {len(unisave)} proteins at their manifest entry version")
        for r in rows:
            if len(seqs[r["accession"]]) != int(r["length"]):
                bench.die(f"{r['accession']}: length {len(seqs[r['accession']])}, manifest says {r['length']}")
        files = {}
        for r in rows:
            rec = f">sp|{r['accession']}|{r['entry_name']} OX={r['taxid']}\n"
            s = seqs[r["accession"]]
            rec += "".join(s[j:j + 60] + "\n" for j in range(0, len(s), 60))
            files.setdefault(f"{r['set']}_L{r['group']}.fa", []).append(rec)
            files.setdefault(f"{r['set']}.fa", []).append(rec)
        for fn, recs in files.items():
            with open(os.path.join(d, fn), "w") as f: f.write("".join(recs))
        with open(os.path.join(d, "MANIFEST.sha256"), "w") as m:
            for fn in sorted(files): m.write(f"{bench.md5(os.path.join(d, fn), 'sha256')}  {fn}\n")
        open(done, "w").close()
        print(f"  {len(rows)} proteins (UniProt release {release}) in {time.time() - t0:.0f}s -> {d}")
    if a.pin:
        with open(PINNED, "w") as f: f.write(open(os.path.join(d, "MANIFEST.sha256")).read())
        print(f"pinned checksums into {PINNED}")
    else:
        want = dict(reversed(l.split()) for l in open(PINNED) if l.strip())
        bad = [f for f, h in want.items() if not os.path.exists(os.path.join(d, f)) or bench.md5(os.path.join(d, f), "sha256") != h]
        if bad: bench.die(f"{len(bad)} file(s) differ from {PINNED}: {', '.join(bad[:6])}")
        print(f"  verified: {len(want)} files match dev/bench/uniprot/data.sha256")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select", help="query UniProt and rewrite the manifest (changes the data)")
    s.add_argument("--pin", action="store_true", help="also record the new checksums in uniprot/data.sha256")
    f = sub.add_parser("fetch", help="download the manifest's proteins and write the FASTA sets")
    f.add_argument("--force", action="store_true", help="download again even if present")
    f.add_argument("--pin", action="store_true", help="record the checksums in uniprot/data.sha256")
    a = ap.parse_args()
    return cmd_select(a) if a.cmd == "select" else cmd_fetch(a)


if __name__ == "__main__":
    sys.exit(main())
