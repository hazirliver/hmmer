#!/usr/bin/env python3
"""hmmsearch's command line on pyhmmer, so dev/perf/series.py and ab.sh can run pyhmmer like a binary.

  pyhmmsearch.py [options] <hmmfile> <seqdb>

Runs the search `hmmsearch` would, with pyhmmer.hmmer.hmmsearch() used the way pyhmmer recommends:
the targets are read into memory once (a DigitalSequenceBlock) and shared by the worker threads, and
the HMMs stream from the file. With fewer HMMs than threads it passes them as a list, so pyhmmer splits
the targets among the threads instead (its `parallel="targets"` mode). The tables are written by
TopHits.write(), which calls HMMER's own p7_tophits_TabularTargets() / TabularDomains(), with the column
header on the first query only, as hmmsearch does. Only the '# Program: ...' block that hmmsearch
appends to a table is missing (series.py --mask-tail drops it before comparing).

Options: the subset dev/perf/vs_pyhmmer.tsv uses. Anything else is an error, not ignored.
  --cpu N          N worker threads; 0 searches in the main thread, as hmmsearch --cpu 0 [2]
  --cut_ga         Pfam GA thresholds (bit_cutoffs="gathering")
  -Z X, -E X       as hmmsearch
  --max            bias filter off, F1 = F2 = F3 = 1 (all hmmsearch --max does for proteins)
  --noali          accepted; only hmmsearch's main output shows alignments
  -o /dev/null     required: pyhmmer can't write hmmsearch's main output
  --tblout F, --domtblout F

Writes where the time went (import, reading the targets, search + writing) to .pyhmmsearch.json in
the current directory; dev/perf/ab.py adds it to the run's sample.

PYHMMSEARCH_CPUS=N in the environment replaces any --cpu above 0 by N, so that ab.sh can compare
thread counts on the same workloads: ab.sh pyhmmer-0.12.3 pyhmmer-0.12.3 --b-env PYHMMSEARCH_CPUS=8

PYHMMSEARCH_FTZ=1 sets flush-to-zero and denormals-are-zero first, as HMMER's impl_Init() does in
hmmsearch. pyhmmer leaves them off on purpose (they would change all of Python's float arithmetic), so
it computes with subnormal numbers, and prints e.g. an E-value of 3.2e-315 where hmmsearch prints 0.
The worker threads inherit the setting. dev/perf/pyhmmer.sh --ftz makes a snapshot with it on.
"""
import time
T0 = time.perf_counter()

import argparse, itertools, json, os, sys     # noqa: E402

import pyhmmer                                  # noqa: E402
from pyhmmer.easel import SequenceFile          # noqa: E402
from pyhmmer.plan7 import HMMFile               # noqa: E402

TIMES = ".pyhmmsearch.json"


def flush_subnormals():
    """FTZ + DAZ in this thread's MXCSR (x86-64), as impl_Init(); threads started later inherit them."""
    import ctypes
    libm = ctypes.CDLL("libm.so.6")
    env = ctypes.create_string_buffer(32)                    # glibc's x86-64 fenv_t: __mxcsr is the last 4 bytes
    if libm.fegetenv(env) != 0: sys.exit("fegetenv() failed")
    mxcsr = int.from_bytes(env.raw[28:32], "little") | 0x8040   # FTZ (bit 15), DAZ (bit 6)
    ctypes.memmove(ctypes.addressof(env) + 28, mxcsr.to_bytes(4, "little"), 4)
    if libm.fesetenv(env) != 0 or float("5e-324") * 1.0 != 0.0: sys.exit("can't set flush-to-zero")


def main():
    ap = argparse.ArgumentParser(prog="hmmsearch (pyhmmer)", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("hmmfile"); ap.add_argument("seqdb")
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--cut_ga", action="store_true")
    ap.add_argument("-Z", type=float); ap.add_argument("-E", type=float)
    ap.add_argument("--max", action="store_true")
    ap.add_argument("--noali", action="store_true")
    ap.add_argument("-o", required=True)
    ap.add_argument("--tblout"); ap.add_argument("--domtblout")
    a = ap.parse_args()
    if a.o != "/dev/null":
        ap.error("-o: only /dev/null (pyhmmer can't write hmmsearch's main output)")
    if a.cpu > 0 and os.environ.get("PYHMMSEARCH_CPUS"): a.cpu = int(os.environ["PYHMMSEARCH_CPUS"])
    cpus = max(1, a.cpu)
    opts = {}
    if a.cut_ga:     opts["bit_cutoffs"] = "gathering"
    if a.Z is not None: opts["Z"] = a.Z
    if a.E is not None: opts["E"] = a.E
    if a.max:        opts.update(bias_filter=False, F1=1.0, F2=1.0, F3=1.0)
    ftz = os.environ.get("PYHMMSEARCH_FTZ") == "1"
    if ftz: flush_subnormals()
    t: dict = dict(pyhmmer=pyhmmer.__version__, cpus=cpus, ftz=ftz, imports=time.perf_counter() - T0)

    with HMMFile(a.hmmfile) as hf:
        first = hf.read()
        if first is None: sys.exit(f"no HMMs in {a.hmmfile}")
        t0 = time.perf_counter()
        with SequenceFile(a.seqdb, digital=True, alphabet=first.alphabet) as sf:   # as hmmsearch: the HMMs' alphabet
            targets = sf.read_block()
        t["read_targets"] = time.perf_counter() - t0
        t["targets"] = len(targets)
        # pyhmmer picks "targets" mode only when it knows there are fewer queries than threads
        head = [first] + list(itertools.islice(hf, cpus - 1))
        few = len(head) < cpus
        queries = head if few else itertools.chain(head, hf)
        t["mode"] = "serial" if cpus == 1 else "targets" if few else "queries"

        t0 = time.perf_counter()
        tbl = open(a.tblout, "wb") if a.tblout else None
        dom = open(a.domtblout, "wb") if a.domtblout else None
        n = 0
        for n, hits in enumerate(pyhmmer.hmmer.hmmsearch(queries, targets, cpus=cpus, **opts), start=1):
            if tbl: hits.write(tbl, format="targets", header=(n == 1))
            if dom: hits.write(dom, format="domains", header=(n == 1))
        for f in (tbl, dom):
            if f: f.close()
        t["search_write"] = time.perf_counter() - t0
        t["queries"] = n
    t["total"] = time.perf_counter() - T0
    with open(TIMES, "w") as f: json.dump(t, f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
