#!/usr/bin/env python3
"""Turn whatever PCAngsd wrote for per-individual inbreeding into one TSV.

    ./pcangsdFisToText.py mcavPcangsd -o mcavFis.txt

For each run NAME it looks for the inbreeding output next to NAME + the
usual PCAngsd suffixes, pairs the values with the bam names in NAME.keep, and
writes one TSV:

    bam <TAB> Fis                           one run, no --lineages
    bam <TAB> run <TAB> Fis                 more than one run
    ... <TAB> lineage <TAB> Fis             with --lineages FILE

--lineages takes a whitespace-separated file with the bam in column 1 and its
lineage in column 2 (e.g. mcavBamsClusters). Column 1 may be a path; the
directory and a trailing .bam are dropped before matching to the .keep names.

PCAngsd's output format has changed across versions and the file name has
changed with it.  Rather than guess, this tries every name and every encoding:

    NAME.inbreed.samples.npy     v1.x  (numpy binary)
    NAME.inbreed.samples         v1.x  built without numpy saving
    NAME.inbreed.npy             some 0.99x builds
    NAME.inbreed                 v0.9x (np.savetxt -- PLAIN TEXT)
    NAME.indf.npy / NAME.indf    very old builds
    ... each also accepted with a .gz suffix

`.npy` is read with numpy if numpy is installed; if it is not, the header is
parsed directly, so this works on a login node with a bare Python.  Plain-text
files are read as whitespace-separated numbers, whether they are one per line
or all on one line.

The only hard requirement is that the number of values equals the number of
lines in NAME.keep -- PCAngsd emits one value per individual in the order of
the beagle it was given, which is the order of the .keep file that produced it.
"""
import argparse
import glob
import gzip
import os
import struct
import sys

SUFFIXES = [".inbreed.samples.npy", ".inbreed.samples",
            ".inbreed.npy", ".inbreed",
            ".indf.npy", ".indf",
            "Fis.inbreed.samples.npy", "Fis.inbreed.samples",
            "Fis.inbreed.npy", "Fis.inbreed"]


def open_maybe_gz(path, binary=False):
    if path.endswith(".gz"):
        return gzip.open(path, "rb" if binary else "rt")
    return open(path, "rb" if binary else "r")


def read_npy(path):
    """Read a 1-D .npy without requiring numpy."""
    try:
        import numpy as np
        return [float(x) for x in np.load(path).ravel()]
    except ImportError:
        pass
    with open_maybe_gz(path, binary=True) as fh:
        magic = fh.read(6)
        if magic != b"\x93NUMPY":
            raise ValueError(f"{path} is not a .npy file")
        major = fh.read(1)[0]
        fh.read(1)                                   # minor
        hlen = struct.unpack("<H" if major == 1 else "<I",
                             fh.read(2 if major == 1 else 4))[0]
        header = fh.read(hlen).decode("latin1")
        if "'descr': '<f8'" in header or '"descr": "<f8"' in header:
            fmt, size = "<d", 8
        elif "'descr': '<f4'" in header or '"descr": "<f4"' in header:
            fmt, size = "<f", 4
        else:
            raise ValueError(f"{path}: unsupported dtype in header {header!r}. "
                             "Install numpy and rerun.")
        if "'fortran_order': True" in header:
            raise ValueError(f"{path}: fortran order; install numpy and rerun.")
        raw = fh.read()
        return [struct.unpack_from(fmt, raw, i)[0]
                for i in range(0, len(raw) - len(raw) % size, size)]


def read_text(path):
    with open_maybe_gz(path) as fh:
        vals = []
        for line in fh:
            for tok in line.split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    raise ValueError(
                        f"{path}: found non-numeric token {tok!r}. If this "
                        "file has a header or more than one column per "
                        "individual it is not the per-sample inbreeding "
                        "output -- check which file PCAngsd actually wrote.")
        return vals


def find_output(name):
    cands = []
    for suf in SUFFIXES:
        cands.extend([name + suf, name + suf + ".gz"])
    # last resort: anything that looks like it
    cands.extend(sorted(glob.glob(name + "*inbreed*")) +
                 sorted(glob.glob(name + "*indf*")))
    seen, out = set(), []
    for c in cands:
        if c not in seen and os.path.isfile(c):
            seen.add(c)
            out.append(c)
    return out


def load(path):
    if ".npy" in path:
        return read_npy(path)
    try:
        return read_text(path)
    except ValueError:
        return read_npy(path)          # mislabelled binary


def bam_name(s):
    s = os.path.basename(s)
    return s[:-4] if s.endswith(".bam") else s


def read_lineages(path):
    lin = {}
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            f = line.split()
            if not f:
                continue
            if len(f) < 2:
                sys.exit(f"{path} line {n}: need bam and lineage columns")
            b = bam_name(f[0])
            if b in lin and lin[b] != f[1]:
                sys.exit(f"{path}: {b} listed as both {lin[b]} and {f[1]}")
            lin[b] = f[1]
    return lin


def main():
    prog = os.path.basename(sys.argv[0])
    ap = argparse.ArgumentParser(
        prog=prog,
        description="Collect PCAngsd per-individual inbreeding (--inbreed-"
                    "samples) into one TSV, paired with bam names from "
                    "NAME.keep.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""examples:
  # one run: bam, Fis
  {prog} mcavPcangsd -o mcavFis.txt

  # add each individual's lineage from the NGSadmix assignment
  {prog} mcavPcangsd --lineages ../angsd/mcavBamsClusters -o mcavFis.txt

  # several runs: bam, run, Fis
  {prog} mcavCluster1Pcangsd mcavCluster2Pcangsd -o mcavLineageFis.txt

output columns:
  bam, then run (only with more than one run), then lineage (only with
  --lineages), then Fis.
""")
    ap.add_argument("groups", nargs="+", metavar="NAME",
                    help="PCAngsd -o prefix; NAME.keep (bam names, one per "
                         "line, in beagle order) must exist alongside")
    ap.add_argument("-o", "--out", required=True, metavar="FILE",
                    help="output TSV")
    ap.add_argument("--lineages", metavar="FILE",
                    help="bam in column 1, lineage in column 2 (e.g. "
                         "mcavBamsClusters); adds a lineage column")
    ap.add_argument("--keep-suffix", default=".keep",
                    help="suffix of the bam-name file (default .keep)")
    a = ap.parse_args()

    lin = read_lineages(a.lineages) if a.lineages else None

    rows = []
    for g in a.groups:
        keepf = g + a.keep_suffix
        if not os.path.isfile(keepf):
            sys.exit(f"missing {keepf}")
        ids = [ln.strip() for ln in open(keepf) if ln.strip()]

        cands = find_output(g)
        if not cands:
            sys.exit(f"no PCAngsd inbreeding output found for {g}. Files "
                     f"present: {sorted(glob.glob(g + '*'))}")

        vals, used, errs = None, None, []
        for c in cands:
            try:
                v = load(c)
            except Exception as e:                       # noqa: BLE001
                errs.append(f"  {c}: {e}")
                continue
            if len(v) == len(ids):
                vals, used = v, c
                break
            errs.append(f"  {c}: {len(v)} values, expected {len(ids)}")
        if vals is None:
            sys.exit(f"could not read a usable file for {g}:\n" +
                     "\n".join(errs))

        print(f"{g:10s} {used:40s} n={len(vals):4d}  "
              f"mean Fis={sum(vals)/len(vals):+.4f}")
        rows.extend((i, g, x) for i, x in zip(ids, vals))

    if lin is not None:
        missing = sorted({i for i, _, _ in rows if i not in lin})
        if missing:
            sys.exit(f"{len(missing)} bams not in {a.lineages}:\n  " +
                     "\n  ".join(missing))

    multi = len(a.groups) > 1
    with open(a.out, "w") as fh:
        fh.write("bam" + ("\trun" if multi else "") +
                 ("\tlineage" if lin is not None else "") + "\tFis\n")
        for i, g, x in rows:
            fh.write(i + (f"\t{g}" if multi else "") +
                     (f"\t{lin[i]}" if lin is not None else "") +
                     f"\t{x:.6f}\n")
    print(f"\nwrote {a.out}: {len(rows)} individuals")


if __name__ == "__main__":
    main()
