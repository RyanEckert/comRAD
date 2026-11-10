#!/usr/bin/env python3
################################################################################
#Script Name    : paralogFilter.py
#Description    : FDR-filter ngsParalog calcLR output into a site keep-list
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
"""Turn ngsParalog calcLR likelihood ratios into a keep-list for ANGSD -sites.

calcLR writes one line per site: scaffold, position, -log null likelihood,
-log alternative likelihood, and the likelihood ratio (LR). It does no test
itself. Under the null of correct mapping the LR is a 50:50 mixture of chi2(1)
and chi2(0) (ngsParalog README), so

    p = 0.5 * P(chi2_1 > LR)

p-values are Benjamini-Hochberg adjusted across ALL sites at once, which is why
chunked runs are combined before the correction rather than corrected one by
one. Sites below --fdr are paralog candidates; everything else in the site
list is kept, in its original order. Sites calcLR never wrote (too few
individuals covered) are kept: later ANGSD -minInd filters remove them anyway.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
from scipy.stats import chi2

prog = os.path.basename(sys.argv[0])


def die(msg):
    sys.exit(f"{prog}: {msg}")


def collect(d, out):
    """Check every chunk in d finished, then concatenate the .lr files in order.

    A chunk is a site list named chunk_NNNN (split -d output) with chunk_NNNN.lr
    and chunk_NNNN.log beside it. Any chunk without 'Finished!' in its log, or
    without an .lr, stops the run: an FDR over part of the sites is not the
    same test, so a partial result is worse than none.
    """
    chunks = sorted(f for f in glob.glob(os.path.join(d, "chunk_*"))
                    if re.fullmatch(r"chunk_\d+", os.path.basename(f)))
    if not chunks:
        die(f"no chunk_NNNN site lists in {d}")
    bad = []
    for c in chunks:
        log = c + ".log"
        ok = False
        if os.path.isfile(log):
            with open(log, errors="replace") as fh:
                ok = any("Finished!" in line for line in fh)
        if not ok:
            bad.append(f"{c}\tno 'Finished!' in {log}")
        elif not os.path.isfile(c + ".lr"):
            bad.append(f"{c}\tno {c}.lr")
    if bad:
        die(f"{len(bad)} of {len(chunks)} chunks did not finish; rerun them "
            f"and collect again:\n  " + "\n  ".join(bad))
    with open(out, "w") as fo:
        for c in chunks:
            with open(c + ".lr") as fi:
                for line in fi:
                    fo.write(line)
    print(f"{prog}: {len(chunks)} chunks finished, combined into {out}",
          file=sys.stderr)
    return [out]


def read_lr(paths):
    keys, lr = [], []
    for p in paths:
        with open(p) as fh:
            for i, line in enumerate(fh, 1):
                f = line.split()
                if not f:
                    continue
                if len(f) < 5:
                    die(f"{p} line {i} has {len(f)} columns; calcLR writes 5")
                keys.append((f[0], f[1]))
                lr.append(float(f[4]))
    if not keys:
        die("no sites in the .lr input")
    return keys, np.array(lr)


def bh(p):
    """Benjamini-Hochberg adjusted p-values (same as R p.adjust 'BH')."""
    n = len(p)
    order = np.argsort(p)[::-1]
    ranked = p[order] * n / np.arange(n, 0, -1)
    adj = np.minimum.accumulate(ranked)
    out = np.empty(n)
    out[order] = np.minimum(adj, 1)
    return out


def main():
    ap = argparse.ArgumentParser(
        prog=prog,
        description="FDR-filter ngsParalog calcLR output and write the sites "
                    "to keep for ANGSD -sites.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""examples:
  # one calcLR run
  {prog} -i mcav.lr -s mcavddSites -o mcav

  # chunked run: check every chunk finished, combine into mcav.lr, then filter
  {prog} --collect paralog -s mcavddSites -o mcav

  # stricter
  {prog} -i mcav.lr -s mcavddSites -o mcav --fdr 0.01

  # then, for ANGSD
  angsd sites index mcavKeepSites

outputs:
  <o>.lr             combined calcLR output (--collect only)
  <o>ParalogSites    flagged sites, scaffold <TAB> position
  <o>KeepSites       -s sites minus the flagged ones, original order
  <o>KeepRegions     scaffolds in <o>KeepSites, for ANGSD -rf

p-values:
  p = 0.5 * P(chi2_1 > LR), Benjamini-Hochberg adjusted across every site
  together. calcLR fits each site on its own reads, so chunking does not change
  the LR; the correction is the only step that needs all sites at once.
""")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("-i", "--input", nargs="+", metavar="LR",
                     help="calcLR output file(s)")
    src.add_argument("--collect", metavar="DIR",
                     help="directory of chunk_NNNN site lists with .lr and "
                          ".log beside them; checks all finished and combines "
                          "them into <o>.lr first")
    ap.add_argument("-s", "--sites", required=True, metavar="FILE",
                    help="site list calcLR was run on, scaffold <TAB> position")
    ap.add_argument("-o", "--out", required=True, metavar="PREFIX",
                    help="output prefix")
    ap.add_argument("--fdr", type=float, default=.05, metavar="Q",
                    help="Benjamini-Hochberg threshold (default .05)")
    a = ap.parse_args()

    if not 0 < a.fdr < 1:
        die("--fdr must be between 0 and 1")
    if not os.path.isfile(a.sites):
        die(f"cannot open {a.sites}")

    paths = collect(a.collect, a.out + ".lr") if a.collect else a.input
    keys, lr = read_lr(paths)

    if len(set(keys)) != len(keys):
        die("some sites appear more than once in the .lr input; were chunks "
            "combined twice?")

    p = .5 * chi2.sf(np.clip(lr, 0, None), df=1)
    q = bh(p)
    flagged = {k for k, v in zip(keys, q) if v < a.fdr}

    sites = []
    with open(a.sites) as fi:
        for line in fi:
            f = line.split()
            if len(f) >= 2:
                sites.append((f[0], f[1]))
    stray = len(set(keys) - set(sites))
    if stray:
        die(f"{stray} tested sites are not in {a.sites}; the .lr and the site "
            f"list do not match")

    keep = [k for k in sites if k not in flagged]
    with open(a.out + "KeepSites", "w") as fo:
        fo.write("".join(f"{s}\t{pos}\n" for s, pos in keep))
    with open(a.out + "KeepRegions", "w") as fo:
        fo.write("".join(r + "\n" for r in dict.fromkeys(s for s, _ in keep)))
    with open(a.out + "ParalogSites", "w") as fo:
        fo.write("".join(f"{k[0]}\t{k[1]}\n" for k, v in zip(keys, q)
                         if v < a.fdr))

    print(f"{prog}: {len(sites)} sites in {a.sites}, {len(keys)} tested, "
          f"{len(flagged)} flagged at FDR {a.fdr}, {len(keep)} kept",
          file=sys.stderr)


if __name__ == "__main__":
    main()
