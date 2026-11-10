#!/usr/bin/env python3
"""Turn a PCAngsd --selection scan into per-SNP p-values and an outlier list.

PCAngsd --selection writes <out>.selection: one row per site that survived its
own MAF filter, one column per PC, each a chi-square statistic with 1 df
(Meisner et al. 2021). --sites-save writes <out>.sites: one 0/1 per beagle
site saying which rows those are. This script lines the statistics back up
with the beagle markers, converts them to p-values, applies Benjamini-Hochberg
within each PC, and calls a SNP an outlier if its q-value is below --fdr on
any PC.
"""
import argparse
import gzip
import math
import os
import sys

prog = os.path.basename(sys.argv[0])


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def markers(beagle):
    """chr and pos of every beagle site, in order. The marker is chr_pos;
    contig names can contain '_', so split on the last one."""
    out = []
    with opener(beagle) as fh:
        fh.readline()
        for line in fh:
            m = line.split("\t", 1)[0]
            c, _, p = m.rpartition("_")
            if not c or not p.isdigit():
                sys.exit(f"{prog}: {beagle}: marker {m!r} is not chr_pos")
            out.append((c, int(p)))
    return out


def bh(p):
    """Benjamini-Hochberg q-values, same order as p."""
    n = len(p)
    order = sorted(range(n), key=lambda i: p[i])
    q = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, p[i] * n / rank)
        q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser(
        prog=prog,
        description="Per-SNP p-values and FDR outliers from a PCAngsd "
                    "--selection scan.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""examples:
  # PCAngsd run with --selection --sites-save -o mcavSelPcangsd
  {prog} -b ../angsd/mcav.beagle.gz -s mcavSelPcangsd -o mcav

  # stricter cutoff
  {prog} -b ../angsd/mcav.beagle.gz -s mcavSelPcangsd --fdr .01 -o mcav

output:
  <out>Selection.tsv   chr, pos, then stat / p / q for each PC, every SNP
                       PCAngsd tested
  <out>Outliers        chr <TAB> pos of SNPs with q < --fdr on any PC

p = P(chi-square, 1 df > stat), as in the PCAngsd README. q is
Benjamini-Hochberg within each PC.
""")
    ap.add_argument("-b", "--beagle", required=True, metavar="FILE",
                    help="the beagle PCAngsd was run on")
    ap.add_argument("-s", "--scan", required=True, metavar="PREFIX",
                    help="PCAngsd -o prefix; needs PREFIX.selection and "
                         "PREFIX.sites (--sites-save)")
    ap.add_argument("--fdr", type=float, default=.05, metavar="Q",
                    help="q-value cutoff (default .05)")
    ap.add_argument("-o", "--out", required=True, metavar="PREFIX",
                    help="output prefix")
    a = ap.parse_args()

    for f in (a.scan + ".selection", a.scan + ".sites"):
        if not os.path.isfile(f):
            sys.exit(f"{prog}: missing {f}")

    sites = markers(a.beagle)
    keep = [int(l) for l in open(a.scan + ".sites") if l.strip()]
    if len(keep) != len(sites):
        sys.exit(f"{prog}: {a.scan}.sites has {len(keep)} lines, "
                 f"{a.beagle} has {len(sites)} sites")
    tested = [s for s, k in zip(sites, keep) if k]

    stats = [[float(x) for x in l.split()] for l in open(a.scan + ".selection")
             if l.strip()]
    if len(stats) != len(tested):
        sys.exit(f"{prog}: {a.scan}.selection has {len(stats)} rows, "
                 f"{a.scan}.sites marks {len(tested)} sites")
    npc = len(stats[0])

    p = [[math.erfc(math.sqrt(max(r[j], 0) / 2)) for r in stats]
         for j in range(npc)]
    q = [bh(pj) for pj in p]
    out = [any(q[j][i] < a.fdr for j in range(npc)) for i in range(len(stats))]

    with open(a.out + "Selection.tsv", "w") as fh:
        fh.write("chr\tpos\t" + "\t".join(f"stat{j+1}\tp{j+1}\tq{j+1}"
                                          for j in range(npc)) + "\n")
        for i, (c, pos) in enumerate(tested):
            fh.write(f"{c}\t{pos}\t" + "\t".join(
                f"{stats[i][j]:.6g}\t{p[j][i]:.6g}\t{q[j][i]:.6g}"
                for j in range(npc)) + "\n")
    with open(a.out + "Outliers", "w") as fh:
        for i, (c, pos) in enumerate(tested):
            if out[i]:
                fh.write(f"{c}\t{pos}\n")

    print(f"{len(sites)} beagle sites, {len(tested)} tested on {npc} PC(s), "
          f"{sum(out)} outliers at q < {a.fdr:g}")
    print(f"wrote {a.out}Selection.tsv, {a.out}Outliers")


if __name__ == "__main__":
    main()
