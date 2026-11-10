#!/usr/bin/env python3
"""Keep one SNP per 2bRAD tag, optionally dropping outlier tags first, and cut
a beagle down to those SNPs.

SNPs on the same 36 bp tag are inherited together, so programs that treat
sites as independent (NGSadmix, evalAdmix, PCAngsd PCA, NGSremix) count them
more than once. Sorted by position, a SNP less than --tag bp after the
previous one on the same contig is on the same tag; the SNP with the highest
minor allele frequency is kept from each tag. With --exclude, every SNP within
--window bp of an excluded SNP (an outlier) is dropped first, which removes
the whole tag it sits on.
"""
import argparse
import bisect
import gzip
import os
import sys

prog = os.path.basename(sys.argv[0])


def opener(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def main():
    ap = argparse.ArgumentParser(
        prog=prog,
        description="One SNP per 2bRAD tag (highest MAF), outlier tags "
                    "removed, and the matching cut-down beagle.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""examples:
  # all individuals: drop outlier tags, one SNP per tag
  {prog} -m mcav.mafs.gz -b mcav.beagle.gz -x mcavOutliers -o mcavUnlinked

  # one lineage, same outliers
  {prog} -m lineages/mcavCluster1.mafs.gz -b lineages/mcavCluster1.beagle.gz \\
      -x mcavOutliers -o lineages/mcavCluster1Unlinked

  # linkage only, no outlier removal
  {prog} -m mcav.mafs.gz -b mcav.beagle.gz -o mcavThin

output:
  <out>.beagle.gz   the beagle, header plus the kept SNPs, in beagle order
  <out>.sites       chr <TAB> pos of the kept SNPs
""")
    ap.add_argument("-m", "--mafs", required=True, metavar="FILE",
                    help="ANGSD .mafs(.gz) from the same run as the beagle")
    ap.add_argument("-b", "--beagle", required=True, metavar="FILE",
                    help="ANGSD .beagle.gz")
    ap.add_argument("-x", "--exclude", metavar="FILE",
                    help="chr <TAB> pos of SNPs whose tags are dropped "
                         "(e.g. pcangsdOutliers.py output)")
    ap.add_argument("--tag", type=int, default=36, metavar="BP",
                    help="tag length; SNPs closer than this are on one tag "
                         "(default 36)")
    ap.add_argument("--window", type=int, default=35, metavar="BP",
                    help="drop SNPs within this many bp of an excluded SNP "
                         "(default 35)")
    ap.add_argument("-o", "--out", required=True, metavar="PREFIX",
                    help="output prefix")
    a = ap.parse_args()

    snps = []
    with opener(a.mafs) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        try:
            ic, ip = head.index("chromo"), head.index("position")
        except ValueError:
            sys.exit(f"{prog}: {a.mafs}: no chromo/position columns")
        maf = next((head.index(c) for c in ("knownEM", "unknownEM", "phat")
                    if c in head), None)
        if maf is None:
            sys.exit(f"{prog}: {a.mafs}: no knownEM/unknownEM/phat column")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            m = float(f[maf])
            snps.append((f[ic], int(f[ip]), min(m, 1 - m)))

    ex = {}
    if a.exclude:
        for line in open(a.exclude):
            f = line.split()
            if len(f) >= 2:
                ex.setdefault(f[0], []).append(int(f[1]))
        for v in ex.values():
            v.sort()

    def near(c, p):
        v = ex.get(c)
        if not v:
            return False
        i = bisect.bisect_left(v, p - a.window)
        return i < len(v) and v[i] <= p + a.window

    kept_ex = [s for s in snps if not near(s[0], s[1])]

    best = {}
    tag = 0
    prev = None
    for c, p, m in sorted(kept_ex, key=lambda s: (s[0], s[1])):
        if prev is None or c != prev[0] or p - prev[1] >= a.tag:
            tag += 1
        prev = (c, p)
        if tag not in best or m > best[tag][2]:
            best[tag] = (c, p, m)
    keep = {(c, p) for c, p, _ in best.values()}

    found = 0
    with opener(a.beagle) as fi, gzip.open(a.out + ".beagle.gz", "wt") as fo:
        fo.write(fi.readline())
        for line in fi:
            c, _, p = line.split("\t", 1)[0].rpartition("_")
            if (c, int(p)) in keep:
                fo.write(line)
                found += 1
    if found != len(keep):
        sys.exit(f"{prog}: {len(keep)} SNPs kept but {found} found in "
                 f"{a.beagle}; are the mafs and beagle from the same run?")

    with open(a.out + ".sites", "w") as fh:
        for c, p in sorted(keep):
            fh.write(f"{c}\t{p}\n")

    print(f"{len(snps)} SNPs, {len(snps) - len(kept_ex)} within {a.window} bp "
          f"of {sum(len(v) for v in ex.values())} excluded SNPs, "
          f"{len(keep)} tags kept")
    print(f"wrote {a.out}.beagle.gz, {a.out}.sites")


if __name__ == "__main__":
    main()
