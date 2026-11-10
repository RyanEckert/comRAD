#!/usr/bin/env python3
"""Subset an ANGSD beagle genotype-likelihood file by individual.

ANGSD `-doGlf 2` output has three leading columns (marker, allele1, allele2)
followed by three columns per individual, in the same order as the bamlist
supplied to ANGSD.

Usage:
    subsetBeagle.py -b xestoNoClones.beagle.gz \
                    -s bams \
                    -k subsample01.keep \
                    -o subsample01.beagle.gz

  -b/--beagle   input beagle(.gz)
  -s/--samples  file listing ALL sample names, one per line, in bamlist order
  -k/--keep     file listing the sample names to retain, one per line
  -o/--out      output beagle(.gz); gzipped if the name ends in .gz

Sample names are matched exactly after stripping whitespace. The kept samples
are written in the order they appear in --keep, and that same order is echoed
to stderr and to <out>.order so downstream row order is unambiguous.
"""
import argparse
import gzip
import io
import os
import sys


def opener(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def read_list(path):
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def main():
    ap = argparse.ArgumentParser(
        description="Subset an ANGSD beagle genotype-likelihood file (-doGlf 2) to some individuals.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""examples:
  subsetBeagle.py -b sintNeutral.beagle.gz -s ../angsd/sintBams -k ../angsd/sintNeutralBams -o sintNeutralNoKin.beagle.gz
      the neutral beagle with kin left out (both lists are bam paths: names must match exactly)
  subsetBeagle.py -b xesto.beagle.gz -s samples.txt -k subsample01.keep -o subsample01.beagle.gz

notes:
  --samples lists every individual in the beagle, in bamlist order; --keep lists the ones to
  keep, in the order to write them. The order written is echoed to stderr and <out>.order.
  The output is gzipped if its name ends in .gz.
""")
    ap.add_argument("-b", "--beagle", required=True, help="input beagle (.gz or plain)")
    ap.add_argument("-s", "--samples", required=True, help="every sample in the beagle, one per line, in bamlist order")
    ap.add_argument("-k", "--keep", required=True, help="the samples to keep, one per line, in output order")
    ap.add_argument("-o", "--out", required=True, help="output beagle; gzipped if it ends in .gz")
    a = ap.parse_args()

    allsamps = read_list(a.samples)
    keep = read_list(a.keep)

    index = {}
    for i, s in enumerate(allsamps):
        if s in index:
            sys.exit(f"ERROR: duplicate sample name in {a.samples}: {s}")
        index[s] = i

    missing = [s for s in keep if s not in index]
    if missing:
        sys.exit(f"ERROR: {len(missing)} sample(s) in --keep not found in "
                 f"--samples, first few: {missing[:5]}")

    # column indices to retain: 0,1,2 then 3+3i .. 5+3i for each kept sample
    cols = [0, 1, 2]
    for s in keep:
        i = index[s]
        cols.extend([3 + 3 * i, 4 + 3 * i, 5 + 3 * i])

    nexp = 3 + 3 * len(allsamps)
    nline = 0
    with opener(a.beagle, "r") as fin, opener(a.out, "w") as fout:
        for line in fin:
            f = line.rstrip("\n").split("\t")
            if len(f) == 1:
                f = line.rstrip("\n").split()
            if nline == 0 and len(f) != nexp:
                sys.exit(f"ERROR: beagle has {len(f)} columns but --samples "
                         f"implies {nexp} ({len(allsamps)} individuals). "
                         "Is the sample list the same one given to ANGSD?")
            fout.write("\t".join(f[c] for c in cols) + "\n")
            nline += 1

    with open(a.out + ".order", "w") as fh:
        fh.write("\n".join(keep) + "\n")

    sys.stderr.write(f"{os.path.basename(a.out)}: {len(keep)} individuals, "
                     f"{nline - 1} sites\n")


if __name__ == "__main__":
    main()

