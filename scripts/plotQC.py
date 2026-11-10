#!/usr/bin/env python3
################################################################################
#Script Name    : plotQC.py
#Description    : Plot ANGSD quality/depth QC and write per-sample coverage lists
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
"""plotQC.py -- Python port of plotQC.R (Fumagalli / Matz).

Same four panels as the R script, laid out 2 x 2, with the counts file streamed
in chunks instead of read as one giant table -- which is what makes the R
version slow, especially with larger sample numbers.

  plotQC.py -p dd -b bams                    # quality5X.txt, 10x marked
  plotQC.py -p dd -b bams -q 3 5 10          # three quality<D>X.txt files
  plotQC.py -p dd -b bams -c 5 10 20         # mark three depths on the panel
  plotQC.py -p dd -b bams --no-counts        # skip the slow bottom row

-q/--quality  writes quality<D>X.txt for each value: the proportion of each
              sample's COVERED sites at depth >= D, worst samples first. This
              is the file filterBams.py -q reads.
-c/--coverage marks each value on the per-sample depth panel with its own
              colour and dash, and reports the median fraction of covered sites
              deeper than it. Nothing is written for these.

Both thresholds are INCLUSIVE: 5 means depth >= 5. The R script used
1-cumul[COV,], i.e. strictly GREATER than COV, so its "coverage of 5" file is
>=6x and corresponds to -q 6 here.

The R script also wrote a bams.qc list, keeping samples whose log coverage
proportion sat within 3 SD of the batch mean. Dropped 2026-09-22: that cut is
relative to whatever else is in the run, so it passes genuinely unusable
samples whenever the spread is wide, and nothing in the pipeline read it. Set
an absolute floor with filterBams.py -p against quality<D>X.txt instead.
"""
import argparse
import gzip
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(
    description="Plot the ANGSD depth/quality QC (a Python port of ngsTools' plotQC.R) and write each sample's "
                "proportion of covered sites at depth >= D, for filtering samples by coverage.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""examples:
  plotQC.py -p sintDepth -b sintAllBams                # sintDepth.pdf, quality5X.txt, 10x marked
  plotQC.py -p sintDepth -b sintAllBams -q 3 5 10      # quality3X.txt, quality5X.txt, quality10X.txt
  plotQC.py -p sintDepth -b sintAllBams -c 5 10 20     # mark three depths on the per-sample panel
  plotQC.py -p sintDepth -b sintAllBams --no-counts    # skip the slow .counts panels

inputs:
  ANGSD run with -doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2 -out <prefix>, which writes
  <prefix>.qs, <prefix>.depthSample, <prefix>.depthGlobal and <prefix>.counts.gz.

thresholds:
  both -q and -c are inclusive: 5 means depth >= 5. quality<D>X.txt lists each sample's
  proportion of its COVERED sites at >= D, worst first; filterBams.py -q reads it.
""")
ap.add_argument("-p", "--prefix", required=True, help="prefix of the angsd -doDepth output")
ap.add_argument("-b", "--bams", required=True, help="bam list, in the order angsd saw them")
ap.add_argument("-q", "--quality", type=int, nargs="+", default=[5], metavar="D",
                help="minimum depths, inclusive; one quality<D>X.txt per value "
                     "(default 5)")
ap.add_argument("-x", "--xl", type=int, default=20, help="max depth to plot (default 20)")
ap.add_argument("-c", "--coverage", type=int, nargs="+", default=[10], metavar="D",
                help="depth cutoffs to mark on the depth panel; one or more (default 10)")
ap.add_argument("--chunk", type=int, default=100000, help="rows per chunk when streaming counts")
ap.add_argument("--no-counts", action="store_true", help="skip the .counts pass (panels 3 and 4)")
a = ap.parse_args()

pref = a.prefix
bams = [l.strip() for l in open(a.bams) if l.strip()]
nind = len(bams)

fig, ax = plt.subplots(2, 2, figsize=(7.5, 7))
ax = ax.ravel()   # 0,1 = top row; 2,3 = bottom row (genotyping rate)

# ---------------------------------------------------------------- panel 1 ---
# base quality score distribution
qs = np.genfromtxt(f"{pref}.qs", names=True, dtype=None, encoding=None)
qscore = np.asarray(qs["qscore"], float)
counts = np.asarray(qs["counts"], float)
perc = 1.0 - np.cumsum(counts) / counts.sum()
ax[0].step(qscore, perc, where="post", color="black")
ax[0].set_xlabel("base Q score cutoff")
ax[0].set_ylabel("fraction of sites remaining")

# ---------------------------------------------------------------- panel 2 ---
# per-sample depth: fraction of COVERED sites remaining at each depth cutoff
deps = np.loadtxt(f"{pref}.depthSample", dtype=np.int64, ndmin=2)
depp = deps[:, 1:]                       # drop the depth-0 column
sums = depp.sum(axis=1).astype(float)    # denominator = covered sites
sums[sums == 0] = np.nan
marks = sorted({d for d in a.coverage if d >= 1})
bad = [d for d in a.coverage if d < 1 or d > depp.shape[1]]
if bad:
    print(f"warning: ignoring --coverage outside 1..{depp.shape[1]}: "
          f"{', '.join(map(str, bad))}", file=sys.stderr)
    marks = [d for d in marks if d <= depp.shape[1]]

# widen the axis if a mark falls beyond it, so every mark is visible
xl = min(max(a.xl, max(marks, default=0)), depp.shape[1])
cum = np.cumsum(depp[:, :xl], axis=1) / sums[:, None]
remain = 1.0 - cum                       # remain[i, d-1] = P(depth > d) for sample i

for i in range(remain.shape[0]):
    ax[1].plot(np.arange(1, xl + 1), remain[i], color=(0, 0, 0, 0.2), lw=0.8)

MARK_C = ["#d7301f", "#0570b0", "#238b45", "#6a51a3", "#d95f0e", "#01665e"]
MARK_S = ["--", ":", "-.", (0, (3, 1, 1, 1)), (0, (5, 1)), (0, (1, 1))]
if marks:
    print("median fraction of covered sites deeper than:")
for i, d in enumerate(marks):
    med = np.nanmedian(remain[:, d - 1])
    c, ls = MARK_C[i % len(MARK_C)], MARK_S[i % len(MARK_S)]
    ax[1].axvline(d,   ls=ls, color=c, lw=1.2)
    ax[1].axhline(med, ls=ls, color=c, lw=1.2, label=f"{d}x: {med:.2f}")
    print(f"  {d}x\t{med:.4f}")
if marks:
    ax[1].legend(fontsize=7, frameon=False, loc="upper right")
ax[1].set_ylim(0, 1)
ax[1].set_xlabel("coverage cutoff")
ax[1].set_ylabel("fraction of sites remaining")

# quality<D>X.txt: proportion of covered sites at depth >= D, per sample, one
# file per --quality value. depp column j (0-based) holds depth j+1, so depth
# >= D is columns D-1 on. Computed from depp directly rather than from
# `remain`, which is capped at xl.
qual = sorted({d for d in a.quality})
bad = [d for d in qual if d < 1 or d > depp.shape[1]]
if bad:
    sys.exit(f"--quality outside 1..{depp.shape[1]}: {', '.join(map(str, bad))}")

for d in qual:
    cX = depp[:, d - 1:].sum(axis=1) / sums
    if len(cX) != nind:
        sys.exit(f"{pref}.depthSample has {len(cX)} rows but {a.bams} has {nind} entries")

    order = np.argsort(cX)
    with open(f"quality{d}X.txt", "w") as fh:
        for i in order:
            fh.write(f"{bams[i]} {cX[i]:.15g}\n")

    print(f"\nproportion of covered sites at depth >= {d}x (worst samples first):")
    for i in order[:10]:
        print(f"  {bams[i]}\t{cX[i]:.6f}")
    print(f"...\nquality{d}X.txt: {nind} samples")

# ------------------------------------------------------------ panels 3, 4 ---
# genotyping rate: per site, the fraction of individuals with any data.
# Streamed as a histogram over 0..nind rather than held as a full vector.
info_lines = []
nsites = None
if not a.no_counts:
    cpath = f"{pref}.counts.gz" if os.path.exists(f"{pref}.counts.gz") else f"{pref}.counts"
    opener = gzip.open if cpath.endswith(".gz") else open
    hist = np.zeros(nind + 1, dtype=np.int64)
    with opener(cpath, "rt") as fh:
        header = fh.readline().split()
        ncol = len(header)
        rows, done = [], 0
        for line in fh:
            rows.append(line)
            if len(rows) >= a.chunk:
                arr = np.array([r.split()[:ncol] for r in rows], dtype=np.int32)
                hist += np.bincount((arr > 0).sum(axis=1), minlength=nind + 1)
                done += len(rows)
                rows = []
                print(f"  counts: {done:,} sites", end="\r", file=sys.stderr)
        if rows:
            arr = np.array([r.split()[:ncol] for r in rows], dtype=np.int32)
            hist += np.bincount((arr > 0).sum(axis=1), minlength=nind + 1)
            done += len(rows)
    print(f"  counts: {done:,} sites", file=sys.stderr)
    nsites = hist.sum()

    # fraction of sites at or above each genotyping rate
    rate = np.arange(nind + 1) / ncol
    surv = hist[::-1].cumsum()[::-1] / nsites
    for k in (2, 3):
        ax[k].step(rate, surv, where="post", color="black")
        ax[k].set_xlabel("genotyping rate cutoff")
        ax[k].set_ylabel("fraction of sites remaining")
    ax[3].set_yscale("log")
    fig.suptitle(f"{nsites} sites")
else:
    for k in (2, 3):
        ax[k].axis("off")

# ------------------------------------------------------------------ .info ---
with open(f"{pref}.info", "w") as fh:
    fh.write("qscore\tcounts\tperc\n")
    for q_, c_, p_ in zip(qscore, counts, perc):
        fh.write(f"{q_:g}\t{c_:g}\t{p_:.6f}\n")
    dg = np.loadtxt(f"{pref}.depthGlobal")
    fh.write("\nGlobal_depth\tpercentile\n")
    for d, v in enumerate(np.cumsum(dg) / dg.sum()):
        fh.write(f"{d}\t{v:.6f}\n")

fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(f"{pref}.pdf")
print(f"{pref}.pdf, {pref}.info written")
