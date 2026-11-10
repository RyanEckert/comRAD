#!/usr/bin/env python3
"""pickQ.py -- pick the best ngsAdmix run within the dominant mode, per K,
and score model fit against K from the evalAdmix residuals.

For each K, all replicate runs are compared pairwise after optimally matching
their clusters (Hungarian assignment, exact at any K). Runs are grouped into
modes by average-linkage clustering on that similarity, and the highest-
likelihood run in the largest mode is chosen.

Similarity is 1 - ||Q1 - Q2||^2 / (2n): 1 for identical solutions, 0 for
maximally different ones.

  pickQ.py --min 1 --max 15 --glob 2bradNoClones
  pickQ.py --min 1 --max 15 --glob 2bradNoClones --eval 2bradNoClones.beagle.gz
  pickQ.py --glob 2bradNoClones --fit
"""
import argparse
import glob
import os
import re
import shutil
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

prog = os.path.basename(sys.argv[0])

p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter,
                            epilog=f"""examples:
  # pass 1: modes, best run per K, and the evalAdmix job list
  {prog} --min 1 --max 10 --glob gommcav --eval mcavNoClones.beagle.gz

  # same, for NGSadmix runs made with ngsAdmixLauncher.py --minMaf 0.05
  {prog} --min 1 --max 10 --glob gommcav --eval mcavNoClones.beagle.gz --minMaf 0.05

  # submit that job list and wait for it
  launcher_creator.py -j gommcav_evalAdmixJobs -n evalAdmix -q shortq7 -t 6:00:00 -e $EMAIL
  sbatch evalAdmix.slurm

  # pass 2: model fit against K, from the residual matrices those jobs wrote
  {prog} --glob gommcav --fit

why two passes:
  --fit needs the *_best.corres files evalAdmix produces, which do not exist
  when the job list is written, so it cannot run in the same invocation. It
  skips the mode analysis entirely -- that work is already on disk -- and
  reads inMode back out of <glob>_bestQ.txt for the support column.

--minMaf:
  evalAdmix filters the beagle with its own -minMaf (built-in default 0.05)
  and stops if the .fopt has a different number of sites. It must match the
  -minMaf the NGSadmix runs used; ngsAdmixLauncher.py defaults to 0, and so
  does this.

the elbow:
  Kneedle on the RMSD of the off-diagonal residual correlations: the K
  furthest from the straight line joining the first and last points. Same
  quantity and same formula as admix_knee() in admixfit.R, so the two cannot
  disagree on the same input. Residual RMSD is used rather than the
  likelihood because it plateaus once K captures the structure, while the
  likelihood keeps climbing with every added cluster.

  It is a geometric heuristic, not a test. It always returns something, it
  moves if the range of K changes, and on a shallow curve it is close to
  arbitrary -- read it next to the inMode column, and next to
  CLUMPAK/StructureSelector.
""")
p.add_argument("--min",  type=int, default=1, help="minimum K (default 1)")
p.add_argument("--max",  type=int, default=10, help="maximum K (default 10)")
p.add_argument("--thr",  type=float, default=0.9,
               help="similarity above which two runs are the same mode (default 0.9)")
p.add_argument("--glob", default="*", help="prefix glob for Q files (default '*')")
p.add_argument("--out",  help="table to write (default <glob>_bestQ.txt, "
                              "or <glob>_fit.txt with --fit)")
p.add_argument("--eval", metavar="BEAGLE",
               help="also write an evalAdmix job list using this beagle file")
p.add_argument("--minMaf", type=float, default=0,
               help="-minMaf for the evalAdmix jobs; must match the NGSadmix "
                    "runs (default 0, as ngsAdmixLauncher.py)")
p.add_argument("--fit", action="store_true",
               help="second pass: score model fit against K from the "
                    "*_best.corres residual matrices and mark the elbow")
p.add_argument("--threads", type=int, default=4, help="-P for the evalAdmix jobs")
p.add_argument("--no-copy", action="store_true",
               help="don't copy winners to *_best.qopt")
a = p.parse_args()

out = a.out or (f"{a.glob}_bestQ.txt" if a.glob != "*" else "bestQ.txt")


def lik(qfile):
    """Final log likelihood from the matching .log file."""
    for line in open(re.sub(r"\.qopt$", ".log", qfile)):
        m = re.search(r"like=(\S+)", line)
        if m:
            return float(m.group(1))
    return float("nan")


def similarity(x, y):
    """1 - ||x - y||^2 / (2n) with clusters optimally matched."""
    n = x.shape[0]
    cost = ((x[:, :, None] - y[:, None, :]) ** 2).sum(axis=0)
    r, c = linear_sum_assignment(cost)
    return 1.0 - cost[r, c].sum() / (2 * n)


def avg_linkage(D, thr):
    """Agglomerative average-linkage clustering, cut at distance thr."""
    cl = [[i] for i in range(D.shape[0])]
    while len(cl) > 1:
        best = min(((np.mean([[D[i, j] for j in cl[b]] for i in cl[aa]]), aa, b)
                    for aa in range(len(cl)) for b in range(aa + 1, len(cl))),
                   key=lambda t: t[0])
        if best[0] > thr:
            break
        cl[best[1]] += cl[best[2]]
        del cl[best[2]]
    lab = np.empty(D.shape[0], dtype=int)
    for k, c in enumerate(cl):
        for i in c:
            lab[i] = k
    return lab


# ----------------------------------------------------------------- --fit ---
def knee(k, y):
    """K furthest from the chord joining the first and last point (Kneedle).

    Both axes are scaled to [0, 1] and the chord of a decreasing curve is then
    y = 1 - x, so the distance to maximise is |yn - (1 - xn)|. This is the
    same expression admix_knee() uses in admixfit.R, deliberately -- two
    implementations of the same estimator should not be able to disagree.

    Returns None when there is nothing to find: fewer than three K, or a
    curve with no spread to measure against.
    """
    if len(k) < 3:
        return None
    k = np.asarray(k, float)
    y = np.asarray(y, float)
    kr, yr = k.max() - k.min(), y.max() - y.min()
    if kr == 0 or yr == 0:
        return None
    xn = (k - k.min()) / kr
    yn = (y - y.min()) / yr
    return int(k[np.argmax(np.abs(yn - (1 - xn)))])


def read_inmode():
    """inMode per K from the first pass, for the support column."""
    src = f"{a.glob}_bestQ.txt" if a.glob != "*" else "bestQ.txt"
    try:
        with open(src) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            ik, im = head.index("K"), head.index("inMode")
            return {int(c[ik]): c[im]
                    for c in (l.rstrip("\n").split("\t") for l in fh)
                    if len(c) > max(ik, im)}
    except (OSError, ValueError) as err:
        print(f"{prog}: no usable {src} ({err}); inMode left blank",
              file=sys.stderr)
        return {}


def fit():
    """Model fit against K from the evalAdmix residual correlation matrices."""
    files = sorted(glob.glob(f"{a.glob}_k*_best.corres"),
                   key=lambda s: int(re.search(r"_k(\d+)_", s).group(1)))
    if not files:
        sys.exit(f"{prog}: no files matched {a.glob}_k*_best.corres -- have the "
                 f"evalAdmix jobs from --eval finished?")

    inmode = read_inmode()
    rows, bad = [], []
    for f in files:
        K = int(re.search(r"_k(\d+)_", f).group(1))
        if not (a.min <= K <= a.max):
            continue
        # one unreadable matrix reports itself and is skipped; the rest of the
        # curve is still worth having
        try:
            r = np.genfromtxt(f)                  # NA / NaN become nan
            if r.ndim != 2 or r.shape[0] != r.shape[1]:
                raise ValueError(f"not a square matrix, shape {r.shape}")
            v = r[np.tril_indices_from(r, -1)]    # strict lower triangle
            if not np.any(np.isfinite(v)):
                raise ValueError("no finite off-diagonal values")
        except (OSError, ValueError) as err:
            print(f"{prog}: {f}: {err}", file=sys.stderr)
            bad.append(f)
            continue
        rows.append([K,
                     float(np.sqrt(np.nanmean(v ** 2))),
                     float(np.nanmean(np.abs(v))),
                     float(np.nanmax(np.abs(v))),
                     int(np.count_nonzero(~np.isfinite(v)))])

    if not rows:
        sys.exit(f"{prog}: no residual matrix could be read")

    rows.sort()
    ks = [r[0] for r in rows]
    rmsd = [r[1] for r in rows]
    kn = knee(ks, rmsd)

    hdr = "K\trmsd\tmeanAbs\tmaxAbs\tnonFinite\tdrop\tinMode\telbow"
    lines = [hdr]
    for i, r in enumerate(rows):
        drop = "" if i == 0 or rmsd[i - 1] == 0 else \
            f"{(rmsd[i - 1] - r[1]) / rmsd[i - 1]:.4f}"
        lines.append("\t".join([str(r[0]), f"{r[1]:.6f}", f"{r[2]:.6f}",
                                f"{r[3]:.6f}", str(r[4]), drop,
                                inmode.get(r[0], ""),
                                "<-" if r[0] == kn else ""]))

    outf = a.out or (f"{a.glob}_fit.txt" if a.glob != "*" else "fit.txt")
    with open(outf, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    print(f"\n{outf}: {len(rows)} K", file=sys.stderr)
    if kn is None:
        print(f"{prog}: fewer than 3 K, or a flat curve -- no elbow",
              file=sys.stderr)
    else:
        print(f"{prog}: elbow at K = {kn}", file=sys.stderr)
    if bad:
        print(f"{prog}: {len(bad)} residual matrix/matrices unreadable and "
              f"absent from {outf}:", file=sys.stderr)
        for f in bad:
            print("  " + f, file=sys.stderr)
        sys.exit(1)


if a.fit:
    fit()
    sys.exit(0)


rows = []
for K in range(a.min, a.max + 1):
    f = sorted(glob.glob(f"{a.glob}_k{K}_run*.qopt"),
               key=lambda s: int(re.search(r"run(\d+)", s).group(1)))
    if not f:
        continue
    L = [lik(x) for x in f]

    if K == 1 or len(f) == 1:
        mode = np.zeros(len(f), dtype=int)
    else:
        Q = [np.loadtxt(x, ndmin=2) for x in f]
        m = len(Q)
        S = np.ones((m, m))
        for i in range(m):
            for j in range(i + 1, m):
                S[i, j] = S[j, i] = similarity(Q[i], Q[j])
        mode = avg_linkage(1 - S, 1 - a.thr)

    sizes = np.bincount(mode)
    big = int(np.argmax(sizes))
    idx = [i for i, mm in enumerate(mode) if mm == big]
    i = max(idx, key=lambda j: L[j])

    rows.append([K, len(f), len(sizes),
                 ",".join(str(s) for s in sorted(sizes)[::-1]),
                 int(sizes[big]), round(sizes[big] / len(f), 2), L[i], f[i]])
    print(f"# K={K}: {len(f)} runs, {len(sizes)} mode(s), "
          f"{sizes[big]} in dominant", file=sys.stderr)

if not rows:
    sys.exit(f"no files matched {a.glob}_k*_run*.qopt")

hdr = "K\truns\tnModes\tmodeSizes\tinMode\tprop\tlik\tbest"
print(hdr)
with open(out, "w") as fh:
    fh.write(hdr + "\n")
    for r in rows:
        line = "\t".join(str(x) for x in r)
        print(line)
        fh.write(line + "\n")

if not a.no_copy:
    for r in rows:
        shutil.copy(r[-1], re.sub(r"_run\d+", "_best", r[-1]))
    print(f"\ncopied {len(rows)} files as *_best.qopt", file=sys.stderr)

if a.eval:
    jobs = f"{a.glob}_evalAdmixJobs" if a.glob != "*" else "evalAdmixJobs"
    with open(jobs, "w") as fh:
        for r in rows:
            stem = re.sub(r"\.qopt$", "", r[-1])
            best = re.sub(r"_run\d+", "_best", stem)
            fh.write(f"evalAdmix -beagle {a.eval} -fname {stem}.fopt.gz "
                     f"-qname {stem}.qopt -o {best}.corres -minMaf {a.minMaf:g} "
                     f"-P {a.threads}\n")
    print(f"{jobs}: {len(rows)} jobs", file=sys.stderr)
