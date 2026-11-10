#!/usr/bin/env python3
################################################################################
#Script Name    : cloneTestScore.py
#Description    : Score clone-step variants (cloneTestRun.py) on how well each
#                 distance separates clones: replicate tightness, gap clarity
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import argparse, gzip, os, sys

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

prog = os.path.basename(sys.argv[0])
ap = argparse.ArgumentParser(
    description="For each variant folder made by cloneTestRun.py, builds two distances between samples: ANGSD's "
                "IBS (1 - IBS) and genotype discordance (the share of sites where two samples' genotype calls "
                "differ, among sites both called with probability above --conf), and scores how well each "
                "separates clones. tight = largest distance within a replicate set / median distance between "
                "samples (lower is better); clarity = widest empty stretch in the average-linkage tree above "
                "the replicates / median distance (higher is better). The cut is at the bottom of that stretch.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=f"""examples:
  cd $koko/2bRAD/nwgom/past/angsd
  {prog} -s past -r "PGM016.1 PGM016.2 PGM016.3 PGM076.1 PGM076.2 PGM076.3 PGM119.1 PGM119.2 PGM119.3"
      every variant in cloneTests/, IBS and genotypes at call confidence 0.90, 0.95 and 0.99
  {prog} -s ofav -r "$(awk -F' *= *' '$1=="TECHREPS" {{print $2}}' ../popgen.config)" --conf 0.95
      replicates taken from the config; one confidence level
  {prog} -s past -r "..." -d cloneTests --variants base post2 -o pastCloneTests.tsv
      two variants, and the table also written to a file
  {prog} -s past -r "..." --conf 0.99 --plot pastCloneTests.pdf
      the table, plus one tree per variant and distance with its cut and the replicates coloured

cuts:
  cut  (red in the plots): the bottom of the widest empty stretch above the replicates.
  cut2 (blue): walking up from the replicates' largest distance, the bottom of the first gap
  where the next join is at least --ratio times higher and the gap at least --min-width of the
  median distance; none if no such gap within --max-mult x the replicates' distance.
  cut3 (green): the widest gap whose lower edge is within --max-mult x the replicates' largest
  distance (the widest gap near the replicates; gaps between relatives start above that).
  cut4: treeTrimmer, the replicates' largest distance x --buffer (comRAD's cut for IBS).

reading it:
  compare rows for the same species. A good distance has low tight, high clarity,
  every replicate set in one group, and (for genotypes) most pairs with at least
  --min-shared shared sites. Clone counts that change between rows show what each
  choice would do to the data.
""")
ap.add_argument("-s", "--sp", required=True, help="species prefix used by cloneTestRun.py")
ap.add_argument("-r", "--techreps", required=True, help="technical replicates, space-separated <sample>.<n>")
ap.add_argument("-d", "--dir", default="cloneTests", help="folder of variant runs (default cloneTests)")
ap.add_argument("--variants", nargs="*", help="variants to score (default every folder in --dir with outputs)")
ap.add_argument("-b", "--bams", help="bam list the runs used (default: the -b from the first variant's command.txt)")
ap.add_argument("--conf", type=float, nargs="+", default=[0.90, 0.95, 0.99], help="genotype call confidence levels (default 0.90 0.95 0.99)")
ap.add_argument("--min-shared", type=int, default=1000, help="shared confident sites a pair needs to count as well measured (default 1000)")
ap.add_argument("-o", "--out", help="also write the table here (tab-separated)")
ap.add_argument("--ratio", type=float, default=1.5, help="first-gap cut: the next join must be at least this many times higher (default 1.5)")
ap.add_argument("--min-width", type=float, default=0.02, help="first-gap cut: and the gap at least this share of the median distance (default 0.02)")
ap.add_argument("--max-mult", type=float, default=3.0, help="first-gap cut: look no higher than this many times the replicates' largest distance (default 3)")
ap.add_argument("--buffer", type=float, default=1.10, help="treeTrimmer cut (cut4): the replicates' largest distance x this (default 1.10)")
ap.add_argument("--truth", metavar="TSV", help="simulated truth (simCloneKin.py): sample1, sample2, relation for every related pair; "
                "adds, per cut, the true clone pairs grouped and the related or unrelated pairs wrongly grouped")
ap.add_argument("--plot", metavar="PDF", help="also draw every tree (one page per variant and distance) into this PDF")
a = ap.parse_args()

variants = a.variants or sorted(v for v in os.listdir(a.dir) if os.path.exists(os.path.join(a.dir, v, a.sp + ".geno.gz")))
if not variants:
    sys.exit(f"{prog}: no {a.dir}/<variant>/{a.sp}.geno.gz yet")
bams = a.bams
if not bams:
    cmd = open(os.path.join(a.dir, variants[0], "command.txt")).read().split()
    bams = cmd[cmd.index("-b") + 1]
names = [os.path.basename(l.split()[0]).removesuffix(".bam") for l in open(bams) if l.strip()]
n = len(names)
pos = {s: i for i, s in enumerate(names)}
sets = {}
for r in a.techreps.split():
    if r in pos:
        sets.setdefault(r.split(".")[0], []).append(pos[r])
sets = {k: v for k, v in sets.items() if len(v) >= 2}
if not sets:
    sys.exit(f"{prog}: no replicate set has two or more samples in {bams}")
iu = np.triu_indices(n, 1)
truth = {}
if a.truth:
    with open(a.truth) as fh:
        next(fh)
        for l in fh:
            x, y, r = l.rstrip("\n").split("\t")
            if x in pos and y in pos:
                truth[frozenset((x, y))] = r
nclone = sum(1 for r in truth.values() if r == "clone")
SHORT = {"parent-offspring": "PO", "full-sibs": "FS", "half-sibs": "HS", "selfed-sibs": "SS", "unrelated": "U"}


def judge(cl):
    """true clone pairs grouped / all true clone pairs, and the non-clone pairs grouped, by relation"""
    if not truth:
        return ""
    from collections import Counter
    got, wrong = 0, Counter()
    for g in np.unique(cl):
        mem = [names[i] for i in np.where(cl == g)[0]]
        for x in range(len(mem)):
            for y in range(x + 1, len(mem)):
                r = truth.get(frozenset((mem[x], mem[y])), "unrelated")
                if r == "clone":
                    got += 1
                else:
                    wrong[r] += 1
    w = " ".join(f"{SHORT.get(k, k)}{v}" for k, v in sorted(wrong.items()))
    return f"{got}/{nclone}" + (f" wrong:{w}" if w else " wrong:0")


def score(M):
    """replicate tightness, the widest empty stretch above the replicates, and the clones at its bottom"""
    scale = float(np.median(M[iu]))
    if not scale > 0:
        return dict(scale=scale, repmax=float("nan"), tight=float("nan"), gaplo=float("nan"), gaphi=float("nan"),
                    clarity=float("nan"), cut=float("nan"), groups="NA", dropped="NA", reps="no data",
                    cut2=float("nan"), groups2="NA", dropped2="NA", cut3=float("nan"), groups3="NA", dropped3="NA",
                    truth1="", truth2="", truth3="", cut4=float("nan"), groups4="NA", dropped4="NA", truth4="")
    repmax = max(M[i, j] for v in sets.values() for i in v for j in v if i < j)
    Z = linkage(squareform(M, checks=False), "average")
    hs = sorted(h for h in np.unique(Z[:, 2]) if repmax < h < scale)
    edges = [repmax] + hs
    gaps = [(edges[k + 1] - edges[k], edges[k], edges[k + 1]) for k in range(len(edges) - 1)]
    g, lo, hi = max(gaps) if gaps else (0.0, repmax, repmax)
    cut = lo + 1e-9
    cl = fcluster(Z, t=cut, criterion="distance")
    ev1 = judge(cl)
    sizes = np.bincount(cl)
    groups = int((sizes > 1).sum())
    dropped = int((sizes[sizes > 1] - 1).sum())
    together = sum(len(set(cl[v])) == 1 for v in sets.values())
    # first-gap cut: walking up from the replicates, the first join at least --ratio times higher than the one
    # below it (and a gap wider than --min-width of the median), no higher than --max-mult x the replicates
    up = [repmax] + sorted(h for h in np.unique(Z[:, 2]) if h > repmax)
    cut2 = None
    for k in range(len(up) - 1):
        lo2, hi2 = up[k], up[k + 1]
        if lo2 > a.max_mult * repmax:
            break
        if hi2 >= a.ratio * max(lo2, 1e-12) and hi2 - lo2 >= a.min_width * scale:
            cut2 = lo2 + 1e-9
            break
    # windowed widest gap: the widest gap between consecutive joins whose lower edge is no higher than
    # --max-mult x the replicates' largest distance (gaps between relatives start far above that)
    win = [(up[k + 1] - up[k], up[k]) for k in range(len(up) - 1) if up[k] <= a.max_mult * repmax]
    cut3 = (max(win)[1] + 1e-9) if win else None
    ev2 = ev3 = ""
    # treeTrimmer: the replicates' largest distance x --buffer (the cut comRAD uses for IBS)
    cut4 = repmax * a.buffer
    cl4 = fcluster(Z, t=cut4, criterion="distance"); sz4 = np.bincount(cl4)
    groups4, dropped4, ev4 = int((sz4 > 1).sum()), int((sz4[sz4 > 1] - 1).sum()), judge(cl4)
    if cut3 is None:
        groups3 = dropped3 = "none"
    else:
        cl3 = fcluster(Z, t=cut3, criterion="distance"); sz3 = np.bincount(cl3); ev3 = judge(cl3)
        groups3, dropped3 = int((sz3 > 1).sum()), int((sz3[sz3 > 1] - 1).sum())
    if cut2 is None:
        groups2 = dropped2 = "none"
    else:
        cl2 = fcluster(Z, t=cut2, criterion="distance"); sz2 = np.bincount(cl2); ev2 = judge(cl2)
        groups2, dropped2 = int((sz2 > 1).sum()), int((sz2[sz2 > 1] - 1).sum())
    return dict(scale=scale, repmax=repmax, tight=repmax / scale, gaplo=lo, gaphi=hi, clarity=g / scale,
                cut=cut, groups=groups, dropped=dropped, reps=f"{together}/{len(sets)}", Z=Z,
                cut2=cut2 if cut2 is not None else float("nan"), groups2=groups2, dropped2=dropped2,
                cut3=cut3 if cut3 is not None else float("nan"), groups3=groups3, dropped3=dropped3,
                truth1=ev1, truth2=ev2, truth3=ev3, cut4=cut4, groups4=groups4, dropped4=dropped4, truth4=ev4)


def read_geno(path):
    """called genotype (0/1/2) and its probability per site and sample; cached next to the file as .npz"""
    cache = path[:-len(".geno.gz")] + ".genocache.npz"
    if os.path.exists(cache) and os.path.getmtime(cache) >= os.path.getmtime(path):
        z = np.load(cache)
        if z["G"].shape[1] == n:
            return z["G"], z["P"]
    G, P = [], []
    with gzip.open(path, "rt") as fh:
        for l in fh:
            q = np.fromstring(l.split("\t", 2)[2], sep="\t").reshape(n, 3)
            G.append(q.argmax(1).astype(np.int8)); P.append(q.max(1).astype(np.float32))
    G = np.array(G, dtype=np.int8); P = np.array(P, dtype=np.float32)
    try:
        np.savez(cache, G=G, P=P)
    except OSError:
        pass
    return G, P


def discordance(G, P, c):
    """share of sites with different calls, among sites both samples called with probability > c;
    counts by matrix products (shared = C'C, matching = sum over genotypes of (C & G==g)'(C & G==g))"""
    C = (P > c)
    Cf = C.astype(np.float32)
    S = Cf.T @ Cf
    same = np.zeros_like(S)
    for g in (0, 1, 2):
        X = (C & (G == g)).astype(np.float32)
        same += X.T @ X
    S = np.rint(S).astype(np.int64); same = np.rint(same).astype(np.int64)
    D = (S - same) / np.maximum(S, 1)
    D = (D + D.T) / 2; np.fill_diagonal(D, 0)
    return D, S


rows = []
for v in variants:
    p = os.path.join(a.dir, v, a.sp)
    if not (os.path.exists(p + ".ibsMat") and os.path.exists(p + ".geno.gz")):
        print(f"{v}: outputs missing, skipped", file=sys.stderr)
        continue
    M = np.loadtxt(p + ".ibsMat", ndmin=2) if os.path.getsize(p + ".ibsMat") else np.empty((0, 0))
    if M.shape != (n, n):
        print(f"{v}: {p}.ibsMat is {M.shape[0]} x {M.shape[1] if M.ndim > 1 else 0}, not {n} x {n} (still running, or failed: see its log); skipped", file=sys.stderr)
        continue
    rows.append(dict(variant=v, metric="ibs", sites="", shared="", wellMeasured="", **score(M)))
    G, P = read_geno(p + ".geno.gz")
    for c in a.conf:
        D, S = discordance(G, P, c)
        sh = S[iu]
        rows.append(dict(variant=v, metric=f"geno{c:g}", sites=len(G), shared=int(np.median(sh)),
                         wellMeasured=f"{(sh >= a.min_shared).mean():.2f}", **score(D)))

cols = ["variant", "metric", "sites", "shared", "wellMeasured", "repmax", "scale", "tight", "gaplo", "gaphi",
        "clarity", "cut", "reps", "groups", "dropped", "cut2", "groups2", "dropped2", "cut3", "groups3", "dropped3",
        "cut4", "groups4", "dropped4"]
if a.truth:
    cols += ["truth1", "truth2", "truth3", "truth4"]
fmt = lambda x: f"{x:.4f}" if isinstance(x, float) else str(x)
w = {c: max(len(c), *(len(fmt(r[c])) for r in rows)) for c in cols}
print("  ".join(c.rjust(w[c]) for c in cols))
for r in rows:
    print("  ".join(fmt(r[c]).rjust(w[c]) for c in cols))
if a.out:
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n" + "".join("\t".join(fmt(r[c]) for c in cols) + "\n" for r in rows))
ok = [r for r in rows if "/" in r["reps"]]
best = max(ok, key=lambda r: (r["reps"].split("/")[0] == r["reps"].split("/")[1], r["clarity"] - r["tight"]))
print(f"\nclearest separation: {best['variant']} / {best['metric']} (tight {best['tight']:.3f}, clarity {best['clarity']:.3f}, "
      f"cut {best['cut']:.4f}: {best['groups']} groups, {best['dropped']} dropped)")

if a.plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from scipy.cluster.hierarchy import dendrogram
    setOf = {names[i]: k for k, v in sets.items() for i in v}
    pal = plt.get_cmap("tab10")
    colour = {k: pal(j % 10) for j, k in enumerate(sorted(sets))}
    with PdfPages(a.plot) as pdf:
        for r in ok:
            fig, ax = plt.subplots(figsize=(max(12, n * 0.12), 7))
            dn = dendrogram(r["Z"], labels=names, ax=ax, color_threshold=0, above_threshold_color="black", leaf_font_size=5)
            for t in ax.get_xticklabels():
                if t.get_text() in setOf:
                    t.set_color(colour[setOf[t.get_text()]]); t.set_fontweight("bold")
            ax.axhline(r["cut"], color="red", ls=(0, (5, 3)), lw=.8)
            if r["cut2"] == r["cut2"]:
                ax.axhline(r["cut2"], color="blue", ls=(0, (2, 2)), lw=.8)
            if r["cut3"] == r["cut3"]:
                ax.axhline(r["cut3"], color="green", ls="-", lw=.8, alpha=.7)
            ax.axhspan(r["gaplo"], r["gaphi"], color="red", alpha=.07, lw=0)
            label = "1 - IBS" if r["metric"] == "ibs" else f"genotype distance (calls > {r['metric'][4:]})"
            ax.set_ylabel(label)
            ax.set_title(f"{a.sp}  {r['variant']} / {r['metric']}:  tight {r['tight']:.3f}, clarity {r['clarity']:.3f}, "
                         f"widest-gap cut (red) {r['cut']:.4f} -> {r['groups']} groups, {r['dropped']} dropped; first-gap cut (blue) "
                         + (f"{r['cut2']:.4f} -> {r['groups2']} groups, {r['dropped2']} dropped" if r["cut2"] == r["cut2"] else "none")
                         + (f"; windowed cut (green) {r['cut3']:.4f} -> {r['groups3']} groups, {r['dropped3']} dropped" if r["cut3"] == r["cut3"] else "")
                         + f"; replicate sets together {r['reps']}"
                         + (f"; {r['sites']} sites" if r["sites"] != "" else ""), fontsize=9)
            for sp_ in ("top", "right"):
                ax.spines[sp_].set_visible(False)
            fig.tight_layout()
            pdf.savefig(fig); plt.close(fig)
    print(f"trees: {a.plot} ({len(ok)} pages)")
