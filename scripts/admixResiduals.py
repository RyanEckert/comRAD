#!/usr/bin/env python3
################################################################################
#Script Name    : admixResiduals.py
#Description    : Summarise evalAdmix residual correlations per K by blocks of
#                 samples, suggest the smallest K with no block left, plot them
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import argparse, csv, math, os, re, sys

import numpy as np

prog = os.path.basename(sys.argv[0])
ap = argparse.ArgumentParser(
    description="evalAdmix leaves positive residual correlations among individuals that share something the "
                "admixture model at that K does not capture. This puts each individual in its main cluster at "
                "--maxK (the finest split NGSadmix found) and averages the residuals within and between those "
                "blocks for every K; a K that lumps two real lineages leaves positive blocks within them. The "
                "suggested K is the smallest whose worst cluster-block mean |residual| is under --cut. Sampling "
                "groups from the metadata are summarised too, as an independent check.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=f"""examples:
  {prog} -s ofav -k 7 -b ../angsd/ofavBams
      reads ofav_k1_best.corres .. ofav_k7_best.corres and ofav_k7_best.qopt here; writes
      ofav_resid.txt, ofav_residBlocks.txt and ofavResiduals.pdf
  {prog} -s ofav -k 7 -b ../angsd/ofavBams --cut 0.02 \\
      --meta ../metadata.csv --meta-sample sample --meta-groups bank depthZone
      a looser cut-off, and sampling groups (bank x depthZone) as the second check

outputs:
  <sp>_resid.txt       K, worstCluster, worstClusterBlock, worstGroup, residK ('<-' marks the suggested K)
  <sp>_residBlocks.txt K, a, b, within, mean, n: every cluster block with >= 3 samples on each side
  <sp>Residuals.pdf    one residual heatmap per K, then the block means for all K
""")
ap.add_argument("-s", "--sp", required=True, help="prefix of the NGSadmix/evalAdmix files (<sp>_k<K>_best.*)")
ap.add_argument("-k", "--maxK", type=int, required=True, help="largest K; its best .qopt sets the blocks")
ap.add_argument("-b", "--bams", required=True, help="bam list, in the order the beagle had the samples")
ap.add_argument("--cut", type=float, default=0.01, help="largest worst-block mean |residual| that counts as fitted (default 0.01)")
ap.add_argument("--meta", help="sample metadata CSV, for the sampling-group check")
ap.add_argument("--meta-sample", default="sample", help="metadata column with the sample names (default sample)")
ap.add_argument("--meta-groups", nargs="*", default=[], help="metadata columns whose combination is a sampling group")
a = ap.parse_args()

sp, maxK = a.sp, a.maxK
ids = [re.sub(r"\.bam$", "", os.path.basename(l.split()[0])) for l in open(a.bams) if l.strip()]
n = len(ids)
grp = a.meta_groups
if grp and a.meta:
    with open(a.meta, newline="", encoding="utf-8-sig") as fh:
        meta = {r[a.meta_sample]: r for r in csv.DictReader(fh)}
    grpOf = ["_".join(meta.get(s, {}).get(g) or "NA" for g in grp) for s in ids]
else:
    grpOf = ["all"] * n


def read(path):
    rows = []
    for l in open(path):
        if l.strip():
            rows.append([float("nan") if x in ("NA", "nan", "NaN", "-nan") else float(x) for x in l.split()])
    return np.array(rows)


qMax = read(f"{sp}_k{maxK}_best.qopt")
if qMax.shape[0] != n:
    sys.exit(f"{prog}: {sp}_k{maxK}_best.qopt has {qMax.shape[0]} rows; {a.bams} has {n}")
blk = [f"C{int(np.argmax(r)) + 1}" for r in qMax]
ordr = sorted(range(n), key=lambda i: (blk[i], grpOf[i]))


def uniq(x):
    return list(dict.fromkeys(x))


def blockMeans(r, lab, minN=3):
    out = []
    for A in uniq(lab):
        for B in uniq(lab):
            if A > B:
                continue
            ia = [i for i, x in enumerate(lab) if x == A]
            ib = [i for i, x in enumerate(lab) if x == B]
            if len(ia) < minN or len(ib) < minN:
                continue
            if A == B:
                s = r[np.ix_(ia, ia)]
                v = s[np.tril_indices(len(ia), -1)]
            else:
                v = r[np.ix_(ia, ib)].ravel()
            v = v[np.isfinite(v)]
            if len(v) >= 3:
                out.append(dict(a=A, b=B, within=A == B, mean=float(v.mean()), n=len(v)))
    return out


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap

cmap = LinearSegmentedColormap.from_list("resid", ["#2166ac", "white", "#b2182b"])
cmap.set_bad("#cccccc")
rows, bm = [], []
pdf = PdfPages(f"{sp}Residuals.pdf")
for K in range(1, maxK + 1):
    f = f"{sp}_k{K}_best.corres"
    if not os.path.exists(f):
        continue
    r = read(f)
    if r.shape != (n, n):
        sys.exit(f"{prog}: {f} is {r.shape[0]} x {r.shape[1]}; the bamlist has {n}")
    b1, b2 = blockMeans(r, blk), blockMeans(r, grpOf)
    w1 = max(b1, key=lambda x: abs(x["mean"])) if b1 else None      # first of any ties, as R's which.max
    rows.append(dict(K=K,
                     worstCluster=round(abs(w1["mean"]), 5) if b1 else None,
                     worstClusterBlock=f"{w1['a']}-{w1['b']}" if b1 else "",
                     worstGroup=round(max(abs(x["mean"]) for x in b2), 5) if b2 else None))
    bm += [dict(K=K, **x) for x in b1]
    o = np.array(ordr)
    m = np.clip(r[np.ix_(o, o)], -.1, .1)
    np.fill_diagonal(m, np.nan)
    fig, ax = plt.subplots(figsize=(10, 9))
    im = ax.imshow(np.ma.masked_invalid(m), cmap=cmap, vmin=-.1, vmax=.1, interpolation="nearest")
    sizes = [sum(1 for i in ordr if blk[i] == c) for c in uniq(blk[i] for i in ordr)]
    for e in np.cumsum(sizes)[:-1]:
        ax.axvline(e - .5, color="black", lw=.3)
        ax.axhline(e - .5, color="black", lw=.3)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{sp} K = {K}: evalAdmix residual correlations\n"
                 f"ordered by main cluster at K = {maxK} (lines), then {' x '.join(grp) if grp else 'sample'}; "
                 f"worst cluster-block mean |r| = {rows[-1]['worstCluster'] if b1 else 'NA'}", fontsize=10)
    fig.colorbar(im, ax=ax, label="residual\ncorrelation", shrink=.6)
    pdf.savefig(fig); plt.close(fig)
if bm:
    Ks = uniq(x["K"] for x in bm)
    cols = math.ceil(math.sqrt(len(Ks)))
    fig, axs = plt.subplots(math.ceil(len(Ks) / cols), cols, figsize=(4 * cols, 4 * math.ceil(len(Ks) / cols)), squeeze=False)
    for ax, K in zip(axs.ravel(), Ks):
        xs = [x for x in bm if x["K"] == K]
        labs = sorted(set([x["a"] for x in xs] + [x["b"] for x in xs]))
        g = np.full((len(labs), len(labs)), np.nan)
        for x in xs:                                 # both halves of each block pair
            g[labs.index(x["a"]), labs.index(x["b"])] = g[labs.index(x["b"]), labs.index(x["a"])] = x["mean"]
        ax.imshow(np.ma.masked_invalid(np.clip(g, -.05, .05)), cmap=cmap, vmin=-.05, vmax=.05)
        for i in range(len(labs)):
            for j in range(len(labs)):
                if np.isfinite(g[i, j]):
                    ax.text(j, i, f"{g[i, j]:.3f}", ha="center", va="center", fontsize=6)
        ax.set_xticks(range(len(labs)), labs, rotation=90, fontsize=7); ax.set_yticks(range(len(labs)), labs, fontsize=7)
        ax.set_title(f"K: {K}", fontsize=9)
    for ax in axs.ravel()[len(Ks):]:
        ax.axis("off")
    fig.suptitle(f"Mean residual correlation per block (main cluster at K = {maxK})")
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)
pdf.close()
if not rows:
    sys.exit(f"{prog}: no {sp}_k*_best.corres to read")

ok = [x["K"] for x in rows if x["worstCluster"] is not None and x["worstCluster"] < a.cut]
rk = min(ok) if ok else None
fmt = lambda v: "NA" if v is None else ("TRUE" if v is True else "FALSE" if v is False else f"{v:.15g}" if isinstance(v, float) else str(v))
cols = ["K", "worstCluster", "worstClusterBlock", "worstGroup", "residK"]
for x in rows:
    x["residK"] = "<-" if x["K"] == rk else ""
txt = "\t".join(cols) + "\n" + "".join("\t".join(fmt(x[c]) for c in cols) + "\n" for x in rows)
open(f"{sp}_resid.txt", "w").write(txt)
if bm:
    bc = ["K", "a", "b", "within", "mean", "n"]
    open(f"{sp}_residBlocks.txt", "w").write("\t".join(bc) + "\n" + "".join("\t".join(fmt(x[c]) for c in bc) + "\n" for x in bm))
print(txt, end="")
