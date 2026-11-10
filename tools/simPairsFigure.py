#!/usr/bin/env python3
################################################################################
#Script Name    : simPairsFigure.py
#Description    : The docs figure of every simulated pair's IBS and genotype
#                 distance, by true relationship, with the clone cuts
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import argparse, gzip, os, sys

import numpy as np

prog = os.path.basename(sys.argv[0])
ap = argparse.ArgumentParser(
    description="Draws docs/examples/img/simPairs.png: each pair of samples in a simCloneKin.py run at its IBS and "
                "genotype distance, sorted by true relationship, with the cuts the clone step chose.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=f"""example (the figure in docs/examples/clone-simulation.md):
  simCloneKin.py -o simP --big-clone 25 --fixed-het 0.3
  cd simP
  cloneGroups.py -i cloneTests/base/sim.ibsMat -g cloneTests/base/sim.geno.gz -b simBams -r "$(cat techreps.txt)" --metric ibs -o ibs --dir cut
  cloneGroups.py -i cloneTests/base/sim.ibsMat -g cloneTests/base/sim.geno.gz -b simBams -r "$(cat techreps.txt)" -o geno --dir cut
  {prog} --ibs-cut $(awk '$1=="cut"{{print $2}}' cut/ibsCloneMetric.txt) --geno-cut $(awk '$1=="cut"{{print $2}}' cut/genoCloneMetric.txt) \\
    --conf 0.95 -o ../comRAD/docs/examples/img/simPairs.png
""")
ap.add_argument("--dir", default=".", help="the simCloneKin.py output folder (default here)")
ap.add_argument("--ibs-cut", type=float, required=True, help="the IBS cut, from CloneMetric.txt")
ap.add_argument("--geno-cut", type=float, required=True, help="the genotype-distance cut, from CloneMetric.txt")
ap.add_argument("--conf", type=float, default=0.95, help="genotype call confidence the cut was made at (default 0.95)")
ap.add_argument("-o", "--out", default="simPairs.png", help="output PNG (default simPairs.png)")
a = ap.parse_args()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = lambda *x: os.path.join(a.dir, *x)
names = [l.split()[0].split("/")[-1].replace(".bam", "") for l in open(d("simBams")) if l.strip()]
n = len(names)
pos = {s: i for i, s in enumerate(names)}
M = np.loadtxt(d("cloneTests", "base", "sim.ibsMat"))
G, P = [], []
for l in gzip.open(d("cloneTests", "base", "sim.geno.gz"), "rt"):
    q = np.array(l.split()[2:], float).reshape(n, 3)
    G.append(q.argmax(1)); P.append(q.max(1))
G, P = np.array(G), np.array(P)
C = P > a.conf
S = C.T.astype(np.float32) @ C.astype(np.float32)
same = sum(((C & (G == g)).astype(np.float32)).T @ (C & (G == g)).astype(np.float32) for g in (0, 1, 2))
D = (S - same) / np.maximum(S, 1)
rel = {}
for l in list(open(d("truth.tsv")))[1:]:
    x, y, r = l.split()
    rel[(pos[x], pos[y])] = rel[(pos[y], pos[x])] = r
big = {i for i, s in enumerate(names) if "big" in s}
cats = ["clone (big heterozygous genet)", "clone / replicate (others)", "parent-offspring", "selfed sibs",
        "full sibs", "half sibs", "unrelated"]


def cat(i, j):
    r = rel.get((i, j))
    if r == "clone":
        return cats[0] if i in big and j in big else cats[1]
    return {"parent-offspring": cats[2], "selfed-sibs": cats[3], "full-sibs": cats[4], "half-sibs": cats[5]}.get(r, cats[6])


iu = np.triu_indices(n, 1)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
rng = np.random.default_rng(1)
cols = ["#d95f02", "#e6ab02", "#7570b3", "#e7298a", "#66a61e", "#1b9e77", "#999999"]
for ax, (mat, lab, cut, title) in zip(axs, [(M, "1 - IBS", a.ibs_cut, "IBS, treeTrimmer cut"),
                                             (D, f"genotype distance (calls > {a.conf:g})", a.geno_cut,
                                              "genotype distance, widest-gap cut")]):
    for k, c in enumerate(cats):
        v = [mat[i, j] for i, j in zip(*iu) if cat(i, j) == c]
        if c == "unrelated" and len(v) > 600:
            v = list(rng.choice(v, 600, replace=False))
        ax.scatter(v, k + rng.uniform(-.25, .25, len(v)), s=4, color=cols[k], alpha=.6, linewidth=0)
    ax.axvline(cut, color="red", ls="--", lw=1)
    ax.set_xlabel(lab)
    ax.set_title(title, loc="left", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
axs[0].set_yticks(range(len(cats)))
axs[0].set_yticklabels(cats)
axs[0].invert_yaxis()
fig.tight_layout()
fig.savefig(a.out, dpi=110)
print(f"wrote {a.out}")
