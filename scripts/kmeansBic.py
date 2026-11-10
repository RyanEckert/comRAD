#!/usr/bin/env python3
################################################################################
#Script Name    : kmeansBic.py
#Description    : K from k-means clustering of PCA scores with BIC (the
#                 find.clusters step of DAPC), from a PCAngsd covariance matrix
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import argparse, os, sys

import numpy as np

prog = os.path.basename(sys.argv[0])
ap = argparse.ArgumentParser(
    description="Clusters samples by k-means on their principal components and scores each K by BIC, as adegenet's "
                "find.clusters does before DAPC (Jombart et al. 2010): BIC = n ln(WSS / n) + K ln(n), with WSS the "
                "within-cluster sum of squares. The suggested K is the lowest BIC (adegenet's 'min' criterion). The PCs come from a PCAngsd "
                "covariance matrix, so genotype uncertainty is accounted for as in the PCA. An independent read on K: "
                "like PCAngsd's own count, it responds to any discrete grouping, families and batches included.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=f"""examples:
  {prog} -c ../pcangsd/ofavPcangsd.cov -k 7 -o ofav_bic.txt
      K = 1..7 on the PCs that explain 80% of the variance (at most n/3 of them)
  {prog} -c ../pcangsd/ofavNeutralPcangsd.cov -k 7 -o ofavNeutral_bic.txt --var 0.9
      the kin-free, neutral run, on PCs explaining 90%
  {prog} -c ../pcangsd/sintPcangsd.cov -k 10 --pcs 20 -o sint_bic.txt
      a fixed number of PCs

output:
  K, BIC, WSS, and '<-' on the lowest BIC; a comment line gives the PCs used.

tested on simulated genotypes with noise (50 per population, 3000 SNPs): one population -> 1; two at FST 0.1
or 0.02 -> 2; three at FST 0.05 -> 3; three at FST 0.01 -> 1 (too weak). Clone groups and close families
also form clusters, so run it on clone-free data and read it beside the kin-free rerun.
""")
ap.add_argument("-c", "--cov", required=True, help="PCAngsd covariance matrix (.cov)")
ap.add_argument("-k", "--maxK", type=int, required=True, help="largest K to test (from 1)")
ap.add_argument("--var", type=float, default=0.8, help="keep the PCs that explain this share of the variance (default 0.8)")
ap.add_argument("--pcs", type=int, help="keep this many PCs instead")
ap.add_argument("--starts", type=int, default=50, help="k-means random starts per K; the best is kept (default 50)")
ap.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
ap.add_argument("-o", "--out", required=True, help="output table")
a = ap.parse_args()

C = np.loadtxt(a.cov, ndmin=2)
n = C.shape[0]
val, vec = np.linalg.eigh((C + C.T) / 2)
o = np.argsort(val)[::-1]
val, vec = val[o], vec[:, o]
pos = val > 0
val, vec = val[pos], vec[:, pos]
if a.pcs:
    m = min(a.pcs, len(val))
else:
    m = int(np.searchsorted(np.cumsum(val) / val.sum(), a.var) + 1)
    m = max(1, min(m, max(1, n // 3)))
X = vec[:, :m] * np.sqrt(val[:m])                      # PC scores
rng = np.random.default_rng(a.seed)


def kmeans(X, k):
    """Lloyd's algorithm from k-means++ starts; the lowest within-cluster sum of squares over --starts runs"""
    best = np.inf
    for _ in range(a.starts):
        c = [X[rng.integers(len(X))]]
        for _ in range(1, k):
            d = np.min([((X - ci) ** 2).sum(1) for ci in c], axis=0)
            c.append(X[rng.choice(len(X), p=d / d.sum()) if d.sum() > 0 else rng.integers(len(X))])
        c = np.array(c)
        for _ in range(300):
            lab = np.argmin(((X[:, None, :] - c[None]) ** 2).sum(2), axis=1)
            new = np.array([X[lab == j].mean(0) if (lab == j).any() else c[j] for j in range(k)])
            if np.allclose(new, c):
                break
            c = new
        wss = ((X - c[lab]) ** 2).sum()
        best = min(best, wss)
    return best


rows = []
for k in range(1, min(a.maxK, n) + 1):
    w = kmeans(X, k)
    rows.append((k, n * np.log(max(w, 1e-300) / n) + k * np.log(n), w))     # as adegenet::find.clusters
kbest = min(rows, key=lambda r: r[1])[0]
with open(a.out, "w") as fh:
    fh.write(f"# pcs {m} of {len(val)} ({100 * val[:m].sum() / val.sum():.1f}% of the variance), n = {n}\n")
    fh.write("K\tBIC\tWSS\tbest\n")
    for k, b, w in rows:
        fh.write(f"{k}\t{b:.3f}\t{w:.6g}\t{'<-' if k == kbest else ''}\n")
print(f"k-means BIC: K = {kbest} on {m} PCs ({n} samples); table in {a.out}")
