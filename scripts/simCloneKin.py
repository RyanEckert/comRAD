#!/usr/bin/env python3
################################################################################
#Script Name    : simCloneKin.py
#Description    : Simulate low-coverage RAD genotypes for a pedigree with clones,
#                 technical replicates and kin, in the format cloneTestScore.py reads
#Author         : Ryan J Eckert
#Email          : ryan.j.eckert@gmail.com
################################################################################
import argparse, gzip, os, sys
from itertools import combinations

import numpy as np

prog = os.path.basename(sys.argv[0])
ap = argparse.ArgumentParser(
    description="Simulates samples from two lineages, each with unrelated founders, full-sib families whose "
                "parents are also sampled (parent-offspring pairs), half-sibs, selfed sibs, clones (extra ramets "
                "of a genet) and technical replicates (the same colony sequenced again). Reads are drawn per "
                "sample and site with the sample's own coverage and a sequencing error rate; genotype "
                "probabilities use a population-frequency prior (as ANGSD -doPost 1) and IBS uses one random read "
                "per sample and site (as -doIBS 1). Writes <out>/<variant>/<sp>.geno.gz, .ibsMat and command.txt, "
                "the bam list, and a truth table, so cloneTestScore.py can score it with --truth.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=f"""examples:
  {prog} -o simTest
      the default design (about 125 samples, 8000 SNPs) in simTest/cloneTests/base; then:
  cd simTest && cloneTestScore.py -s sim -r "$(cat techreps.txt)" --truth truth.tsv --conf 0.95 0.99 --plot sim.pdf
  {prog} -o simPast --big-clone 25 --fixed-het 0.3
      lineage B gets one genet with 25 ramets and extra fixed heterozygosity (past's comb)
  {prog} -o simLow --cov 2.5
      everything at about 2.5x: where genotype calls run short and IBS should do better
  {prog} -o simTest --variant noisy --err 0.02
      a second run of the same design (same samples) with more sequencing error, scored next to base;
      a different design (counts, --big-clone, --seed) needs its own --outdir

design (per lineage; --scale multiplies the counts):
  12 unrelated founders; 4 full-sib families (2 parents sampled + 3 offspring);
  2 half-sib groups (1 parent, 2 mates, 2 offspring each); 1 selfed family (parent + 3 selfed offspring);
  4 genets with 1-3 extra clone ramets; 2 colonies sequenced 3 times as technical replicates
""")
ap.add_argument("-o", "--outdir", required=True, help="folder for the simulation (created)")
ap.add_argument("--variant", default="base", help="name of this run under <outdir> (default base)")
ap.add_argument("-s", "--sp", default="sim", help="file prefix (default sim)")
ap.add_argument("--sites", type=int, default=8000, help="SNPs (default 8000)")
ap.add_argument("--fst", type=float, default=0.1, help="differentiation between the two lineages (default 0.1)")
ap.add_argument("--cov", type=float, default=7.0, help="median per-sample coverage (default 7)")
ap.add_argument("--cov-sd", type=float, default=0.5, help="spread of coverage between samples, log scale (default 0.5)")
ap.add_argument("--weak", type=int, default=6, help="extra weak samples per lineage, at about 2x, all clones or replicates of others (default 6)")
ap.add_argument("--err", type=float, default=0.005, help="sequencing error per base (default 0.005)")
ap.add_argument("--somatic", type=float, default=2e-4, help="somatic mutation rate between clone ramets, per site (default 2e-4)")
ap.add_argument("--big-clone", type=int, default=0, help="lineage B: one genet with this many ramets (default 0)")
ap.add_argument("--fixed-het", type=float, default=0.0, help="with --big-clone: share of sites forced heterozygous in that genet (default 0)")
ap.add_argument("--scale", type=int, default=1, help="multiply the design counts (default 1)")
ap.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
S = a.sites

# allele frequencies: an ancestral frequency, then one per lineage (Balding-Nichols)
p0 = np.clip(rng.beta(0.6, 0.6, S), 0.02, 0.98)
k = (1 - a.fst) / a.fst
pl = {L: np.clip(rng.beta(p0 * k, (1 - p0) * k), 0.005, 0.995) for L in "AB"}

genets, samples = {}, []          # genet -> (hap S x 2, lineage, mother, father, kind); sample = (name, genet, role)


def founder(L):
    return (rng.random((S, 2)) < pl[L][:, None]).astype(np.int8)


def child(m, f):
    hm, hf = genets[m][0], genets[f][0]
    return np.stack([hm[np.arange(S), rng.integers(0, 2, S)], hf[np.arange(S), rng.integers(0, 2, S)]], 1).astype(np.int8)


def add_genet(name, hap, L, mother=None, father=None):
    genets[name] = (hap, L, mother, father)


def add_sample(name, g, role):
    samples.append((name, g, role))


for L in "AB":
    c = 0
    def nm(tag):
        global c
        c += 1
        return f"{L}{tag}{c:03d}"
    for _ in range(12 * a.scale):
        g = nm("u"); add_genet(g, founder(L), L); add_sample(g, g, "founder")
    for _ in range(4 * a.scale):                     # full-sib families with both parents sampled
        m, f = nm("p"), nm("p"); add_genet(m, founder(L), L); add_genet(f, founder(L), L)
        add_sample(m, m, "parent"); add_sample(f, f, "parent")
        for _ in range(3):
            o = nm("o"); add_genet(o, child(m, f), L, m, f); add_sample(o, o, "offspring")
    for _ in range(2 * a.scale):                     # half-sibs: one parent (not sampled), two mates
        p = nm("h"); add_genet(p, founder(L), L)
        for _ in range(2):
            mate = nm("h"); add_genet(mate, founder(L), L)
            for _ in range(2):
                o = nm("o"); add_genet(o, child(p, mate), L, p, mate); add_sample(o, o, "offspring")
    for _ in range(1 * a.scale):                     # selfed family
        p = nm("s"); add_genet(p, founder(L), L); add_sample(p, p, "parent")
        for _ in range(3):
            o = nm("o"); add_genet(o, child(p, p), L, p, p); add_sample(o, o, "offspring")
    pool = [s[1] for s in samples if genets[s[1]][1] == L]
    for g in rng.choice(pool, 4 * a.scale, replace=False):     # clones: extra ramets
        for r in range(int(rng.integers(1, 4))):
            add_sample(f"{g}r{r + 1}", g, "ramet")
    for g in rng.choice([x for x in pool if not any(s[0].startswith(x + "r") for s in samples)], 2 * a.scale, replace=False):
        idx = [i for i, s in enumerate(samples) if s[0] == g][0]
        samples[idx] = (f"{g}.1", g, "techrep")
        add_sample(f"{g}.2", g, "techrep"); add_sample(f"{g}.3", g, "techrep")
    clonable = [s[1] for s in samples if genets[s[1]][1] == L and s[2] in ("ramet", "techrep")]
    for w in range(a.weak):                          # weak samples: extra ramets at about 2x
        g = clonable[w % len(clonable)]
        add_sample(f"{g}w{w + 1}", g, "weak")
if a.big_clone:
    g = "Bbig"
    hap = founder("B")
    if a.fixed_het:
        fx = rng.random(S) < a.fixed_het
        hap[fx, 0], hap[fx, 1] = 0, 1
    add_genet(g, hap, "B")
    for r in range(a.big_clone):
        add_sample(f"Bbig{r + 1:02d}", g, "bigclone")

n = len(samples)
names = [s[0] for s in samples]
# genotypes per sample: the genet's, plus somatic mutations for every ramet but the first
G = np.zeros((S, n), np.int8)
seen = set()
for i, (nmx, g, role) in enumerate(samples):
    hap = genets[g][0].copy()
    if g in seen and role != "techrep" and a.somatic:
        mut = rng.random(S) < a.somatic
        hap[mut, 0] = 1 - hap[mut, 0]
    seen.add(g)
    G[:, i] = hap.sum(1)

# coverage: lognormal around --cov; weak samples about 2x
cov = np.exp(rng.normal(np.log(a.cov), a.cov_sd, n))
for i, s in enumerate(samples):
    if s[2] == "weak":
        cov[i] = 2.0
depth = rng.poisson(cov[None, :], (S, n))
pa = G / 2.0                                     # chance a read carries the alternative allele
pa = pa * (1 - a.err) + (1 - pa) * a.err
alt = rng.binomial(depth, pa)
ref = depth - alt

# genotype likelihoods and posteriors with a frequency prior from the pooled sample (as -doPost 1)
e = a.err
lik = np.stack([((1 - e) ** ref) * (e ** alt),
                (0.5 ** depth) * np.ones_like(depth, float),
                (e ** ref) * ((1 - e) ** alt)], 2)
q = np.clip((alt.sum(1) + 0.5) / (depth.sum(1) + 1.0), 0.01, 0.99)
prior = np.stack([(1 - q) ** 2, 2 * q * (1 - q), q ** 2], 1)[:, None, :]
post = lik * prior
post /= post.sum(2, keepdims=True)

# IBS: one random read per sample and site
has = depth > 0
draw = rng.random((S, n)) < np.where(has, alt / np.maximum(depth, 1), 0)
A1 = (draw & has).astype(np.float32); A0 = ((~draw) & has).astype(np.float32); H = has.astype(np.float32)
both = H.T @ H
same = A1.T @ A1 + A0.T @ A0
ibs = 1 - same / np.maximum(both, 1)
np.fill_diagonal(ibs, 0)

# outputs
bams = os.path.join(a.outdir, f"{a.sp}Bams")
blist = "".join(f"../bams/{x}.bam\n" for x in names)
if os.path.exists(bams) and open(bams).read() != blist:
    sys.exit(f"{prog}: {a.outdir} already holds a different design (other samples); use a new --outdir for this one. "
             f"Variants in one folder must share the design (counts, --big-clone, --weak, --seed) and change only "
             f"read-level settings such as --cov, --err or --somatic")
d = os.path.join(a.outdir, "cloneTests", a.variant)
os.makedirs(d, exist_ok=True)
open(bams, "w").write(blist)
np.savetxt(os.path.join(d, f"{a.sp}.ibsMat"), ibs, fmt="%.6f", delimiter="\t")
with gzip.open(os.path.join(d, f"{a.sp}.geno.gz"), "wt") as fh:
    for s_ in range(S):
        fh.write(f"chr1\t{s_ + 1}\t" + "\t".join(f"{x:.4f}" for x in post[s_].ravel()) + "\n")
open(os.path.join(d, "command.txt"), "w").write(f"simulated: angsd -b {os.path.abspath(bams)}\n")
open(os.path.join(d, "settings.txt"), "w").write(" ".join(sys.argv) + "\n")
techreps = " ".join(s[0] for s in samples if s[2] == "techrep")
open(os.path.join(a.outdir, "techreps.txt"), "w").write(techreps + "\n")


def relation(gi, gj):
    if gi == gj:
        return "clone"
    mi, fi = genets[gi][2], genets[gi][3]
    mj, fj = genets[gj][2], genets[gj][3]
    if gi in (mj, fj) or gj in (mi, fi):
        return "parent-offspring"
    if mi is not None and mj is not None:
        if mi == fi and mj == fj and mi == mj:
            return "selfed-sibs"
        if {mi, fi} == {mj, fj}:
            return "full-sibs"
        if {mi, fi} & {mj, fj}:
            return "half-sibs"
    return "unrelated"


with open(os.path.join(a.outdir, "truth.tsv"), "w") as fh:
    fh.write("sample1\tsample2\trelation\n")
    for i, j in combinations(range(n), 2):
        r = relation(samples[i][1], samples[j][1])
        if r != "unrelated":
            fh.write(f"{names[i]}\t{names[j]}\t{r}\n")
ng = len({s[1] for s in samples})
from collections import Counter
cnt = Counter(s[1] for s in samples)
print(f"{n} samples, {S} sites, median coverage {np.median(cov):.1f}x; {ng} genets: "
      f"{sum(1 for v in cnt.values() if v > 1)} with copies, {n - ng} copies to drop (true answer)")
print(f"wrote {d}/{a.sp}.geno.gz, .ibsMat; {bams}; {a.outdir}/truth.tsv; technical replicates in {a.outdir}/techreps.txt")
