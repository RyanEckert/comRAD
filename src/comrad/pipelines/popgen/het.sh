inskip && skip "in SKIP"
cd "$WD/angsd" && mkdir -p ../het
>${SP}HetJobs
while read b; do s=$(basename "$b" .bam)
job "test -s ../het/$s.ml" "angsd -i $b -sites ${SP}SFSsitesToDo -rf ${SP}SFSregions -GL 1 -P $CFG_het_cpus -anc $GENOME_FASTA -ref $GENOME_FASTA -doSaf 1 -doCounts 1 $CFG_het_args -out ../het/$s && realSFS ../het/$s.saf.idx -fold 1 >../het/$s.ml" >>${SP}HetJobs
done <${SP}Bams
comrad-jobs ${SP}HetJobs -n ${SP}Het --res het || stop "heterozygosity failed"
cd ../het
{ printf 'sample\tlineage\tnSites\thet\n'
awk 'FILENAME ~ /BamsClusters$/ {n=$1; sub(/.*\//,"",n); sub(/\.bam$/,"",n); lin[n]=$2; next}
  {s=FILENAME; sub(/.*\//,"",s); sub(/\.ml$/,"",s); t=$1+$2+$3; printf "%s\t%s\t%d\t%.6f\n", s, (s in lin ? lin[s] : "NA"), t, $2/t}' ../angsd/${SP}BamsClusters *.ml
} >${SP}.het
equal $(( $(wc -l <${SP}.het) - 1 )) $(wc -l <../angsd/${SP}Bams) "${SP}.het rows vs bams"
na=$(awk -F'\t' 'NR>1 && $2 == "NA"' ${SP}.het | wc -l); [ $na -eq 0 ] || stop "$na samples in ${SP}.het have no lineage"
{ printf 'lineage\tnInd\tHo\n'
awk -F'\t' 'NR>1 {s[$2]+=$4; c[$2]++} END {for (l in s) printf "%s\t%d\t%.6f\n", l, c[l], s[l]/c[l]}' ${SP}.het | sort
} >${SP}Ho.txt
cat ${SP}Ho.txt
# does heterozygosity follow coverage? the 10x floor removes the per-site bias, but
# low-coverage samples reach 10x at fewer (and different) sites
QDEPTH=$QDEPTH python3 - ${SP}.het ../angsd/${SP}Quality${QDEPTH}X.txt ${SP}HetDepth.txt <<'PY' || stop "heterozygosity depth check failed; see above"
import os, random, sys
het, q5, out = sys.argv[1:4]


def ranks(v):
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def pearson(x, y):
    n = len(x); mx = sum(x) / n; my = sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sx = sum((a - mx) ** 2 for a in x) ** .5; sy = sum((b - my) ** 2 for b in y) ** .5
    return None if sx == 0 or sy == 0 else sxy / (sx * sy)


def spearmanr(x, y, perms=9999):
    """Spearman's rho, two-sided permutation p (no scipy needed)"""
    rx, ry = ranks(x), ranks(y)
    rho = pearson(rx, ry)
    if rho is None:
        return None, None
    rnd = random.Random(1)
    hits = 0
    for _ in range(perms):
        p = ry[:]
        rnd.shuffle(p)
        r = pearson(rx, p)
        if r is not None and abs(r) >= abs(rho) - 1e-12:
            hits += 1
    return rho, (hits + 1) / (perms + 1)
rows = [l.rstrip("\n").split("\t") for l in open(het)][1:]
cov = {}
for l in open(q5):
    f = l.split()
    if len(f) >= 2:
        cov[os.path.basename(f[0]).removesuffix(".bam")] = float(f[1])
h = {r[0]: float(r[3]) for r in rows}
n = {r[0]: float(r[2]) for r in rows}
lines, warn = ["against\trho\tp\tn"], []
for name, x in (("nSites", n), (f"propSites{os.environ['QDEPTH']}X", cov)):
    s = [k for k in h if k in x]
    if len(s) < 4:
        lines.append(f"{name}\tNA\tNA\t{len(s)}")
        continue
    rho, p = spearmanr([h[k] for k in s], [x[k] for k in s])
    if rho is None:
        lines.append(f"{name}\tNA\tNA\t{len(s)}")
        print(f"NOTE: heterozygosity vs {name}: no variation to test")
        continue
    lines.append(f"{name}\t{rho:.3f}\t{p:.3g}\t{len(s)}")
    (warn if p < 0.05 else []).append(f"{name} rho = {rho:.2f}, p = {p:.2g}")
    print(f"NOTE: heterozygosity vs {name}: Spearman rho = {rho:.3f}, p = {p:.3g}, n = {len(s)}")
open(out, "w").write("\n".join(lines) + "\n")
if warn:
    print("NOTE: WARNING: heterozygosity correlates with coverage (" + "; ".join(warn) +
          "): compare groups with coverage as a covariate, or downsample the bams to one depth")
PY

plots het
