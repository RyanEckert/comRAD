inskip && skip "in SKIP"
cd "$WD/angsd" && mkdir -p ../het
# each sample on the SFS sites where it has at least the [het] -setMinDepthInd reads (10 by default, as in
# Matz's 2bRAD scripts). With HET_DEPTH set, every sample's reads are first thinned at random (fixed seed)
# to that mean depth, so all samples are scored at the same depth; samples below it are left out
rm -f ${SP}HetBams ${SP}HetExcluded ../het/${SP}HetFraction
if [ -n "$HET_DEPTH" ]; then
  # mean depth per sample over the SFS sites (filters as [het], no depth floor)
  CARGS=$(echo " $CFG_het_args " | sed -E 's/ -set(Min|Max)Depth(Ind)? +[0-9]+/ /g')
  dkey=$( { cat ${SP}Bams; md5sum <${SP}SFSsitesToDo; echo "$CARGS"; } | md5sum | cut -c1-32)
  [ "$(cat ../het/${SP}HetCounts.key 2>/dev/null)" = "$dkey" ] || { rm -f ../het/${SP}HetCounts.*; echo "$dkey" >../het/${SP}HetCounts.key; }
  job "test -s ../het/${SP}HetCounts.counts.gz" "angsd -b ${SP}Bams -sites ${SP}SFSsitesToDo -rf ${SP}SFSregions -ref $GENOME_FASTA -doCounts 1 -dumpCounts 2 -P $CFG_het_cpus $CARGS -out ../het/${SP}HetCounts" >${SP}HetCountJobs
  comrad-jobs ${SP}HetCountJobs -n ${SP}HetCount --res het || stop "per-sample depth (angsd -dumpCounts) failed"
  zcat ../het/${SP}HetCounts.counts.gz | awk 'NR > 1 {for (i = 1; i <= NF; i++) s[i] += $i; n++} END {for (i = 1; i in s; i++) printf "%.4f\n", s[i] / n}' |
    paste ${SP}Bams - >../het/${SP}HetMeanDepth.txt
  equal $(wc -l <../het/${SP}HetMeanDepth.txt) $(wc -l <${SP}Bams) "mean depths vs bams"
  awk -F'\t' -v t=$HET_DEPTH '$2 >= t {print $1 > "'${SP}HetBams'"; printf "%s\t%.6f\n", $1, t / $2 > "'../het/${SP}HetFraction'"; next}
    {print $1 "\t" $2 > "'${SP}HetExcluded'"}' ../het/${SP}HetMeanDepth.txt
  [ -s ${SP}HetExcluded ] || : >${SP}HetExcluded
  nb=$(wc -l <${SP}HetBams 2>/dev/null || echo 0)
  [ $nb -ge 2 ] || stop "fewer than 2 samples reach HET_DEPTH = ${HET_DEPTH}x mean depth (het/${SP}HetMeanDepth.txt)"
  note "het: reads thinned to ${HET_DEPTH}x mean depth over the SFS sites; $nb of $(wc -l <${SP}Bams) samples reach it$([ -s ${SP}HetExcluded ] && echo "; $(wc -l <${SP}HetExcluded) below it left out (angsd/${SP}HetExcluded)")"
else
  cp ${SP}Bams ${SP}HetBams; : >${SP}HetExcluded
  note "het: each sample on the SFS sites where it has at least the [het] -setMinDepthInd reads (no thinning; HET_DEPTH sets a common depth)"
fi
# a new sample set, site list, filters or depth target redoes everything below
key=$( { cat ${SP}HetBams ../het/${SP}HetFraction 2>/dev/null; md5sum <${SP}SFSsitesToDo; echo "$CFG_het_args $HET_DEPTH"; } | md5sum | cut -c1-32)
if [ "$(cat ../het/${SP}HetKey 2>/dev/null)" != "$key" ]; then
  rm -f ../het/*.ml ../het/*.saf.* ../het/*.arg ../het/${SP}HetSites*
  echo "$key" >../het/${SP}HetKey
fi
>${SP}HetJobs
while read b; do s=$(basename "$b" .bam)
  ds=$(awk -F'\t' -v b="$b" '$1 == b && $2 < 1 {print "-downSample " $2 " -seed 1"}' ../het/${SP}HetFraction 2>/dev/null)
  job "test -s ../het/$s.ml" "angsd -i $b -sites ${SP}SFSsitesToDo -rf ${SP}SFSregions -GL 1 -P $CFG_het_cpus -anc $GENOME_FASTA -ref $GENOME_FASTA -doSaf 1 -doCounts 1 $CFG_het_args $ds -out ../het/$s && localSaf.py realSFS ../het/$s.saf.idx -fold 1 >../het/$s.ml" >>${SP}HetJobs
done <${SP}HetBams
comrad-jobs ${SP}HetJobs -n ${SP}Het --res het || stop "heterozygosity failed"
cd ../het
{ printf 'sample\tlineage\tnSites\thet\n'
awk 'FILENAME ~ /BamsClusters$/ {n=$1; sub(/.*\//,"",n); sub(/\.bam$/,"",n); lin[n]=$2; next}
  {s=FILENAME; sub(/.*\//,"",s); sub(/\.ml$/,"",s); t=$1+$2+$3; printf "%s\t%s\t%d\t%.6f\n", s, (s in lin ? lin[s] : "NA"), t, $2/t}' ../angsd/${SP}BamsClusters $(sed 's#.*/##; s#\.bam$#.ml#' ../angsd/${SP}HetBams)
} >${SP}.het
equal $(( $(wc -l <${SP}.het) - 1 )) $(wc -l <../angsd/${SP}HetBams) "${SP}.het rows vs het samples"
awk -F'\t' 'NR > 1 && $3 == 0 {print $1; bad++} END {exit bad > 0}' ${SP}.het || stop "samples above have no site at the depth floor"
note "het: $(awk -F'\t' 'NR > 1 {print $3}' ${SP}.het | sort -n | awk '{a[NR] = $1} END {printf "%d to %d sites per sample (median %d)", a[1], a[NR], a[int((NR + 1) / 2)]}')"
na=$(awk -F'\t' 'NR>1 && $2 == "NA"' ${SP}.het | wc -l); [ $na -eq 0 ] || stop "$na samples in ${SP}.het have no lineage"
{ printf 'lineage\tnInd\tHo\n'
awk -F'\t' 'NR>1 {s[$2]+=$4; c[$2]++} END {for (l in s) printf "%s\t%d\t%.6f\n", l, c[l], s[l]/c[l]}' ${SP}.het | sort
} >${SP}Ho.txt
cat ${SP}Ho.txt
# does heterozygosity follow coverage? if it does, set HET_DEPTH (or raise it)
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
lines, warn = ["against\trho\tp\tr\tR2\tn"], []
for name, x in (("nSites", n), (f"propSites{os.environ['QDEPTH']}X", cov)):
    s = [k for k in h if k in x]
    if len(s) < 4:
        lines.append(f"{name}\tNA\tNA\tNA\tNA\t{len(s)}")
        continue
    rho, p = spearmanr([h[k] for k in s], [x[k] for k in s])
    if rho is None:
        lines.append(f"{name}\tNA\tNA\tNA\tNA\t{len(s)}")
        print(f"NOTE: heterozygosity vs {name}: no variation to test")
        continue
    r = pearson([h[k] for k in s], [x[k] for k in s])        # linear fit of het on coverage: R2 = r^2
    lines.append(f"{name}\t{rho:.3f}\t{p:.3g}\t{r:.3f}\t{r * r:.3f}\t{len(s)}")
    (warn if p < 0.05 else []).append(f"{name} rho = {rho:.2f}, p = {p:.2g}")
    print(f"NOTE: heterozygosity vs {name}: Spearman rho = {rho:.3f}, p = {p:.3g}; linear R2 = {r * r:.3f}; n = {len(s)}")
open(out, "w").write("\n".join(lines) + "\n")
if warn:
    print("NOTE: WARNING: heterozygosity correlates with coverage (" + "; ".join(warn) +
          "): set HET_DEPTH (or raise it) so every sample is scored at one depth (mean depths: het/" + het.replace(".het", "HetMeanDepth.txt") + " with HET_DEPTH set)")
PY

plots het
