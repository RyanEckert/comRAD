inskip && skip "in SKIP"
cd "$WD/angsd"
K=$(getK)
pruneSnps.py -m ${SP}.mafs.gz -b ${SP}.beagle.gz -x ${SP}Outliers -o ${SP}Neutral || stop "pruneSnps.py failed"
# kin left out, so this is also the kin-free structure run
awk 'FILENAME==ARGV[1] {d[$1]; next} !($1 in d)' ../relatedness/${SP}KinDrop ${SP}Bams >${SP}NeutralBams
B=${SP}Neutral.beagle.gz
if [ -s ../relatedness/${SP}KinDrop ]; then
subsetBeagle.py -b ${SP}Neutral.beagle.gz -s ${SP}Bams -k ${SP}NeutralBams -o ${SP}NeutralNoKin.beagle.gz || stop "subsetBeagle.py failed"
B=${SP}NeutralNoKin.beagle.gz
fi
note "neutral run: $(wc -l <${SP}Neutral.sites) SNPs, $(wc -l <${SP}NeutralBams) individuals"
mkdir -p ../ngsAdmix
job "ls ../pcangsd/${SP}NeutralPcangsd.admix.*.Q >/dev/null 2>&1" "pcangsd -b $B --admix --threads $CFG_neutral_cpus -o ../pcangsd/${SP}NeutralPcangsd" >${SP}NeutralPcaJobs
comrad-jobs ${SP}NeutralPcaJobs -n ${SP}NeutralPca --res neutral || stop "neutral PCAngsd failed"

# K again, without kin or outlier loci: NGSadmix K = 1..MAXK, evalAdmix and the same reads on K as the K step
cd ../ngsAdmix
ngsAdmixLauncher.py -f ../angsd/$B --maxK $MAXK -r $ADMIX_REPS -n ${SP}Neutral -P $CFG_admix_cpus || stop "ngsAdmixLauncher.py failed"
comrad-jobs ${SP}NeutralNgsAdmix -n ${SP}NeutralNgsAdmix --res admix -c 'test -s {1}.qopt' -x='-o (\S+)' || stop "neutral NGSadmix failed"
pickQ.py --min 1 --max $MAXK --glob ${SP}Neutral --eval ../angsd/$B --threads $CFG_evaladmix_cpus || stop "pickQ.py failed"
comrad-jobs ${SP}Neutral_evalAdmixJobs -n ${SP}NeutralEvalAdmix --res evaladmix -c 'test -s {1}' -x='-o (\S+)' || stop "neutral evalAdmix failed"
pickQ.py --glob ${SP}Neutral --fit || stop "pickQ.py --fit failed"
admixResiduals.py -s ${SP}Neutral -k $MAXK -b ../angsd/${SP}NeutralBams --cut $RESID_MAX --meta "$META" --meta-sample "$META_SAMPLE" \
  --meta-groups $META_GROUPS || stop "neutral residual summary failed; see above"
kmeansBic.py -c ../pcangsd/${SP}NeutralPcangsd.cov -k $MAXK -o ${SP}Neutral_bic.txt || note "WARNING: neutral k-means BIC failed; see above"
kevidence ${SP}Neutral ${SP}NeutralPcangsd
kf=kin-free; [ -s ../relatedness/${SP}KinDrop ] || kf="kin kept"
note "$kf, neutral K: evalAdmix elbow ${ev:-none}, PCAngsd ${pc:-?}, residual-block ${rk:-none below $RESID_MAX}, k-means BIC ${bic:-?}; runs in the dominant mode per K: $sup (chosen K = $K; plots: ngsAdmix/${SP}NeutralAdmix.pdf, ngsAdmix/${SP}NeutralResiduals.pdf)"
[ "$ev" = "$K" ] || [ "$pc" = "$K" ] || [ "$rk" = "$K" ] || [ "$bic" = "$K" ] ||
  note "WARNING: $kf, without outlier loci, no read on K gives the chosen K = $K (elbow ${ev:-none}, PCAngsd ${pc:-?}, residual-block ${rk:-none}, k-means BIC ${bic:-?}); look at ngsAdmix/${SP}NeutralAdmix.pdf before going on"
echo; cat ${SP}Neutral_bestQ.txt; echo; cat ${SP}Neutral_fit.txt; echo; cat ${SP}Neutral_resid.txt; echo; cat ${SP}Neutral_bic.txt 2>/dev/null
cp ${SP}Neutral_k${K}_best.qopt ${SP}Neutral_best.qopt

# lineage calls from the neutral, kin-free run against the lineage assignment
python3 - ../angsd/${SP}BamsClusters ../angsd/${SP}NeutralBams ${SP}Neutral_best.qopt $ASSIGN ${SP}NeutralCompare.txt <<'PY' || stop "neutral comparison failed; see above"
import sys
from collections import Counter
cl, nb, qf, t, out = sys.argv[1:6]
t = float(t)
old = dict(l.split()[:2] for l in open(cl) if l.strip())
bams = [l.strip() for l in open(nb) if l.strip()]
Q = [[float(x) for x in l.split()] for l in open(qf) if l.strip()]
if len(Q) != len(bams):
    sys.exit(f"{qf} has {len(Q)} rows, {nb} {len(bams)}")
new = {}
for x, q in zip(bams, Q):
    m = max(q)
    new[x] = "Admixed" if m < t else f"neutral{q.index(m) + 1}"
tab = ["full\tneutral\tn"] + [f"{o}\t{n}\t{c}" for (o, n), c in sorted(Counter((old[x], new[x]) for x in bams).items())]
open(out, "w").write("\n".join(tab) + "\n")
print("\n".join(tab))
bad, hit = [], {}
for o in sorted({old[x] for x in bams} - {"Admixed"}):
    ns = {new[x] for x in bams if old[x] == o and new[x] != "Admixed"}
    if len(ns) != 1:
        bad.append(f"{o} -> {sorted(ns) or 'nothing pure'}")
        continue
    n = ns.pop()
    if n in hit:
        bad.append(f"{o} and {hit[n]} -> {n}")
    hit[n] = o
moved = [x for x in bams if (old[x] == "Admixed") != (new[x] == "Admixed")]
for x in moved:
    print(f"moved: {x}\t{old[x]} -> {new[x]}")
if bad:
    print("NOTE: WARNING: the neutral run does not match the lineages one to one (" + "; ".join(bad) +
          f"); {len(moved)} individuals moved to or from Admixed. Later steps use the lineages from the full SNP set; "
          "structure may be driven by outlier loci or kin (table: ngsAdmix/" + out + ")")
else:
    print(f"NOTE: neutral run matches the lineages; {len(moved)} individuals moved to or from Admixed")
PY

plots neutral
plots admix --prefix Neutral
