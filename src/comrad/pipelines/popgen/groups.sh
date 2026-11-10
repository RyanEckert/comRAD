inskip && skip "in SKIP"
[ -n "$META_GROUPS" ] || skip "META_GROUPS is empty"
cd "$WD/angsd" && mkdir -p groups
python3 - <<'PY' || stop "groups failed; see above"
import csv, os, re, sys
from collections import Counter
e = os.environ
sp, minN = e["SP"], int(e["MIN_GROUP"])
grp = e["META_GROUPS"].split()
kin = {l.strip() for l in open(os.path.join("..", "relatedness", sp + "KinDrop")) if l.strip()}
with open(e["META"], newline="", encoding="utf-8-sig") as fh:
    meta = {r[e["META_SAMPLE"]]: r for r in csv.DictReader(fh)}
g = []
for l in open(sp + "BamsClusters"):
    if not l.strip():
        continue
    bam, lin = l.split()[:2]
    if lin == "Admixed" or bam in kin:
        continue
    s = re.sub(r"\.bam$", "", os.path.basename(bam))
    vals = [(meta.get(s) or {}).get(c, "").strip() for c in grp]
    g.append(dict(bam=bam, lineage=lin, sample=s, vals=vals))
na = [x["sample"] for x in g if not all(x["vals"])]
if na:
    sys.exit(f"{len(na)} samples without {'/'.join(grp)} in META, e.g. {na[0]}")
for x in g:
    x["group"] = re.sub(r"[^A-Za-z0-9_]", "", "_".join(x["vals"]))
n = Counter((x["lineage"], x["group"]) for x in g)
small = sorted((k, v) for k, v in n.items() if v < minN)
if small:
    print(f"NOTE: {len(small)} lineage x group combinations below MIN_GROUP ({minN}) left out: " +
          ", ".join(f"{l} {gr} (n = {v})" for (l, gr), v in small))
g = [x for x in g if n[(x["lineage"], x["group"])] >= minN]
print("lineage\tgroup\tn")
for (l, gr), v in sorted(Counter((x["lineage"], x["group"]) for x in g).items()):
    print(f"{l}\t{gr}\t{v}")
with open(os.path.join("groups", sp + "GroupBams"), "w") as fh:
    fh.write("".join(f"{x['bam']}\t{x['lineage']}\t{x['group']}\n" for x in g))
with open(sp + "Groups", "w") as fh:
    fh.write("".join(f"{x['sample']}\t{x['lineage']}\t{x['group']}\n" for x in g))
PY
cd groups
cut -f2,3 ${SP}GroupBams | sort -u >${SP}GroupList
while read L G; do awk -F'\t' -v l=$L -v g=$G '$2==l && $3==g {print $1}' ${SP}GroupBams >${SP}${L}_${G}Bams; done <${SP}GroupList
awk '{g[$1] = g[$1] " " $2} END {for (l in g) {n = split(g[l], a, " "); for (i = 1; i < n; i++) for (j = i + 1; j <= n; j++) print l, a[i], a[j]}}' ${SP}GroupList >${SP}GroupPairs
[ -s ${SP}GroupPairs ] || skip "no lineage has two groups of at least $MIN_GROUP"
cd ..
FILTERS="$CFG_saf_args"
TODO="-doCounts 1 -doSaf 1 -anc $GENOME_FASTA -ref $GENOME_FASTA"
>${SP}GroupSFS
while read L G; do
job "test -s groups/${SP}${L}_${G}.saf.idx" "angsd -b groups/${SP}${L}_${G}Bams -sites ${SP}SFSsitesToDo -rf ${SP}SFSregions -GL 1 -P $CFG_saf_cpus $FILTERS $TODO -out groups/${SP}${L}_${G}" >>${SP}GroupSFS
done <groups/${SP}GroupList
comrad-jobs ${SP}GroupSFS -n ${SP}GroupSFS --res saf || stop "group SAFs failed"
cd groups
>${SP}GroupFstJobs
while read L a b; do
job "test -s ${SP}${L}_${a}.${b}.fst.idx" "realSFS ${SP}${L}_${a}.saf.idx ${SP}${L}_${b}.saf.idx -maxIter 1000 -P $CFG_fst_cpus >${SP}${L}_${a}.${b}.sfs && realSFS fst index ${SP}${L}_${a}.saf.idx ${SP}${L}_${b}.saf.idx -sfs ${SP}${L}_${a}.${b}.sfs -whichFst 1 -fstout ${SP}${L}_${a}.${b}" >>${SP}GroupFstJobs
done <${SP}GroupPairs
comrad-jobs ${SP}GroupFstJobs -n ${SP}GroupFst --res fst || stop "group FST failed"
converged $(rtlogs ${SP}GroupFst)
{ printf 'lineage\tgroup1\tgroup2\tunweightedFst\tweightedFst\n'
while read L a b; do v=$(fststats ${SP}${L}_${a}.${b}.fst.idx) || stop "no FST for $L $a-$b; see above"; printf '%s\t%s\t%s\t%s\n' "$L" "$a" "$b" "$v"; done <${SP}GroupPairs
} >${SP}GroupFst.out
cat ${SP}GroupFst.out

plots fst
