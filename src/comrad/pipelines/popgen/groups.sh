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
job "test -s ${SP}${L}_${a}.${b}.fst.idx" "localSaf.py realSFS ${SP}${L}_${a}.saf.idx ${SP}${L}_${b}.saf.idx -maxIter 1000 -P $CFG_fst_cpus >${SP}${L}_${a}.${b}.sfs && localSaf.py realSFS fst index ${SP}${L}_${a}.saf.idx ${SP}${L}_${b}.saf.idx -sfs ${SP}${L}_${a}.${b}.sfs -whichFst 1 -fstout ${SP}${L}_${a}.${b}" >>${SP}GroupFstJobs
done <${SP}GroupPairs
comrad-jobs ${SP}GroupFstJobs -n ${SP}GroupFst --res fst || stop "group FST failed"
converged $(rtlogs ${SP}GroupFst)
{ printf 'lineage\tgroup1\tgroup2\tunweightedFst\tweightedFst\n'
while read L a b; do v=$(fststats ${SP}${L}_${a}.${b}.fst.idx) || stop "no FST for $L $a-$b; see above"; printf '%s\t%s\t%s\t%s\n' "$L" "$a" "$b" "$v"; done <${SP}GroupPairs
} >${SP}GroupFst.out
cat ${SP}GroupFst.out

# diversity per group: pi, Watterson's theta and Tajima's D, every group on the same sites (those in all the
# group SAFs), from all its samples and from the smallest group's n, since theta W and Tajima's D depend on n
NG=$(cut -f2,3 ${SP}GroupBams | sort | uniq -c | awk 'NR == 1 || $1 < m {m = $1} END {print m}')
mkdir -p sub
>${SP}GroupSubSAF
while read L G; do
B=sub/${SP}${L}_${G}SubBams      # bam paths as in GroupBams, relative to angsd/, so the job runs from there
# the same individuals on every run (ordered by an md5 of the bam name)
while read -r b; do printf '%s\t%s\n' "$(printf '%s' "$(basename "$b")" | md5sum | cut -c1-32)" "$b"; done <${SP}${L}_${G}Bams |
  sort -k1,1 | awk -F'\t' -v n=$NG 'NR <= n {print $2}' | sort >$B.new
if ! cmp -s $B.new $B; then rm -f sub/${SP}${L}_${G}Sub.* sub/${SP}${L}_${G}Sub_*; mv $B.new $B; else rm $B.new; fi
job "gzip -t sub/${SP}${L}_${G}Sub.saf.gz && test -s sub/${SP}${L}_${G}Sub.saf.idx" "cd .. && angsd -b groups/$B -sites ${SP}SFSsitesToDo -rf ${SP}SFSregions -GL 1 -P $CFG_saf_cpus $FILTERS $TODO -out groups/sub/${SP}${L}_${G}Sub" >>${SP}GroupSubSAF
done <${SP}GroupList
comrad-jobs ${SP}GroupSubSAF -n ${SP}GroupSubSAF --res saf || stop "subsampled group SAFs failed"
SAFS=$(while read L G; do echo ${SP}${L}_${G}.saf.idx sub/${SP}${L}_${G}Sub.saf.idx; done <${SP}GroupList)
commonsites ${SP}GroupThetaPos ${SP}GroupThetaSites $SAFS &&
  rm -f ${SP}*_theta.thetas.* ${SP}*.folded.sfs sub/*_theta.thetas.* sub/*.folded.sfs
NT=$(wc -l <${SP}GroupThetaSites)
[ $NT -gt 0 ] || stop "no site is in every group's SAF"
note "group diversity: $NT sites in every group's SAF; theta per group from all its samples and from $NG per group (the smallest group)"
ST="-sites ${SP}GroupThetaSites"
>${SP}GroupThetas
while read L G; do for S in ${SP}${L}_${G} sub/${SP}${L}_${G}Sub; do
job "test -s ${S}_theta.thetas.idx.pestPG" "localSaf.py realSFS $S.saf.idx $ST -fold 1 -maxIter 1000 -P $CFG_theta_cpus >$S.folded.sfs && localSaf.py realSFS saf2theta $S.saf.idx $ST -sfs $S.folded.sfs -outname ${S}_theta -fold 1 && localSaf.py thetaStat do_stat ${S}_theta.thetas.idx" >>${SP}GroupThetas
done; done <${SP}GroupList
comrad-jobs ${SP}GroupThetas -n ${SP}GroupThetas --res theta || stop "group thetas failed"
converged $(rtlogs ${SP}GroupThetas)
{ printf 'lineage\tgroup\tn\tpi\tthetaW\ttajimaD\tnSites\tnSub\tpiSub\tthetaWSub\ttajimaDSub\tnSitesSub\n'
while read L G; do
printf '%s\t%s\t%d\t%s\t%d\t%s\n' $L $G $(wc -l <${SP}${L}_${G}Bams) "$(tstat ${SP}${L}_${G}_theta.thetas.idx.pestPG $(wc -l <${SP}${L}_${G}Bams))" \
  $NG "$(tstat sub/${SP}${L}_${G}Sub_theta.thetas.idx.pestPG $NG)"
done <${SP}GroupList
} >${SP}GroupTheta.out
cat ${SP}GroupTheta.out

plots fst --update
