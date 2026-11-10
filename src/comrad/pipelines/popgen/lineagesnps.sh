cd "$WD/angsd" && mkdir -p lineages
FILTERS="$CFG_lineagesnps_args"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 15 -doPost 1 -doGlf 2 -doHWE 1"
>${SP}LineageJobs
for POP in $POPS; do
job "test -s lineages/${SP}${POP}.beagle.gz -a -s lineages/${SP}${POP}.mafs.gz" "angsd -b ${SP}${POP}Bams -sites ${SP}KeepSites -rf ${SP}KeepRegions -GL 1 -P $CFG_lineagesnps_cpus $FILTERS $TODO -minInd $(mind ${SP}${POP}Bams ${CFG_lineagesnps_min_ind:-$MIN_IND}) -out lineages/${SP}${POP}" >>${SP}LineageJobs
done
comrad-jobs ${SP}LineageJobs -n ${SP}Lineage --res lineagesnps || stop "within-lineage ANGSD failed"
cd lineages
# heterozygosity per site within each lineage: sites where most samples are called heterozygous are
# likely collapsed paralogs; reported, and removed only when LINEAGE_MAXHET is set
for POP in $POPS; do
python3 - ${SP}${POP} "$LINEAGE_MAXHET" <<'PY' || stop "lineage heterozygosity check failed for $POP; see above"
import gzip, os, sys
import numpy as np
pre, maxhet = sys.argv[1], sys.argv[2]
share = []
with gzip.open(pre + ".geno.gz", "rt") as fh:
    for l in fh:
        t = l.split()[4:]                         # -doGeno 15: chr pos major minor, then 5 per sample
        q = np.array(t, dtype=object).reshape(-1, 5)[:, 2:].astype(float)
        ok = q.max(1) > 0.95
        share.append((q.argmax(1)[ok] == 1).mean() if ok.sum() >= 5 else np.nan)
share = np.array(share)
m = ~np.isnan(share)
hi = m & (share > 0.5)
line = (f"{pre}: {len(share)} sites; median share of confident calls heterozygous {np.median(share[m]):.3f}; "
        f"{hi.sum()} sites ({hi.mean():.1%}) heterozygous in more than half the samples")
if not maxhet:
    print("NOTE:", line, "(kept; LINEAGE_MAXHET removes them)")
    sys.exit(0)
drop = m & (share > float(maxhet))
for ext, head in ((".beagle.gz", 1), (".mafs.gz", 1)):
    src = pre + ext
    if not os.path.exists(src + ".all") or os.path.getmtime(src) > os.path.getmtime(src + ".all"):
        os.replace(src, src + ".all")             # a fresh ANGSD output (the filtered copy keeps .all's time)
    with gzip.open(src + ".all", "rt") as fi, gzip.open(src, "wt") as fo:
        for k, l in enumerate(fi):
            if k < head or not drop[k - head]:
                fo.write(l)
    t = os.path.getmtime(src + ".all"); os.utime(src, (t, t))
print("NOTE:", line + f"; {drop.sum()} with more than {maxhet} removed from {pre}.beagle.gz and .mafs.gz (originals: .all)")
PY
done
>${SP}SelectionJobs
for POP in $POPS; do
job "test -s ${SP}${POP}SelPcangsd.selection" "pcangsd -b ${SP}${POP}.beagle.gz --selection --sites-save --iter 1000 --threads $CFG_selection_cpus -o ${SP}${POP}SelPcangsd" >>${SP}SelectionJobs
done
comrad-jobs ${SP}SelectionJobs -n ${SP}LineageSel --res selection || stop "selection scan failed"
