if inskip; then : >"$WD/angsd/${SP}Outliers"; skip "in SKIP: no outlier SNPs removed (empty angsd/${SP}Outliers)"; fi
cd "$WD/angsd/lineages"
for POP in $POPS; do
grep -E 'Loaded|MAF filtering|converge' ${SP}${POP}SelPcangsd.log
if grep -qi 'not converge' ${SP}${POP}SelPcangsd.log; then stop "$POP selection scan did not converge"; fi
pcangsdOutliers.py -b ${SP}${POP}.beagle.gz -s ${SP}${POP}SelPcangsd -o ${SP}${POP} || stop "pcangsdOutliers.py failed for $POP"
done
cat $(for POP in $POPS; do echo ${SP}${POP}Outliers; done) | sort -k1,1 -k2,2n -u >../${SP}Outliers
wc -l ${SP}*Outliers ../${SP}Outliers
note "$(wc -l <../${SP}Outliers) outlier SNPs (union over lineages)"
# one SNP per tag within each lineage, outlier tags removed
>${SP}LineagePcaJobs
for POP in $POPS; do
pruneSnps.py -m ${SP}${POP}.mafs.gz -b ${SP}${POP}.beagle.gz -x ../${SP}Outliers -o ${SP}${POP}Unlinked || stop "pruneSnps.py failed for $POP"
job "test -s ../../pcangsd/${SP}${POP}Pcangsd.cov" "pcangsd -b ${SP}${POP}Unlinked.beagle.gz --threads $CFG_outliers_cpus -o ../../pcangsd/${SP}${POP}Pcangsd" >>${SP}LineagePcaJobs
done
comrad-jobs ${SP}LineagePcaJobs -n ${SP}LineagePca --res outliers || stop "within-lineage PCA failed"

plots outliers
