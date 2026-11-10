cd "$WD/angsd"
FILTERS="$CFG_sfssites_args"
TODO="-doCounts 1 -doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 11 -doHWE 1"
>${SP}Discover
for POP in $POPS; do
job "test -s ${SP}${POP}Disc.geno.gz" "angsd -b ${SP}${POP}Bams -sites ${SP}KeepSites -rf ${SP}KeepRegions -GL 1 -P $CFG_sfssites_cpus -ref $GENOME_FASTA $FILTERS $TODO -minInd $(mind ${SP}${POP}Bams ${CFG_sfssites_min_ind:-$MIN_IND}) -out ${SP}${POP}Disc" >>${SP}Discover
done
comrad-jobs ${SP}Discover -n ${SP}Discover --res sfssites || stop "site discovery failed"
for POP in $POPS; do zcat ${SP}${POP}Disc.geno.gz | awk '{print $1"\t"$2}' >${SP}${POP}DiscSites; done
awk '(++c[$0])==(ARGC-1)' $(for POP in $POPS; do echo ${SP}${POP}DiscSites; done) >${SP}SFSallSites
# drop sites within 35 bp of an outlier SNP (the outlier tags)
awk -v W=35 'FILENAME==ARGV[1] {o[$1]=o[$1]" "$2; next}
  {n=split(o[$1], v, " "); for (i=1; i<=n; i++) if ($2>=v[i]-W && $2<=v[i]+W) next; print}' \
  ${SP}Outliers ${SP}SFSallSites >${SP}SFSsitesToDo
need ${SP}SFSsitesToDo
cut -f 1 ${SP}SFSsitesToDo | uniq >${SP}SFSregions
note "SFS sites: $(wc -l <${SP}SFSsitesToDo) neutral of $(wc -l <${SP}SFSallSites)"
sitesindex ${SP}SFSsitesToDo
