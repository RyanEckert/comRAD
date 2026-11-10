cd "$WD/angsd" && mkdir -p ../sfs
FILTERS="$CFG_saf_args"
TODO="-doCounts 1 -doSaf 1 -anc $GENOME_FASTA -ref $GENOME_FASTA"
>${SP}SFS
for POP in $POPS; do
job "gzip -t ../sfs/${SP}${POP}SFS.saf.gz && test -s ../sfs/${SP}${POP}SFS.saf.idx" "angsd -b ${SP}${POP}Bams -sites ${SP}SFSsitesToDo -rf ${SP}SFSregions -GL 1 -P $CFG_saf_cpus $FILTERS $TODO -out ../sfs/${SP}${POP}SFS" >>${SP}SFS
done
comrad-jobs ${SP}SFS -n ${SP}SFS --res saf || stop "SAFs failed"
