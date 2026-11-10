inskip && skip "in SKIP"
inskip fst && skip "fst is in SKIP"
inskip outliers && skip "outliers is in SKIP: no outlier-free panel to compare with"
[ $(echo $POPS | wc -w) -ge 2 ] || skip "one lineage (K = 1): no FST between lineages; FST among sampling groups is in the groups step"
cd "$WD/angsd"
if [ ! -s ${SP}Outliers ]; then
cp ../sfs/${SP}LineageFst.out ../sfs/${SP}LineageFstAll.out
skip "no outliers: the full panel is the neutral panel; LineageFst.out copied"
fi
cut -f 1 ${SP}SFSallSites | uniq >${SP}SFSallRegions
angsd sites index ${SP}SFSallSites
FILTERS="$CFG_saf_args"
TODO="-doCounts 1 -doSaf 1 -anc $GENOME_FASTA -ref $GENOME_FASTA"
>${SP}SFSall
for POP in $POPS; do
job "gzip -t ../sfs/${SP}${POP}SFSall.saf.gz && test -s ../sfs/${SP}${POP}SFSall.saf.idx" "angsd -b ${SP}${POP}Bams -sites ${SP}SFSallSites -rf ${SP}SFSallRegions -GL 1 -P $CFG_saf_cpus $FILTERS $TODO -out ../sfs/${SP}${POP}SFSall" >>${SP}SFSall
done
comrad-jobs ${SP}SFSall -n ${SP}SFSall --res saf || stop "full-panel SAFs failed"
cd ../sfs
>${SP}2dSFSall
while read a b; do
job "test -s ${SP}${a}.${b}.all.sfs -a -s ${SP}${a}.${b}.all.fst.idx" "realSFS ${SP}${a}SFSall.saf.idx ${SP}${b}SFSall.saf.idx -maxIter 1000 -P $CFG_fst_cpus >${SP}${a}.${b}.all.sfs && realSFS fst index ${SP}${a}SFSall.saf.idx ${SP}${b}SFSall.saf.idx -sfs ${SP}${a}.${b}.all.sfs -whichFst 1 -fstout ${SP}${a}.${b}.all" >>${SP}2dSFSall
done < <(pairs)
comrad-jobs ${SP}2dSFSall -n ${SP}2dSFSall --res fst || stop "full-panel 2D SFS failed"
converged $(rtlogs ${SP}2dSFSall)
{ printf 'pop1\tpop2\tunweightedFst\tweightedFst\n'
while read a b; do v=$(fststats ${SP}${a}.${b}.all.fst.idx) || stop "no full-panel FST for $a-$b; see above"; printf '%s\t%s\t%s\n' "$a" "$b" "$v"; done < <(pairs)
} >${SP}LineageFstAll.out
# neutral vs full panel, weighted FST
paste <(cut -f 1,2,4 ${SP}LineageFst.out) <(cut -f 4 ${SP}LineageFstAll.out) | sed '1s/weightedFst\tweightedFst/neutralFst\tallFst/'

plots fst
