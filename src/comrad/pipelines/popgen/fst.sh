inskip && skip "in SKIP"
[ $(echo $POPS | wc -w) -ge 2 ] || skip "one lineage (K = 1): no FST between lineages; FST among sampling groups is in the groups step"
cd "$WD/sfs"
>${SP}2dSFS
while read a b; do
job "test -s ${SP}${a}.${b}.sfs -a -s ${SP}${a}.${b}.fst.idx" "localSaf.py realSFS ${SP}${a}SFS.saf.idx ${SP}${b}SFS.saf.idx -maxIter 1000 -P $CFG_fst_cpus >${SP}${a}.${b}.sfs && localSaf.py realSFS fst index ${SP}${a}SFS.saf.idx ${SP}${b}SFS.saf.idx -sfs ${SP}${a}.${b}.sfs -whichFst 1 -fstout ${SP}${a}.${b}" >>${SP}2dSFS
done < <(pairs)
comrad-jobs ${SP}2dSFS -n ${SP}2dSFS --res fst || stop "2D SFS failed"
converged $(rtlogs ${SP}2dSFS)
{ printf 'pop1\tpop2\tunweightedFst\tweightedFst\n'
while read a b; do v=$(fststats ${SP}${a}.${b}.fst.idx) || stop "no FST for $a-$b; see above"; printf '%s\t%s\t%s\n' "$a" "$b" "$v"; done < <(pairs)
} >${SP}LineageFst.out
cat ${SP}LineageFst.out
while IFS=$'\t' read a b u w; do note "FST $a-$b weighted $w"; done < <(tail -n +2 ${SP}LineageFst.out)

plots fst
