cd "$WD/angsd"
[ -e ${SP}CloneBams ] || mv clones/${SP}Clone* .     # a rerun after the move below
filterBams.py -b ${SP}CloneBams -q ${SP}Quality${QDEPTH}X.txt -p $MINQ -c ${SP}CloneGroups -o ${SP}Bams || stop "filterBams.py failed"
need ${SP}Bams
note "$(wc -l <${SP}Bams) of $(wc -l <${SP}AllBams) samples kept ($(wc -l <${SP}CloneLowCoverage.txt) below CLONE_MINQ before the clone step; the rest dropped as clones or below MINQ)"
mkdir -p clones
mv ${SP}Clone* clones/
