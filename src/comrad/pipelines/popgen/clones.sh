cd "$WD/angsd"
# near-empty samples look close to everything in IBS and distort the clone tree; leave
# them out of the clone step (they fall below MINQ later in any case)
rm -f ${SP}CloneLowCoverage.txt
awk -v m="$CLONE_MINQ" 'NR==FNR {q[$1]=$2; next}
  ($1 in q) && q[$1]+0 < m+0 {printf "%s\t%s\n", $1, q[$1] >"'${SP}CloneLowCoverage.txt'"; next} {print}' \
  ${SP}Quality${QDEPTH}X.txt ${SP}AllBams >${SP}CloneBams
touch ${SP}CloneLowCoverage.txt
nl=$(wc -l <${SP}CloneLowCoverage.txt)
if [ $nl -gt 0 ]; then
note "$nl of $(wc -l <${SP}AllBams) samples below CLONE_MINQ ($CLONE_MINQ of sites at >= ${QDEPTH}x) left out of the clone step: $(awk -F'\t' '{s=$1; sub(/.*\//,"",s); sub(/\.bam$/,"",s); printf "%s (%.3f) ", s, $2}' ${SP}CloneLowCoverage.txt)"
else
note "no samples below CLONE_MINQ ($CLONE_MINQ of sites at >= ${QDEPTH}x)"
fi
FILTERS="$CFG_clones_args -minInd $(mind ${SP}CloneBams ${CFG_clones_min_ind:-$MIN_IND})"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"
job "test -s ${SP}Clones.ibsMat" "angsd -b ${SP}CloneBams -GL 1 $FILTERS $TODO -P $CFG_clones_cpus -out ${SP}Clones" >${SP}CloneJobs
comrad-jobs ${SP}CloneJobs -n ${SP}Clones --res clones || stop "clone ANGSD run failed"
