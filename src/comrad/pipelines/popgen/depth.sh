mkdir -p "$WD/angsd" && cd "$WD/angsd"
# the bam lists (from the clean step), remade if angsd/ was cleared
if [ ! -s ${SP}AllBams ]; then
  ls ../bams/*.bam 2>/dev/null | grep -v '\.sym\.' >${SP}AllBams
  hasSym && ls ../bams/*.sym.bam >${SP}SymBams
  [ -s ${SP}AllBams ] || stop "no bams in bams/; run  comrad from popgen map"
  note "bam lists remade from bams/ ($(wc -l <${SP}AllBams) samples)"
fi
NB=$(wc -l <${SP}AllBams)
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth $(( NB * 50 )) -minInd $(mind ${SP}AllBams $CFG_depth_min_ind)"
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
job "test -s ${SP}Depth.depthSample" "angsd -b ${SP}AllBams -GL 1 $FILTERS $TODO -P $CFG_depth_cpus -out ${SP}Depth" >${SP}DepthJobs
comrad-jobs ${SP}DepthJobs -n ${SP}Depth --res depth || stop "depth pass failed"
QS=$(printf "%s\n" 3 5 10 $QDEPTH | sort -nu | tr "\n" " ")
job "test -s ${SP}Quality${QDEPTH}X.txt" "plotQC.py -p ${SP}Depth -b ${SP}AllBams -q $QS -c 3 5 10 && for q in $QS; do mv quality\${q}X.txt ${SP}Quality\${q}X.txt; done && if [ -e ${SP}Depth.counts ]; then gzip -9f ${SP}Depth.counts; fi" >${SP}QCJobs
comrad-jobs ${SP}QCJobs -n ${SP}QC --res depth || stop "plotQC failed"
head -5 ${SP}Quality${QDEPTH}X.txt
