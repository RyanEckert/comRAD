cd "$WD/angsd"
FILTERS="$CFG_genotype_args -minInd $(mind ${SP}Bams ${CFG_genotype_min_ind:-$MIN_IND})"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 15 -doBcf 1 -doPost 1 -doGlf 2"
# bams copied to node-local disk first: beegfs reads gave bgzf errors (claude/beegfs-bgzf-reads.md)
cat >${SP}Angsd.sh <<SH
cd $WD/angsd
rm -f ${SP}.FAILED
L=\$TMPDIR/bams; mkdir -p \$L
cp \$(sed 's/\$/*/' ${SP}Bams) \$L/
sed "s#.*/#\$L/#" ${SP}Bams >\$L/bams
angsd -b \$L/bams -sites ${SP}KeepSites -rf ${SP}KeepRegions -GL 1 $FILTERS $TODO -P $CFG_genotype_cpus -out ${SP} >${SP}.log 2>&1
n=\$(grep -c bgzf ${SP}.log)
if [ \$n -gt 0 ]; then echo "\$n bgzf errors" >${SP}.FAILED; exit 1; fi
SH
job "test -s ${SP}.mafs.gz -a ! -e ${SP}.FAILED" "bash ${SP}Angsd.sh" >${SP}Job
comrad-jobs ${SP}Job -n ${SP}Angsd --res genotype || stop "genotyping failed$( [ -e ${SP}.FAILED ] && echo ": $(cat ${SP}.FAILED)")"
grep "filtering:" ${SP}.log
n=$(zcat ${SP}.mafs.gz | tail -n +2 | wc -l); [ $n -gt 0 ] || stop "no SNPs in ${SP}.mafs.gz"
note "$n SNPs"

# depth at the SNPs, per kept sample
paste <(zcat ${SP}Depth.pos.gz) <(zcat ${SP}Depth.counts.gz) | awk -v OFS='\t' -v sp=$SP '
  FILENAME==ARGV[1] {s=$0; sub(/.*\//,"",s); sub(/\.bam$/,"",s); cl[FNR]=s; nc=FNR; next}
  FILENAME==ARGV[2] {s=$0; sub(/.*\//,"",s); sub(/\.bam$/,"",s); ord[++nk]=s; keep[s]=1; next}
  FILENAME==ARGV[3] {if (FNR>1) {snp[$1":"$2]=1; ns++}; next}
  FNR==1 {if (NF-3 != nc) {printf "counts has %d columns, AllBams %d lines\n", NF-3, nc >"/dev/stderr"; err=1; exit 1}; next}
  ($1":"$2) in snp {m++; for (i=1; i<=nc; i++) if (cl[i] in keep) {d[cl[i]]+=$(i+3); if ($(i+3)>0) c[cl[i]]++}}
  END {
    if (err) exit 1
    if (m != ns) {printf "%d of %d SNPs found in %sDepth.pos.gz\n", m, ns, sp >"/dev/stderr"; exit 1}
    for (i=1; i<=nc; i++) seen[cl[i]]=1
    for (j=1; j<=nk; j++) if (!(ord[j] in seen)) {print ord[j]" not in AllBams" >"/dev/stderr"; exit 1}
    print "bam", "nSNPs", "meanDepth", "propCovered"
    for (j=1; j<=nk; j++) printf "%s\t%d\t%.4f\t%.4f\n", ord[j], m, d[ord[j]]/m, c[ord[j]]/m
  }' ${SP}AllBams ${SP}Bams <(zcat ${SP}.mafs.gz) - >${SP}SnpDepth.txt || stop "depth at SNPs failed; see above"
# a sample with little data at the final SNPs can sit apart in the PCA for technical reasons
med=$(awk -F'\t' 'NR>1 {print $4}' ${SP}SnpDepth.txt | sort -g | awk '{v[NR]=$1} END {printf "%.4f", NR%2 ? v[(NR+1)/2] : (v[NR/2]+v[NR/2+1])/2}')
low=$(awk -F'\t' -v m=$med 'NR>1 && $4 < m/2 {printf "%s (%.2f) ", $1, $4}' ${SP}SnpDepth.txt)
note "median proportion of SNPs covered per sample: $med"
[ -z "$low" ] || note "WARNING: samples covering under half the median proportion of SNPs: $low(check them in the PCA)"

# clone-free check: the kept samples' tree with the clone cut (same distance); pairs below it are flagged, not removed
CM=$(ls ${SP}CloneMetric.txt clones/${SP}CloneMetric.txt 2>/dev/null | head -1)
if [ -n "$CM" ]; then
  met=$(awk '$1=="metric"{print $2}' $CM); conf=$(awk '$1=="confidence"{print $2}' $CM); cut=$(awk '$1=="cut"{print $2}' $CM)
  [ "$conf" != NA ] || conf=
  cloneGroups.py -i ${SP}.ibsMat -g ${SP}.geno.gz -b ${SP}Bams --metric $met ${conf:+--conf $conf} --check $cut \
    --quality ${SP}Quality${QDEPTH}X.txt --meta "$META" --meta-sample "$META_SAMPLE" --meta-groups $META_GROUPS ${PLOT_COLOR:+--color "$PLOT_COLOR"} ${PLOT_SHAPE:+--shape "$PLOT_SHAPE"} --levels "$PLOT_LEVELS" --colors "$PLOT_COLORS" --shapes "$PLOT_SHAPES" -o $SP \
    || note "WARNING: clone-free check failed; see above"
fi
