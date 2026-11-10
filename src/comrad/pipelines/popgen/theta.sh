inskip && skip "in SKIP"
mkdir -p "$WD/theta"
cd "$WD/angsd"
# THETA_N: also estimate theta from the same number of individuals in every lineage, since
# Watterson's theta and Tajima's D depend on n (pi much less so)
NSUB=
if [ -n "$THETA_N" ]; then
min=$(for POP in $POPS; do wc -l <${SP}${POP}Bams; done | sort -n | head -1)
case $THETA_N in
min) NSUB=$min ;;
*[!0-9]*) stop "THETA_N = $THETA_N: leave it empty, or give a whole number or min" ;;
*) NSUB=$THETA_N ;;
esac
[ $NSUB -ge 2 ] || stop "THETA_N = $THETA_N gives $NSUB individuals; theta needs at least 2"
[ $NSUB -le $min ] || stop "THETA_N = $NSUB is more than the smallest lineage ($min individuals); set it to min or lower"
mkdir -p ../theta/sub
>${SP}ThetaSubSAF
for POP in $POPS; do
B=../theta/sub/${SP}${POP}SubBams
# the same individuals on every run and machine (ordered by an md5 of the bam name);
# a new sample size or lineage list redoes that lineage
while read -r b; do printf '%s\t%s\n' "$(printf '%s' "$(basename "$b")" | md5sum | cut -c1-32)" "$b"; done <${SP}${POP}Bams |
  sort -k1,1 | awk -F'\t' -v n=$NSUB 'NR <= n {print $2}' | sort >$B.new
[ $(wc -l <$B.new) -eq $NSUB ] || stop "subsampling $POP to $NSUB individuals gave $(wc -l <$B.new); see above"
if ! cmp -s $B.new $B; then rm -f ../theta/sub/${SP}${POP}Sub.* ../theta/sub/${SP}${POP}Sub_*; mv $B.new $B; else rm $B.new; fi
job "gzip -t ../theta/sub/${SP}${POP}Sub.saf.gz && test -s ../theta/sub/${SP}${POP}Sub.saf.idx" "angsd -b $B -sites ${SP}SFSsitesToDo -rf ${SP}SFSregions -GL 1 -P $CFG_saf_cpus $CFG_saf_args -doCounts 1 -doSaf 1 -anc $GENOME_FASTA -ref $GENOME_FASTA -out ../theta/sub/${SP}${POP}Sub" >>${SP}ThetaSubSAF
done
comrad-jobs ${SP}ThetaSubSAF -n ${SP}ThetaSubSAF --res saf || stop "subsampled SAFs failed"
note "theta also from $NSUB individuals per lineage (THETA_N = $THETA_N); individuals listed in theta/sub/${SP}<lineage>SubBams"
else
rm -rf ../theta/sub     # from an earlier run with THETA_N set
fi

cd "$WD/sfs"
# every lineage (and subsample) on the same sites: those in all of their SAFs. A SAF keeps only sites with
# reads in that lineage, so without this the lineages' site sets differ slightly
SAFS=$(for POP in $POPS; do echo ${SP}${POP}SFS.saf.idx; [ -z "$NSUB" ] || echo ../theta/sub/${SP}${POP}Sub.saf.idx; done)
commonsites ${SP}ThetaPos ../theta/${SP}ThetaSites $SAFS &&
  rm -f ../theta/${SP}*.thetas.* ../theta/sub/*_theta.thetas.* ${SP}*.folded.sfs ../theta/sub/*.folded.sfs
NT=$(wc -l <../theta/${SP}ThetaSites)
[ $NT -gt 0 ] || stop "no site is in every lineage's SAF"
note "theta: $NT sites in every lineage's SAF$([ -z "$NSUB" ] || echo " and every subsample's"); all theta from these"
ST="-sites ../theta/${SP}ThetaSites"
>${SP}Thetas
for POP in $POPS; do
job "test -s ../theta/${SP}${POP}.thetas.idx.pestPG" "localSaf.py realSFS ${SP}${POP}SFS.saf.idx $ST -fold 1 -maxIter 1000 -P $CFG_theta_cpus >${SP}${POP}.folded.sfs && localSaf.py realSFS saf2theta ${SP}${POP}SFS.saf.idx $ST -sfs ${SP}${POP}.folded.sfs -outname ../theta/${SP}${POP} -fold 1 && localSaf.py thetaStat do_stat ../theta/${SP}${POP}.thetas.idx" >>${SP}Thetas
[ -z "$NSUB" ] || { S=../theta/sub/${SP}${POP}Sub
job "test -s ${S}_theta.thetas.idx.pestPG" "localSaf.py realSFS $S.saf.idx $ST -fold 1 -maxIter 1000 -P $CFG_theta_cpus >$S.folded.sfs && localSaf.py realSFS saf2theta $S.saf.idx $ST -sfs $S.folded.sfs -outname ${S}_theta -fold 1 && localSaf.py thetaStat do_stat ${S}_theta.thetas.idx" >>${SP}Thetas; }
done
comrad-jobs ${SP}Thetas -n ${SP}Thetas --res theta || stop "thetas failed"
converged $(rtlogs ${SP}Thetas)
cd "$WD"
{ printf 'lineage\tn\tpi\tthetaW\ttajimaD\tnSites'
[ -z "$NSUB" ] || printf '\tnSub\tpiSub\tthetaWSub\ttajimaDSub\tnSitesSub'
printf '\n'
for POP in $POPS; do
n=$(wc -l <angsd/${SP}${POP}Bams)
printf '%s\t%d\t%s' $POP $n "$(tstat theta/${SP}${POP}.thetas.idx.pestPG $n)"
[ -z "$NSUB" ] || printf '\t%d\t%s' $NSUB "$(tstat theta/sub/${SP}${POP}Sub_theta.thetas.idx.pestPG $NSUB)"
printf '\n'
done
} >theta/${SP}Theta.out
cat theta/${SP}Theta.out
awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) if ($i ~ /^nSites/) c[i]; next} {for (i in c) v[$i]} END {exit length(v) > 1}' theta/${SP}Theta.out ||
  note "WARNING: theta site counts differ between lineages or subsamples (nSites columns above), though all used theta/${SP}ThetaSites"

plots theta
