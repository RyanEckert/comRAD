cd "$WD/filteredReads"
mkdir -p ../bams
n=$(ls *.filtered.fastq.gz | wc -l)
if hasSym; then
splitHostSym.py scan -H '*.host.bam' -S '*.symbiont.bam' --outdir ../bams/bed -n ${SP}ScanJobs || stop "splitHostSym.py scan failed"
comrad-jobs ${SP}ScanJobs -n ${SP}Scan --res split -c 'test -s logs/{1}.scanlog' -x='--sample (\S+)' || stop "scan failed"
{ printf 'sample\thostAligned\tsymAligned\tbothRefs\thostLoci\thostBp\tsymLoci\tsymBp\n'; grep -hP '^\S+\t\d+\t' logs/*.scanlog | sort -k1,1; } >${SP}CrossMapCounts

splitHostSym.py split -H '*.host.bam' -S '*.symbiont.bam' --bed ../bams/bed --outdir ../bams --sort -n ${SP}SplitJobs || stop "splitHostSym.py split failed"
comrad-jobs ${SP}SplitJobs -n ${SP}Split --res split -c 'test -s logs/{1}.splitlog && test -s ../bams/{1}.bam' -x='--sample (\S+)' || stop "split failed"

# logs left from an earlier run whose samples no longer exist
awk '{for (i=1; i<=NF; i++) if ($i == "--sample") print $(i+1)}' ${SP}SplitJobs | sort >${SP}WantedSamples
ls logs/*.splitlog | sed 's|logs/||; s|\.splitlog$||' | sort | comm -13 ${SP}WantedSamples - >${SP}StaleSamples
while read s; do rm -f logs/$s.splitlog ../bams/$s.bam ../bams/$s.bam.bai ../bams/$s.sym.bam ../bams/$s.sym.bam.bai; done <${SP}StaleSamples
[ -s ${SP}StaleSamples ] && note "removed $(wc -l <${SP}StaleSamples) stale samples"

{ printf 'sample\thostAligned\tsymAligned\tbothRefs\thostWritten\tsymWritten\thostMasked\tsymMasked\tmaskLoci\tmaskBp\n'; grep -hP '^\S+\t\d+\t' logs/*.splitlog | sort -k1,1; } >${SP}SplitCounts
equal $(awk 'NR>1' ${SP}SplitCounts | wc -l) $n "split logs vs samples"
awk 'BEGIN{print "sample\tinputReads\thostAligned%\tsymAligned%\thostKept\thostKept%\tsymKept\tsymKept%"}
     NR==FNR{if(FNR>1){s=$1; sub(/\.filtered\.fastq\.gz$/,"",s); tot[s]=$2} next}
     FNR>1 {n=tot[$1]; if(!n){print $1"\tno read count" > "/dev/stderr"; next}
      printf "%s\t%d\t%.1f\t%.1f\t%d\t%.1f\t%d\t%.1f\n", $1, n, 100*$2/n, 100*$3/n, $5, 100*$5/n, $6, 100*$6/n}' \
    ${SP}ReadCounts ${SP}SplitCounts >${SP}RetainedRates
else
>${SP}SortJobs
for f in *.host.bam; do s=${f%.host.bam}
job "test -s ../bams/$s.bam.bai" "samtools sort -l 9 -o ../bams/$s.bam $f && samtools index ../bams/$s.bam" >>${SP}SortJobs
done
comrad-jobs ${SP}SortJobs -n ${SP}Sort --res split || stop "sort failed"
{ printf 'sample\taligned\n'; for i in ../bams/*.bam; do s=$(basename $i .bam); printf '%s\t%s\n' $s "$(samtools view -c $i)"; done; } >${SP}AlignedCounts
awk 'BEGIN{print "sample\tinputReads\taligned\taligned%"}
     NR==FNR{if(FNR>1){s=$1; sub(/\.filtered\.fastq\.gz$/,"",s); tot[s]=$2} next}
     FNR>1 {n=tot[$1]; if(!n){print $1"\tno read count" > "/dev/stderr"; next}
      printf "%s\t%d\t%d\t%.1f\n", $1, n, $2, 100*$2/n}' ${SP}ReadCounts ${SP}AlignedCounts >${SP}RetainedRates
fi
