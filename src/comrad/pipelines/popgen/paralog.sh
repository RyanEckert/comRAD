cd "$WD/angsd"
zcat ${SP}Depth.pos.gz | awk 'NR>1 {print $1"\t"$2}' >${SP}DepthSites
mkdir -p paralog
rm -f paralog/chunk_*
split -l 20000 -d -a 4 ${SP}DepthSites paralog/chunk_
MI=$(mind ${SP}Bams ${CFG_paralog_min_ind:-0.5})
>${SP}ParalogJobs
for c in paralog/chunk_*; do
job "grep -q 'Finished!' $c.log" "( samtools mpileup -b ${SP}Bams -l $c -q 0 -Q 0 -B -d 100000 --ff UNMAP,DUP | ngsParalog calcLR -infile - -outfile $c.lr -minQ 20 -minind $MI -mincov 1 ) > $c.log 2>&1" >>${SP}ParalogJobs
done
comrad-jobs ${SP}ParalogJobs -n ${SP}Paralog --res paralog || stop "ngsParalog failed"
paralogFilter.py --collect paralog -s ${SP}DepthSites -o ${SP} || stop "paralogFilter.py failed; see above"
wc -l ${SP}DepthSites ${SP}ParalogSites ${SP}KeepSites
note "$(wc -l <${SP}ParalogSites) paralog sites removed of $(wc -l <${SP}DepthSites)"
sitesindex ${SP}KeepSites
