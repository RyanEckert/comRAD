cd "$WD/filteredReads"
n=$(ls *.filtered.fastq.gz | wc -l)
rates(){ # TAG OUT
{ printf 'sample\talignRate\n'
grep -H 'overall alignment rate' logs/*.$1.bowtie2.log | sed "s|logs/||; s|\.$1\.bowtie2\.log:|\t|; s|% overall alignment rate||"
} >$2
equal $(awk 'NR>1' $2 | wc -l) $n "$2 rows vs samples"
awk -F'\t' 'NR>1 {s+=$2; n++} END {printf "NOTE: '"$2"' mean %.1f%%\n", s/n}' $2
}
rates host ${SP}HostAlignmentRates
hasSym && rates symbiont ${SP}SymbiontAlignmentRates

plots reads
