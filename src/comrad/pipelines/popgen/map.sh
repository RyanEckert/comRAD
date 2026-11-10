cd "$WD/filteredReads"
ls *.filtered.fastq.gz >/dev/null 2>&1 || stop "no .filtered.fastq.gz in $PWD"
# bt GENOME TAG: <sample>.<tag>.bam and logs/<sample>.<tag>.bowtie2.log
# bowtie2's memory is mostly its index, so with [map] mem empty each task asks for the index size + 50% + 2 GB
idxgb(){ du -cbL "$1".*.bt2* 2>/dev/null | awk 'END {printf "%dGB", $1 / 1e9 * 1.5 + 2.999}'; }
bt(){
local T=${2^}
bowtie2Launcher.py -g "$1" -f .filtered.fastq.gz -n ${SP}${T}MapJobs -s $2 --plain --bam -p $CFG_map_cpus --args "$CFG_map_args" || return 1
comrad-jobs ${SP}${T}MapJobs -n ${SP}${T}Map --res map --mem ${CFG_map_mem:-$(idxgb "$1")} -c "grep -q 'overall alignment rate' logs/{1}.$2.bowtie2.log" -x='-U (\S+)\.filtered\.fastq\.gz'
}
if hasSym; then
bt "$SYMGENOME" symbiont & ps=$!
bt "$GENOME_FASTA" host & ph=$!
wait $ps || stop "symbiont mapping failed"
wait $ph || stop "host mapping failed"
else
bt "$GENOME_FASTA" host || stop "host mapping failed"
fi
