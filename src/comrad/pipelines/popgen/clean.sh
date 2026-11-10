cd "$WD/bams"
n=$(ls ../filteredReads/*.filtered.fastq.gz | wc -l)
equal $(ls *.bam | grep -v '\.sym\.' | wc -l) $n "host bams vs samples"
hasSym && equal $(ls *.sym.bam | wc -l) $n "symbiont bams vs samples"
rm -f ../filteredReads/*.host.bam ../filteredReads/*.symbiont.bam
mkdir -p "$WD/angsd" && cd "$WD/angsd"
ls ../bams/*.bam | grep -v '\.sym\.' >${SP}AllBams
hasSym && ls ../bams/*.sym.bam >${SP}SymBams
note "$(wc -l <${SP}AllBams) samples mapped"
