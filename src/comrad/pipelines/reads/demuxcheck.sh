cd "$(abs "$RAW")"
n=$(ls *.fastq.gz | grep -vc '\.demux\.')
m=$(grep -l 'reads in' logs/*.demux.log | wc -l)
equal $m $n "demux logs vs pools"
{ printf 'pool\trawReads\tkept\tpctKept\n'
grep -h 'reads in' logs/*.demux.log | sed 's/: /\t/; s/ reads in, /\t/; s/ kept (/\t/; s/%)//' | sort -k1,1
} >RawCounts
cat RawCounts
awk -F'\t' 'NR>1 {r+=$2; k+=$3} END {printf "NOTE: %d pools, %d raw, %d kept (%.1f%%)\n", NR-1, r, k, 100*k/r}' RawCounts
low=$(awk -F'\t' -v m=$MIN_KEPT 'NR>1 && $4+0 < m {print $1 " " $4 "%"}' RawCounts)
[ -z "$low" ] || stop "pools below $MIN_KEPT% kept: $(echo $low)"
note "$(ls *.demux.fastq.gz | wc -l) .demux.fastq.gz (pools x barcodes)"

plotStep.py pools --dir "$WD" --pipeline reads || note "WARNING: pool plots failed"
