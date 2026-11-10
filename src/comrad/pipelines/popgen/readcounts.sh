cd "$WD/demuxReads"
miss=$(for f in *.demux.fastq.gz; do [ -s ${f%.demux.fastq.gz}.filtered.fastq.gz ] || echo $f; done | wc -l)
[ $miss -eq 0 ] || stop "$miss .demux.fastq.gz have no .filtered.fastq.gz"
if [ -n "$EXTRA_READS" ]; then
for f in "$EXTRA_READS"/*.filtered.fastq.gz; do [ -e "$(basename "$f")" ] || ln -s "$f" .; done
note "$(ls "$EXTRA_READS"/*.filtered.fastq.gz | wc -l) extra .filtered.fastq.gz from $EXTRA_READS"
fi
readCounts.py -f filtered.fastq.gz -j ${SLURM_CPUS_ON_NODE:-8} -o ${SP}Read || stop "readCounts.py failed"
awk 'NR>1 {n++; t+=$2} END {printf "NOTE: %d samples, %d reads, %.0f mean\n", n, t, t/n}' ${SP}ReadCounts
mkdir -p ../filteredReads/logs
mv *.filtered.fastq.gz ${SP}ReadCounts ../filteredReads/
mv logs/*.cutadapt.log ../filteredReads/logs/ 2>/dev/null
rm -f *.demux.fastq.gz
