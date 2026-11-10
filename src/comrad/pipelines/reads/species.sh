cd "$(abs "$RAW")"
for p in $SPECIES; do
sp=${p%%=*}; pre=${p#*=}
mkdir -p "$WD/$sp/demuxReads"
n=$(ls ${pre}*.demux.fastq.gz 2>/dev/null | wc -l)
[ $n -eq 0 ] || mv -n ${pre}*.demux.fastq.gz "$WD/$sp/demuxReads/"
note "$sp: $n samples -> $sp/demuxReads"
done
left=$(ls *.demux.fastq.gz 2>/dev/null | wc -l)
[ $left -eq 0 ] || stop "$left .demux.fastq.gz match no prefix in SPECIES; they are still in $PWD"
