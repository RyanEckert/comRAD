cd "$WD/demuxReads" 2>/dev/null || stop "no $WD/demuxReads"
ls *.demux.fastq.gz >/dev/null 2>&1 || stop "no .demux.fastq.gz in $PWD"
mkdir -p logs
>${SP}CutadaptJobs
[ "${CFG_cutadapt_cpus:-1}" -gt 1 ] && CUTJ="-j $CFG_cutadapt_cpus " || CUTJ=""
for f in *.demux.fastq.gz; do s=${f%.demux.fastq.gz}
job "grep -q 'Total written' logs/$s.cutadapt.log" "cutadapt $CFG_cutadapt_args ${CUTJ}--compression-level 6 -o $s.filtered.fastq.gz $f >logs/$s.cutadapt.log" >>${SP}CutadaptJobs
done
comrad-jobs ${SP}CutadaptJobs -n ${SP}Cutadapt --res cutadapt || stop "cutadapt failed"
