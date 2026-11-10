cd "$WD/sfs"
>${SP}BootJobs
while read a b; do nf=$(nfields $a $b)
for s in $(seq 1 $BOOT_JOBS); do
job "[ \$(awk -v nf=$nf 'NF == nf' ${SP}$a.$b.boot$s.sfs 2>/dev/null | wc -l) -ge $BOOT_REPS ]" \
    "localSaf.py realSFS ${SP}${a}SFS.saf.idx ${SP}${b}SFS.saf.idx -bootstrap $BOOT_REPS -resample_chr 1 -maxIter 1000 -seed $s -P $CFG_boot_cpus >${SP}$a.$b.boot$s.sfs" >>${SP}BootJobs
done
done < <(pairs)
comrad-jobs ${SP}BootJobs -n ${SP}Boot --res boot || stop "bootstrap failed"
converged $(rtlogs ${SP}Boot "logs/${SP}Boot*.e*")
