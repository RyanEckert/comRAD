cd "$WD/sfs"
while read a b; do nf=$(nfields $a $b)
for s in $(seq 1 $BOOT_JOBS); do
n=$(awk -v nf=$nf 'NF == nf' ${SP}$a.$b.boot$s.sfs 2>/dev/null | wc -l); [ $n -ge $BOOT_REPS ] || stop "$a.$b seed $s has $n of $BOOT_REPS replicates"
done
cat $(seq -f "${SP}$a.$b.boot%g.sfs" 1 $BOOT_JOBS) | awk -v nf=$nf 'NF == nf' >${SP}$a.$b.boot.sfs
note "$a.$b: $(grep -c . ${SP}$a.$b.boot.sfs) bootstrap spectra"
done < <(pairs)
