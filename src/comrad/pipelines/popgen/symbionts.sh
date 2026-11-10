inskip && skip "in SKIP"
hasSym || skip "host only (no SYMGENOME)"
cd "$WD/bams"
{ printf 'sample\tchr\treads\n'
for i in *.sym.bam; do
samtools view -q 20 $i | awk -v s=${i%.sym.bam} 'BEGIN{OFS="\t"} {n[$3]++} END {for (c in n) print s, c, n[c]}'
done; } >${SP}SymReads
{ printf 'sample\tgenus\treads\n'
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {g[$1]=$2; next}
  FNR>1 {n[$1 FS (($2 in g) ? g[$2] : "UNMAPPED_" $2)] += $3}
  END {for (k in n) {split(k, a, FS); print a[1], a[2], n[k]}}' "$SYMGENUS" ${SP}SymReads | sort -k1,1 -k2,2
} >${SP}SymGenusReads
u=$(awk '/UNMAPPED_/' ${SP}SymGenusReads | wc -l)
[ $u -eq 0 ] || stop "$u rows with a chromosome missing from SYMGENUS"

plots symbionts
