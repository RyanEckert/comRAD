cd "$WD/ngsAdmix"
K=$(getK)
awk -v t=$ASSIGN '
  NR==FNR {b[FNR]=$0; nb=FNR; next}
  {m=0; j=0; for(i=1;i<=NF;i++) if($i>m){m=$i; j=i}
   lin[FNR] = (m < t ? "Admixed" : "Cluster" j); nq=FNR}
  END {if (nb != nq) {printf "bamlist %d rows, qopt %d rows\n", nb, nq >"/dev/stderr"; exit 1}
       for(i=1;i<=nb;i++) print b[i], lin[i]}
' ../angsd/${SP}Bams ${SP}_k${K}_best.qopt >../angsd/${SP}BamsClusters || stop "bamlist and qopt differ"
awk '{print $2}' ../angsd/${SP}BamsClusters | sort | uniq -c
note "lineages: $(awk '{print $2}' ../angsd/${SP}BamsClusters | sort | uniq -c | awk '{printf "%s %s  ", $2, $1}')"
awk -v m=$MINPURE '$2 != "Admixed" {n[$2]++} END {for (c in n) if (n[c] < m) {print c " has " n[c]; bad=1}; exit bad}' ../angsd/${SP}BamsClusters ||
  stop "a lineage has fewer than $MINPURE pure individuals"
if [ "$K" -eq 1 ]; then
  note "K = 1: every sample is one lineage (Cluster1); lineage FST is skipped and demography will not run"
else
  [ $(awk '$2 != "Admixed" {print $2}' ../angsd/${SP}BamsClusters | sort -u | wc -l) -ge 2 ] || stop "fewer than 2 lineages with pure individuals at K = $K"
fi
