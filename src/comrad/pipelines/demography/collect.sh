cd "$WD/demo"
n=$(ls */*.s[0-9]*.best.tsv 2>/dev/null | wc -l)
equal $n $(( $(pairs | wc -l) * $(echo $MODELS | wc -w) * FIT_SEEDS )) "moments results vs pairs x models x seeds"
H=$(head -1 $(ls */*.s[0-9]*.best.tsv | head -1))
# best seed per pair and model; its run is labelled s<seed>.<run>
for f in */*.s[0-9]*.best.tsv; do k=${f%.best.tsv}; k=${k##*.}; awk -F'\t' -v OFS='\t' -v k=$k 'FNR>1 {$4 = k "." $4; print}' $f; done |
  sort -t$'\t' -k1,1 -k2,2 -k3,3 -k5,5gr | awk -F'\t' '!s[$1 FS $2 FS $3]++' >best.tmp
{ echo "$H"
awk -F'\t' -v OFS='\t' -v h="$H" '
  BEGIN {n = split(h, c, "\t"); for (i = 1; i <= n; i++) col[c[i]] = i}
  {r[NR] = $0; p = $1 FS $2; a = $col["AIC"] + 0; cl = $col["CLAIC"]
   if (!(p in ma) || a < ma[p]) ma[p] = a
   if (cl != "" && (!(p in mc) || cl + 0 < mc[p])) mc[p] = cl + 0}
  END {for (i = 1; i <= NR; i++) {
         $0 = r[i]; p = $1 FS $2
         $col["deltaAIC"] = sprintf("%.4f", $col["AIC"] - ma[p])
         if ($col["CLAIC"] == "") {$col["deltaCLAIC"] = ""; key = 1e300}
         else {$col["deltaCLAIC"] = sprintf("%.4f", $col["CLAIC"] - mc[p]); key = $col["deltaCLAIC"]}
         print key, $0}}' best.tmp | sort -t$'\t' -k2,2 -k3,3 -k1,1g | cut -f2-
} >${SP}MomentsBest.tsv
rm best.tmp
cut -f1-5,7-11 ${SP}MomentsBest.tsv
while IFS=$'\t' read p1 p2 M rest; do note "best model $p1-$p2: $M"; done < <(awk -F'\t' 'NR>1 && !s[$1 FS $2]++' ${SP}MomentsBest.tsv)

# convergence: runs and seeds within 1 ll unit of the best; CHECK = fewer than 2 seeds
for f in */*.s[0-9]*.tsv; do case $f in *.best.tsv|*.godambe.tsv) continue;; esac
k=${f%.tsv}; k=${k##*.}; awk -F'\t' -v k=$k 'FNR>1 && $5 ~ /^-?[0-9]/ {print $1"."$2, $3, k, $5}' $f; done |
awk -v S=$FIT_SEEDS '{p = $1 " " $2; ll[NR] = $4; key[NR] = p; sd[NR] = $3; n[p]++; if (!(p in b) || $4 > b[p]) b[p] = $4}
  END {for (i = 1; i <= NR; i++) if (ll[i] >= b[key[i]] - 1) {w[key[i]]++; if (!((key[i], sd[i]) in ss)) {ss[key[i], sd[i]]; ns[key[i]]++}}
       for (p in n) printf "%s\t%s\tbest %.2f\t%d of %d runs within 1\t%d of %d seeds\n", (ns[p] < 2 ? "CHECK" : "ok   "), p, b[p], w[p], n[p], ns[p], S}' | sort -k2,2 -k3,3 >${SP}Convergence.txt
cat ${SP}Convergence.txt
c=$(grep -c '^CHECK' ${SP}Convergence.txt); [ $c -eq 0 ] || note "$c models found their best by only one seed (${SP}Convergence.txt); consider more seeds"

# nesting: a model must fit at least as well as one it contains
awk -F'\t' 'NR==FNR {ll[$1 "." $2, $3] = $5; pr[$1 "." $2]; next}
  {for (p in pr) if (((p, $1) in ll) && ((p, $2) in ll)) {d = ll[p, $2] - ll[p, $1]
     printf "%s\t%s\t%s in %s\t%+.2f\n", (d < -.1 ? "CHECK" : "ok   "), p, $1, $2, d}}' <(tail -n +2 ${SP}MomentsBest.tsv) - <<'NEST' | sort -k2,2 >${SP}Nesting.txt
SI	IM
IM	IMa
IM	AM
IM	SC
IM	IMghost
SC	SCghost
SI	SI2N
IM	IM2N
IMa	IMa2N
AM	AM2N
SC	SC2N
IM	IM2m
AM	AM2m
SC	SC2m
IM2N	IM2m2N
IM2m	IM2m2N
AM2N	AM2m2N
AM2m	AM2m2N
SC2N	SC2m2N
SC2m	SC2m2N
NEST
cat ${SP}Nesting.txt
c=$(grep -c '^CHECK' ${SP}Nesting.txt); [ $c -eq 0 ] || note "$c nested pairs where the larger model fits worse (${SP}Nesting.txt): its optimiser stopped early"

# bound and Godambe warnings from the kept seed's log
tail -n +2 ${SP}MomentsBest.tsv | while IFS=$'\t' read p1 p2 M run rest; do
grep -h "at a bound\|Godambe failed" $p1.$p2/${SP}$p1.$p2.$M.${run%%.*}.log 2>/dev/null | sed "s/^/$p1.$p2\t/"
done | tee ${SP}Warnings.txt
{ echo -e "pop1\tpop2\tmodel\tparam\tvalue\tse"
tail -n +2 ${SP}MomentsBest.tsv | while IFS=$'\t' read p1 p2 M run rest; do
awk -F'\t' -v m=$M 'FNR>1 && $3 == m' $p1.$p2/${SP}$p1.$p2.$M.${run%%.*}.godambe.tsv 2>/dev/null
done; } >${SP}MomentsGodambe.tsv
