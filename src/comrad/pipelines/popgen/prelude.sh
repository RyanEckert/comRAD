# lineages, once assigned
_bc=$WD/angsd/${SP}BamsClusters
if [ -s "$_bc" ]; then
export POPS=$(awk '$2 != "Admixed" {print $2}' "$_bc" | sort -u | tr '\n' ' ')
pairs(){ set -- $POPS; while [ $# -gt 1 ]; do a=$1; shift; for b in "$@"; do echo "$a $b"; done; done; }
fi
hasSym(){ [ -n "$SYMGENOME" ]; }
getK(){ [[ $K =~ ^[0-9]+$ ]] || stop "K is not set in $CONFIG (got '$K')"; echo "$K"; }
GENOME_FASTA=$(abs "$GENOME_FASTA"); [ -z "$SYMGENOME" ] || SYMGENOME=$(abs "$SYMGENOME")
[ -z "$SYMGENUS" ] || SYMGENUS=$(abs "$SYMGENUS"); META=$(abs "$META")
[ -z "$EXTRA_READS" ] || EXTRA_READS=$(abs "$EXTRA_READS")
export GENOME_FASTA SYMGENOME SYMGENUS META EXTRA_READS
# plots KIND: the step's quick-look PDF (plotStep.py); a failed plot warns, never stops the run
plots(){ plotStep.py "$1" -s "$SP" --dir "$WD" --meta "$META" --meta-sample "$META_SAMPLE" --meta-groups $META_GROUPS \
  ${PLOT_COLOR:+--color "$PLOT_COLOR"} ${PLOT_SHAPE:+--shape "$PLOT_SHAPE"} ${PLOT_ORDER:+--order $PLOT_ORDER} \
  --levels "$PLOT_LEVELS" --colors "$PLOT_COLORS" --shapes "$PLOT_SHAPES" ${PLOT_KCOLOR:+--kcolor "$PLOT_KCOLOR"} --qdepth "$QDEPTH" --kin "$KIN" ${K:+-k "$K"} "${@:2}" || note "WARNING: $1 plots failed; see run/popgen/$RT_STEP.log"; }
# kevidence PREFIX PCAPREFIX (from ngsAdmix/): sets ev, pc, rk, bic and sup, the reads on K, from
# PREFIX_fit.txt, PREFIX_resid.txt, PREFIX_bic.txt, PREFIX_bestQ.txt and ../pcangsd/PCAPREFIX.admix.<K>.Q (or its log)
kevidence(){
  ev=$(awk -F'\t' '$NF=="<-" {print $1}' $1_fit.txt)
  # PCAngsd's K: its eigenvectors (Velicer's MAP test) + 1, as its --admix uses (Meisner & Albrechtsen 2018)
  pc=$(ls ../pcangsd/$2.admix.*.Q 2>/dev/null | sed 's/.*\.admix\.\([0-9]*\)\.Q$/\1/' | head -1)
  [ -n "$pc" ] || pc=$(awk 'match(tolower($0), /[0-9]+ eigenvector/) {k = substr($0, RSTART, RLENGTH) + 0 + 1} END {if (k != "") print k}' ../pcangsd/$2.log 2>/dev/null)
  bic=$(awk -F'\t' '!/^#/ && $NF=="<-" {print $1}' $1_bic.txt 2>/dev/null)
  rk=$(awk -F'\t' 'NR==1 {for (i = 1; i <= NF; i++) if ($i == "residK") c = i; next} c && $c == "<-" {print $1}' $1_resid.txt 2>/dev/null)
  sup=$(awk -F'\t' 'NR==1 {for (i = 1; i <= NF; i++) {if ($i == "K") a = i; if ($i == "inMode") b = i; if ($i == "runs") c = i}; next}
    {printf "%s%s:%s/%s", (n++ ? " " : ""), $a, $b, $c}' $1_bestQ.txt)
}
