mkdir -p "$WD/demo" && cd "$WD/demo"
for POP in $POPS; do printf '%s\t%s\n' "$POP" "$(wc -l <../angsd/${SP}${POP}Bams)"; done >${SP}PopN.txt
cat ${SP}PopN.txt
note "pairs: $(pairs | tr ' ' '-' | tr '\n' ' ')"
