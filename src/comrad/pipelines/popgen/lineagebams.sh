cd "$WD/angsd"
for POP in $POPS Admixed; do
awk -v p=$POP 'FILENAME==ARGV[1] {d[$1]; next} $2 == p && !($1 in d) {print $1}' ../relatedness/${SP}KinDrop ${SP}BamsClusters >${SP}${POP}Bams
done
for POP in $POPS Admixed; do wc -l ${SP}${POP}Bams; done
t=$(cat $(for POP in $POPS Admixed; do echo ${SP}${POP}Bams; done) | wc -l)
equal $(( t + $(wc -l <../relatedness/${SP}KinDrop) )) $(wc -l <${SP}Bams) "lineage bamlists + kin vs ${SP}Bams"
for POP in $POPS; do [ $(wc -l <${SP}${POP}Bams) -ge $MINPURE ] || stop "$POP has fewer than $MINPURE pure individuals once kin are dropped"; done
kf=kin-free; [ -s ../relatedness/${SP}KinDrop ] || kf="kin kept"
note "pure, $kf: $(for POP in $POPS; do printf '%s %s  ' $POP $(wc -l <${SP}${POP}Bams); done)"
