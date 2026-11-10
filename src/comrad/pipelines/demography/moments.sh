cd "$WD/demo"
>${SP}MomentsJobs
while read a b; do
na=$(awk -v p=$a '$1==p {print $2}' ${SP}PopN.txt); nb=$(awk -v p=$b '$1==p {print $2}' ${SP}PopN.txt)
mkdir -p $a.$b
for M in $MODELS; do for k in $(seq 1 $FIT_SEEDS); do
o=$a.$b/${SP}$a.$b.$M.s$k
job "test -s $o.best.tsv" "export OMP_NUM_THREADS=$CFG_moments_cpus; momentsModels.py -i ../sfs/${SP}$a.$b.sfs -n $na $nb -p $a $b -m $M -r $FIT_RUNS --maxiter $MAXITER --seed $k --mu $MU --boot ../sfs/${SP}$a.$b.boot.sfs ${PROJECT:+--project $PROJECT} --plot -o $o >$o.log 2>&1" >>${SP}MomentsJobs
done; done
done < <(pairs)
note "$(wc -l <${SP}MomentsJobs) moments jobs"
comrad-jobs ${SP}MomentsJobs -n ${SP}Mom --res moments || stop "moments jobs failed"
