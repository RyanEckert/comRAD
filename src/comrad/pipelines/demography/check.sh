tools comrad-jobs realSFS momentsModels.py
python3 -c 'import moments, demes' 2>/dev/null || stop "python cannot import moments and demes (they install with comRAD; reinstall it with the python you run it with)"
[ -e "$WD/run/popgen/fst.done" ] || stop "popgen has not finished its fst step in $WD"
grep -qs '^skipped' "$WD/run/popgen/kin.done" && echo "WARNING: popgen skipped kin: close kin are in the spectra, which skews allele frequencies toward their families"
grep -qs '^skipped' "$WD/run/popgen/outliers.done" && echo "WARNING: popgen skipped outliers: SNPs under selection are in the spectra; the models assume neutral sites"
! grep -q '^skipped' "$WD/run/popgen/fst.done" || stop "popgen skipped fst (SKIP or one lineage): there are no 2D spectra to fit"
[ -n "$(pairs)" ] || stop "one lineage in $WD/angsd/${SP}BamsClusters: the demographic models here need two lineages"
while read a b; do need "$WD/sfs/${SP}$a.$b.sfs" "$WD/sfs/${SP}${a}SFS.saf.idx" "$WD/sfs/${SP}${b}SFS.saf.idx"; done < <(pairs)
