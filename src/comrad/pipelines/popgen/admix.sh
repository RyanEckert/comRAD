mkdir -p "$WD/ngsAdmix" && cd "$WD/ngsAdmix"
ngsAdmixLauncher.py -f ../angsd/${SP}Unlinked.beagle.gz --maxK $MAXK -r $ADMIX_REPS -n ${SP} -P $CFG_admix_cpus || stop "ngsAdmixLauncher.py failed"
comrad-jobs ${SP}NgsAdmix -n ${SP}NgsAdmix --res admix -c 'test -s {1}.qopt' -x='-o (\S+)' || stop "NGSadmix failed"
pickQ.py --min 1 --max $MAXK --glob ${SP} --eval ../angsd/${SP}Unlinked.beagle.gz --threads $CFG_evaladmix_cpus || stop "pickQ.py failed"
comrad-jobs ${SP}_evalAdmixJobs -n ${SP}EvalAdmix --res evaladmix -c 'test -s {1}' -x='-o (\S+)' || stop "evalAdmix failed"
pickQ.py --glob ${SP} --fit || stop "pickQ.py --fit failed"

# residual correlations by the numbers: evalAdmix leaves positive residuals among
# individuals that share something the model at that K does not capture. Blocks are
# each individual's main cluster at MAXK (the finest split NGSadmix found), so two
# lineages lumped at a low K show up as their own positive blocks; the sampling groups
# (META_GROUPS) are reported too, as an independent check. Suggested K: the smallest K
# whose worst block mean |residual| is below RESID_MAX.
admixResiduals.py -s $SP -k $MAXK -b ../angsd/${SP}Bams --cut $RESID_MAX --meta "$META" --meta-sample "$META_SAMPLE" \
  --meta-groups $META_GROUPS || stop "residual summary failed; see above"

# k-means on the PCA scores with BIC (DAPC's find.clusters), another read on K for the K step
kmeansBic.py -c ../pcangsd/${SP}Pcangsd.cov -k $MAXK -o ${SP}_bic.txt || note "WARNING: k-means BIC failed; see above"

plots admix
