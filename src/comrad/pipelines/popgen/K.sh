cd "$WD/ngsAdmix"
kevidence ${SP} ${SP}Pcangsd
note "evalAdmix elbow K = ${ev:-none} (residual RMSD against K)"
note "PCAngsd K = ${pc:-?} (eigenvectors kept by its MAP test + 1; admixture at that K in pcangsd/${SP}Pcangsd.admix.${pc:-K}.Q and ngsAdmix/${SP}Admix.pdf)"
note "residual-block K = ${rk:-none below RESID_MAX $RESID_MAX} (smallest K whose worst cluster-block mean |residual| < $RESID_MAX; plots: ngsAdmix/${SP}Residuals.pdf)"
note "k-means BIC K = ${bic:-?} (lowest BIC of k-means on the PCA scores, as adegenet find.clusters; table in ngsAdmix/${SP}_bic.txt)"
note "runs in the dominant mode per K: $sup"
echo; cat ${SP}_bestQ.txt; echo; cat ${SP}_fit.txt; echo; cat ${SP}_resid.txt 2>/dev/null; echo; cat ${SP}_bic.txt 2>/dev/null
# the evalAdmix elbow is never K = 1 or K = MAXK; whether one cluster fits is read
# from the K = 1 residuals above
if [ -z "$K" ]; then
  [ ! -s K ] || stop "ngsAdmix/K is no longer read; put  K = $(cat K)  in $CONFIG and delete ngsAdmix/K"
  pause "choose K from run/popgen/K.log, pcangsd/${SP}Pca.pdf, ngsAdmix/${SP}Admix.pdf and ngsAdmix/${SP}Residuals.pdf, then set  K = <n>  in $CONFIG"
fi
[[ $K =~ ^[0-9]+$ ]] && [ "$K" -ge 1 ] && [ "$K" -le "$MAXK" ] || stop "K = $K in $CONFIG; it must be a whole number from 1 to MAXK ($MAXK)"
note "K = $K"
