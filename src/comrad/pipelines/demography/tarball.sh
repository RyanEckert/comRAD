cd "$WD"
comrad report demography --dir "$WD" >/dev/null || note "WARNING: comrad report failed; the tarball has no report"
rm -rf ${SP}Demography ${SP}Demography.tar.gz
mkdir -p ${SP}Demography
for f in demo/${SP}MomentsBest.tsv demo/${SP}MomentsGodambe.tsv demo/${SP}Convergence.txt demo/${SP}Nesting.txt demo/${SP}PopN.txt $CONFIG run/demography/demographySummary.txt; do
need $f; ln -s $PWD/$f ${SP}Demography/
done
ln -s $PWD/demo/${SP}Warnings.txt ${SP}Demography/      # may be empty
for f in run/demography/demographyReport.html run/demography/versions.txt run/demography/demography.conf*.*; do [ -e "$f" ] && ln -s $PWD/$f ${SP}Demography/; done
mkdir -p ${SP}Demography/plots
for f in demo/*/${SP}*.pdf demo/*/${SP}*.png; do [ -e "$f" ] && ln -s $PWD/$f ${SP}Demography/plots/; done
tar -czhf ${SP}Demography.tar.gz ${SP}Demography && rm -r ${SP}Demography
note "results: $(readlink -f ${SP}Demography.tar.gz)"
