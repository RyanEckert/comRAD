cd "$WD"
comrad report popgen --dir "$WD" >/dev/null || note "WARNING: comrad report failed; the tarball has no report"
o(){ local s; for s; do inskip $s && { echo -; return; }; done; }   # files from a step in SKIP are left out (an older run may have left copies)
# folder in the tarball (. = top level), then the files that go there; a leading + marks a file that may be absent, - one left out
FILES=". angsd/${SP}*Bams angsd/${SP}BamsClusters +angsd/${SP}Groups angsd/${SP}SnpDepth.txt $(o kin)relatedness/${SP}Kin relatedness/${SP}KinDrop $(o symbionts)+bams/${SP}SymGenusReads filteredReads/${SP}*AlignmentRates filteredReads/${SP}RetainedRates filteredReads/${SP}ReadCounts
qc angsd/${SP}Depth.pdf angsd/${SP}Depth.info angsd/${SP}Quality*X.txt
clones angsd/clones/${SP}CloneBams angsd/clones/${SP}Clones.ibsMat angsd/clones/${SP}CloneGroups angsd/clones/${SP}CloneDend.pdf +angsd/clones/${SP}CloneMetric.txt +angsd/clones/${SP}CloneBorderline.txt +angsd/clones/${SP}CloneReplicates.txt +angsd/clones/${SP}CloneDendIBS.pdf +angsd/${SP}CloneFreeDend.pdf
plots +filteredReads/${SP}Reads*.pdf +bams/${SP}Symbionts*.pdf +pcangsd/${SP}Pca*.pdf +ngsAdmix/${SP}Admix*.pdf +relatedness/${SP}Kin*.pdf +angsd/lineages/${SP}Outliers*.pdf +ngsAdmix/${SP}Neutral*.pdf +sfs/${SP}Fst*.pdf +het/${SP}Het*.pdf +theta/${SP}Theta*.pdf
structure angsd/${SP}.beagle.gz angsd/${SP}Outliers $(o outliers)angsd/lineages/${SP}*Selection.tsv angsd/${SP}Unlinked.sites angsd/${SP}.ibsMat angsd/lineages/${SP}*.ibsMat pcangsd/${SP}*Pcangsd.cov pcangsd/${SP}Fis.txt ngsAdmix/${SP}_k*_best.qopt ngsAdmix/${SP}_k*_best.corres ngsAdmix/${SP}_bestQ.txt ngsAdmix/${SP}_fit.txt $(o neutral)ngsAdmix/${SP}Neutral_best.qopt $(o neutral)ngsAdmix/${SP}NeutralCompare.txt $(o neutral)+ngsAdmix/${SP}Neutral_bestQ.txt $(o neutral)+ngsAdmix/${SP}Neutral_fit.txt $(o neutral)+ngsAdmix/${SP}Neutral_resid.txt $(o neutral)+ngsAdmix/${SP}Neutral_bic.txt ngsAdmix/${SP}_resid.txt +ngsAdmix/${SP}_bic.txt +ngsAdmix/${SP}_residBlocks.txt ngsAdmix/${SP}Residuals.pdf
relatedness $(o kin)relatedness/${SP}Remix
diversity $(o het)het/${SP}.het $(o het)het/${SP}Ho.txt $(o het)het/${SP}HetDepth.txt $(o het)+angsd/${SP}HetExcluded $(o het)+het/${SP}HetMeanDepth.txt $(o theta)theta/${SP}Theta.out $(o theta)+theta/sub/${SP}*SubBams $(o fst)+sfs/${SP}LineageFst.out $(o fstall fst outliers)+sfs/${SP}LineageFstAll.out $(pairs | while read a b; do printf "$(o fst)sfs/${SP}$a.$b.sfs "; done)$(o groups)+angsd/groups/${SP}GroupFst.out $(o groups)+angsd/groups/${SP}GroupTheta.out
run run/popgen/popgenSummary.txt $CONFIG +run/popgen/popgenReport.html +run/popgen/versions.txt +run/popgen/popgen.conf*.*"
rm -rf ${SP}Results ${SP}Results.tar.gz
bad=0
while read -r dir files; do
mkdir -p ${SP}Results/$dir
for g in $files; do
case $g in -*) continue;; esac
opt=0; case $g in +*) opt=1; g=${g#+};; esac
for f in $g; do
if [ ! -e "$f" ]; then [ $opt -eq 1 ] && continue; echo "missing: $f"; bad=1; continue; fi
if [ -e ${SP}Results/$dir/$(basename $f) ]; then echo "name clash: $f"; bad=1; continue; fi
ln -s $PWD/$f ${SP}Results/$dir/
done
done
done <<<"$FILES"
[ $bad -eq 0 ] || stop "files missing from the tarball; see above"
tar -czhf ${SP}Results.tar.gz ${SP}Results && rm -r ${SP}Results
tar -tzf ${SP}Results.tar.gz | grep -v '/$' | sed 's|^[^/]*/||' | awk -F/ '{print (NF>1 ? $1 : ".")}' | sort | uniq -c
note "results: $(readlink -f ${SP}Results.tar.gz)"
