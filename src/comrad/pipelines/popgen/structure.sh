cd "$WD/angsd"
pruneSnps.py -m ${SP}.mafs.gz -b ${SP}.beagle.gz -o ${SP}Unlinked || stop "pruneSnps.py failed"
note "$(wc -l <${SP}Unlinked.sites) unlinked SNPs"
mkdir -p ../pcangsd && cd ../pcangsd
{ job "ls ${SP}Pcangsd.admix.*.Q >/dev/null 2>&1" "pcangsd -b ../angsd/${SP}Unlinked.beagle.gz --admix --threads $CFG_structure_cpus -o ${SP}Pcangsd"
  job "ls ${SP}FisPcangsd*inbreed* >/dev/null 2>&1" "pcangsd -b ../angsd/${SP}Unlinked.beagle.gz --inbreed-samples --threads $CFG_structure_cpus -o ${SP}FisPcangsd"; } >${SP}PcangsdJobs
comrad-jobs ${SP}PcangsdJobs -n ${SP}Pcangsd --res structure || stop "PCAngsd failed"
grep -E 'Loaded|MAF filtering|converge|eigenvector' ${SP}Pcangsd.log
sed 's|.*/||; s|\.bam$||' ../angsd/${SP}Bams >${SP}Pcangsd.keep
cp ${SP}Pcangsd.keep ${SP}FisPcangsd.keep
pcangsdFisToText.py ${SP}FisPcangsd -o ${SP}Fis.txt || stop "pcangsdFisToText.py failed"

plots pca
