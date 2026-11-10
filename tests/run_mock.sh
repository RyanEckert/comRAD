#!/bin/bash
# end-to-end mock run of reads -> popgen -> demography. Usage: bash run_mock.sh
H=$(cd "$(dirname "$0")" && pwd)
[ -z "$COMRAD_VENV" ] || export PATH=$COMRAD_VENV/bin:$PATH
source $H/setup_mock.sh
P=$koko/2bRAD/proj; mkdir -p $P && cd $P
export FAKE_SBATCH_LOG=$T/sbatch.log; >$FAKE_SBATCH_LOG
comrad init cluster >/dev/null
sed -i 's/^queue = .*/queue = shortq7/; s/^long_queue = .*/long_queue = longq7/; s/^email = .*/email = $EMAIL/; s/^slurm = .*/slurm = constraint=amd/' cluster.config
[ -n "$MOCK_LOCALCOPY" ] && sed -i "s/^local_copy = .*/local_copy = yes/" cluster.config
# sample sheet: 16 barcodes -> 14 sint samples (one admixed, two technical replicates) + 2 of another species
i=0; >sampleID.csv
for p in pool1 pool2 pool3 pool4; do for bc in AAAA CCCC GGGG TTTT; do i=$((i+1))
case $i in 1) s=SA01;; 13) s=ST1.1;; 14) s=ST1.2;; 15) s=X01;; 16) s=X02;; *) s=$(printf 'S%02d' $i);; esac
echo "${p}_$bc,$s" >>sampleID.csv; done; done
comrad init reads >/dev/null
sed -i 's/^BS_PROJECTS = .*/BS_PROJECTS = JA1 JA2/; s/^SAMPLE_SHEET = .*/SAMPLE_SHEET = sampleID.csv/; s/^SPECIES = .*/SPECIES = sint=S xsp=X/' reads.config
comrad check reads && comrad run reads --here | grep -E "STOP|PAUSE|note|FINISHED|skip"

# ---- popgen
cd $P/sint
{ echo "sample,bank,depthZone"; for f in demuxReads/*.demux.fastq.gz; do s=$(basename $f .demux.fastq.gz); n=${s//[^0-9]/}; n=${n:0:2}
  echo "$s,$([ $((10#$n % 3)) -eq 0 ] && echo EFGB || echo WFGB),$([ $((10#$n % 2)) -eq 0 ] && echo Shallow || echo Meso)"; done; } >meta.csv
[ -n "$MOCK_BOM" ] && sed -i '1s/^/\xEF\xBB\xBF/' meta.csv
comrad init popgen >/dev/null
sed -i "s|^SP = .*|SP = sint|; s|^GENOME_FASTA = .*|GENOME_FASTA = $koko/ref/host.fa|; s|^SYMGENOME = .*|SYMGENOME = $koko/ref/sym.fa|;
  s|^SYMGENUS = .*|SYMGENUS = $koko/ref/symChrGenus|; s|^META = .*|META = meta.csv|; s|^TECHREPS = .*|TECHREPS = ST1.1 ST1.2|; s|^MAXK = .*|MAXK = 4|; s|^ADMIX_REPS = .*|ADMIX_REPS = 3|; s|^META_GROUPS = .*|META_GROUPS = bank depthZone|" popgen.config
[ -n "$MOCK_HOSTONLY" ] && sed -i 's|^SYMGENOME = .*|SYMGENOME =|; s|^SYMGENUS = .*|SYMGENUS =|' popgen.config
[ -n "$MOCK_CLONEPAUSE" ] || sed -i "s|^CLONE_PAUSE = .*|CLONE_PAUSE = no|" popgen.config
[ -n "$MOCK_QDEPTH" ] && sed -i "s|^QDEPTH = .*|QDEPTH = $MOCK_QDEPTH|" popgen.config
[ -n "$MOCK_OLDKEYS" ] && sed -i "s|^MINQ = .*|MINQ5 = 0.3|; s|^CLONE_MINQ = .*|CLONE_MINQ5 = 0.4|" popgen.config
[ -n "$MOCK_KINKEEP" ] && sed -i "s|^KIN_DROP = .*|KIN_DROP = no|" popgen.config
[ -n "$MOCK_HETDEPTH" ] && sed -i "s|^HET_DEPTH = .*|HET_DEPTH = $MOCK_HETDEPTH|" popgen.config
[ -n "$MOCK_MINGROUP" ] && sed -i "s|^MIN_GROUP = .*|MIN_GROUP = $MOCK_MINGROUP|" popgen.config
[ -n "$MOCK_THETAN" ] && sed -i "s|^THETA_N = .*|THETA_N = $MOCK_THETAN|" popgen.config
[ -n "$MOCK_METRIC" ] && sed -i "s|^CLONE_METRIC = .*|CLONE_METRIC = $MOCK_METRIC|" popgen.config
[ -n "$MOCK_CLONEH" ] && sed -i "s|^CLONE_CUT = .*|CLONEH = $MOCK_CLONEH|" popgen.config
comrad check popgen && comrad run popgen --here | grep -E "STOP|PAUSE|note|FINISHED|skip"

# ---- K, then the rest of popgen
sed -i "s|^K =.*|K = ${MOCK_K:-2}|" popgen.config
comrad run popgen --here | grep -E "STOP|PAUSE|FINISHED"
comrad summary popgen | grep -E "^(STOP|PAUSE)" && exit 1
[ -s sintResults.tar.gz ] || { echo "no sintResults.tar.gz"; exit 1; }

# ---- demography
[ "${MOCK_K:-2}" -ge 2 ] || exit 0
comrad init demography >/dev/null
sed -i "s|^SP = .*|SP = sint|; s|^MODELS = .*|MODELS = SI IM|; s|^BOOT_JOBS = .*|BOOT_JOBS = 2|; s|^FIT_SEEDS = .*|FIT_SEEDS = 2|" demography.config
comrad check demography >/dev/null && comrad run demography --here | grep -E "STOP|PAUSE|FINISHED"
[ -s sintDemography.tar.gz ] || { echo "no sintDemography.tar.gz"; exit 1; }
echo "MOCK OK"
