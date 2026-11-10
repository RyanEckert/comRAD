# source this: builds a mock HPC in $T (default $TMPDIR/comrad-mock) and puts it first on PATH
T=${T:-${TMPDIR:-/tmp}/comrad-mock}
rm -rf $T; mkdir -p $T/bin $T/koko/ref $T/tmp
H=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
for t in sbatch squeue launcher_creator.py mind fetchReads.py trim2bRAD.py sampleRename.py cutadapt readCounts.py \
  bowtie2Launcher.py bowtie2 splitHostSym.py samtools ngsParalog paralogFilter.py angsd realSFS thetaStat plotQC.py \
  filterBams.py pruneSnps.py subsetBeagle.py pcangsd pcangsdFisToText.py ngsAdmixLauncher.py NGSadmix pickQ.py \
  evalAdmix NGSremix pcangsdOutliers.py momentsModels.py; do ln -sf $H/fake.py $T/bin/$t; done
printf '#!/bin/bash\n[ "$1" = whoami ] && echo me\n' >$T/bin/bs; chmod +x $T/bin/bs
export PATH=$T/bin:$PATH koko=$T/koko EMAIL=test@example.org USER=me TMPDIR=$T/tmp RTJOBS_POLL=1
for g in host sym; do echo ">c1" >$T/koko/ref/$g.fa; touch $T/koko/ref/$g.fa.1.bt2; echo "c1 4 4 60 61" >$T/koko/ref/$g.fa.fai; done
printf 'Sym1\tSymbiodinium\nSym2\tBreviolum\nSym3\tCladocopium\n' >$T/koko/ref/symChrGenus
# momentsModels.py is faked; demography's check only needs these to import
python3 -c 'import moments, demes' 2>/dev/null || { mkdir -p $T/py; touch $T/py/moments.py $T/py/demes.py; export PYTHONPATH=$T/py; }
