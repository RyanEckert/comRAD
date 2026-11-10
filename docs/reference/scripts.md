# Scripts

comRAD installs these scripts on your `PATH`. The pipelines call them, and each one also works on its
own. Below is every script's `--help`, as `tools/makeDocs.py` captured it.

| script | what it does |
|---|---|
| [`trim2bRAD.py`](#trim2bradpy) | Trim, deduplicate and demultiplex 2bRAD reads by in-line barcode, or write a list of trim commands to run in parallel. |
| [`sampleRename.py`](#samplerenamepy) | Rename files based on columns of data from a specified .csv file. |
| [`fetchReads.py`](#fetchreadspy) | Download pooled Illumina reads from BaseSpace, merge lane parts without decompressing, and flatten to one file per pool. |
| [`readCounts.py`](#readcountspy) | Prints a tab-separated file of Illumina read counts for files with a specified extension. |
| [`bowtie2Launcher.py`](#bowtie2launcherpy) | Prints a list of bowtie2 commands, one per read file (for comrad-jobs or launcher_creator.py). |
| [`splitHostSym.py`](#splithostsympy) | Split host and symbiont 2bRAD alignments of the same trimmed reads into two disjoint sets, dropping every read that aligned to both references. |
| [`packMetagenome.py`](#packmetagenomepy) | Pack genomes into a metagenome of balanced, per-genome named chromosomes. |
| [`plotQC.py`](#plotqcpy) | Plot the ANGSD depth/quality QC (a Python port of ngsTools' plotQC.R) and write each sample's proportion of covered sites at depth &gt;= D, for filtering samples by coverage. |
| [`filterBams.py`](#filterbamspy) | Drop low-quality bams and/or prune clone groups to their best member. |
| [`paralogFilter.py`](#paralogfilterpy) | FDR-filter ngsParalog calcLR output and write the sites to keep for ANGSD -sites. |
| [`pruneSnps.py`](#prunesnpspy) | One SNP per 2bRAD tag (highest MAF), outlier tags removed, and the matching cut-down beagle. |
| [`pcangsdFisToText.py`](#pcangsdfistotextpy) | Collect PCAngsd per-individual inbreeding (--inbreed-samples) into one TSV, paired with bam names from NAME.keep. |
| [`ngsAdmixLauncher.py`](#ngsadmixlauncherpy) | Prints a list of NGSadmix commands for a range of K (for comrad-jobs or launcher_creator.py). |
| [`pickQ.py`](#pickqpy) | pickQ.py -- pick the best ngsAdmix run within the dominant mode, per K, and score model fit against K from the evalAdmix residuals. |
| [`pcangsdOutliers.py`](#pcangsdoutlierspy) | Per-SNP p-values and FDR outliers from a PCAngsd --selection scan. |
| [`subsetBeagle.py`](#subsetbeaglepy) | Subset an ANGSD beagle genotype-likelihood file (-doGlf 2) to some individuals. |
| [`momentsModels.py`](#momentsmodelspy) | Fit two-population divergence models to a folded realSFS spectrum with moments; optional Godambe CLAIC and SEs from block-bootstrap spectra. |
| [`zipper.py`](#zipperpy) | Prints a list of gzip/gunzip commands to be used with launcher_creator.py. |
| [`mind`](#mind) | mind &lt;bamlist&gt; [proportion]  -- prints ceiling(nlines * proportion). |
| [`cloneGroups.py`](#clonegroupspy) | Builds an average-linkage (UPGMA) tree of the samples from one of two distances, cuts it into clone groups, checks the technical replicates and draws the tree. |
| [`admixResiduals.py`](#admixresidualspy) | evalAdmix leaves positive residual correlations among individuals that share something the admixture model at that K does not capture. |
| [`cloneTestRun.py`](#clonetestrunpy) | Submits one ANGSD job per variant of the clone step, each in its own folder under --outdir, with the IBS matrix and genotype probabilities that cloneTestScore.py compares. |
| [`cloneTestScore.py`](#clonetestscorepy) | For each variant folder made by cloneTestRun.py, builds two distances between samples: ANGSD's IBS (1 - IBS) and genotype discordance (the share of sites where two samples' genotype calls differ, among sites both called with probability above --conf), and scores how well each separates clones. |
| [`simCloneKin.py`](#simclonekinpy) | Simulates samples from two lineages, each with unrelated founders, full-sib families whose parents are also sampled (parent-offspring pairs), half-sibs, selfed sibs, clones (extra ramets of a genet) and technical replicates (the same colony sequenced again). |
| [`plotStep.py`](#plotsteppy) | Draws comRAD's quick-look plots for one step from the files that step wrote, as one PDF per step (several panels or pages where there is more than one thing to show), plus PNGs of each page in run/&lt;pipeline&gt;/plots/ for comrad report. |
| [`kmeansBic.py`](#kmeansbicpy) | Clusters samples by k-means on their principal components and scores each K by BIC, as adegenet's find.clusters does before DAPC (Jombart et al. |

## trim2bRAD.py

```text
usage: trim2bRAD.py [-h] [-i FASTQ] [-f EXT] [--site RE] [--bc2 RE]
                    [--adaptor RE] [--degen RE] [--clip N] [--sample-id N]
                    [--min-bc-count N] [--gzip] [--out-ext EXT] [--level 1-9]
                    [--no-dedup] [-n FILE] [--logdir DIR] [--log-ext EXT]
                    [--launcher] [-N N] [-w N] [-q QUEUE] [-t TIME] [-e ADDR]

Trim, deduplicate and demultiplex 2bRAD reads by in-line barcode, or write a list of trim commands to run in parallel.

options:
  -h, --help         show this help message and exit

mode (give exactly one):
  -i, --input FASTQ  trim this one fastq (plain or .gz)
  -f, --suffix EXT   write one trim command per file ending in EXT in the
                     working directory, e.g. 'fastq' or '.fastq.gz'

trimming:
  --site RE          restriction site pattern (default: BcgI)
  --bc2 RE           in-line barcode following the fragment (default:
                     [ATGC]{4})
  --adaptor RE       far-end adaptor (default: AGATC?)
  --degen RE         degenerate 5' leader (default: [ATGC][ATGC][AG][AT]CC)
  --clip N           bases to remove from each end of the fragment (default:
                     0)
  --sample-id N      1-based position of the name-deriving field in the
                     underscore-separated filename; default is the filename up
                     to its first dot
  --min-bc-count N   write a barcode's file only if it has at least N reads
                     (default: 100000)
  --gzip             gzip the output (<out-ext>.gz); cutadapt reads it
                     directly
  --out-ext EXT      output name ending after <sample>_<barcode> (default:
                     .demux.fastq; '.tr0' gives the 2bRAD_denovo names)
  --level 1-9        gzip compression level for --gzip (default: 5). 9 costs
                     about 6x the CPU of 5 for ~5% less size
  --no-dedup         keep PCR duplicates; output then retains the degenerate
                     leader, as in the Perl

job list:
  -n, --name FILE    commands file to write (default: trims)
  --logdir DIR       per-file logs, <sample><log-ext> (default: logs), '' to
                     disable
  --log-ext EXT      log name ending (default: .demux.log)
  --launcher         also build the .slurm with launcher_creator.py
  -N, --nodes N      nodes for launcher
  -w, --wayness N    commands per node for launcher
  -q, --queue QUEUE  SLURM partition
  -t, --time TIME    walltime
  -e, --email ADDR   address for SLURM notifications

examples:
  # trim one file; sample name is the first underscore field
  trim2bRAD.py -i Sample_L001_R1.fastq --sample-id 1

  # gzipped in and out, keep any barcode with at least 1000 reads
  trim2bRAD.py -i Sample_L001_R1.fastq.gz --sample-id 1 --min-bc-count 1000 --gzip

  # AlfI instead of BcgI
  trim2bRAD.py -i Sample_L001_R1.fastq --site '.{12}GCA.{6}TGC.{12}'

  # clip 2 bases off each end, keep PCR duplicates
  trim2bRAD.py -i Sample_L001_R1.fastq --clip 2 --no-dedup

  # write one trim command per .fastq here, then submit them
  trim2bRAD.py -f fastq --sample-id 1 -n trims --launcher -e you@fau.edu

  # write the commands file only, submit it yourself
  trim2bRAD.py -f .fastq.gz --sample-id 1 -n trims --gzip

  # collect the per-sample yields afterwards
  grep 'total goods' logs/*.demux.log
```

## sampleRename.py

```text
usage: sampleRename.py [-h] [-f EXT] [-i CSV] [-n PREFIX] [--header]
                       [--dry-run]

Rename files based on columns of data from a specified .csv file.

options:
  -h, --help  show this help message and exit
  -f EXT      file extension, not including the leading '.'; two-part
              extensions such as tr0.gz are fine [default: fastq]
  -i CSV      input .csv, with or without the .csv suffix, in the working
              directory [default: sampleList]
  -n PREFIX   common name given to all samples, preceding the sample number
              from the .csv
  --header    skip the first row of the .csv
  --dry-run   report the plan and rename nothing

examples:
  # rename trimmed, demultiplexed reads to sample names
  sampleRename.py -i sampleNames -f tr0.gz -n MGM_

  # see the whole plan first, change nothing
  sampleRename.py -i sampleNames -f tr0.gz -n MGM_ --dry-run

  # csv has a header row
  sampleRename.py -i sampleNames -f tr0.gz -n MGM_ --header

csv:
  Two columns, no quoting needed: the existing name without its extension, and
  the sample number.

    FG6-6_ACCA,001
    FG6-6_TGGT,002

  With -f tr0.gz -n MGM_ that renames FG6-6_ACCA.tr0.gz to MGM_001.tr0.gz.
  -f takes the whole extension without the leading dot, so two-part
  extensions like 'tr0.gz' or 'trim.gz' work.

  A new name may be another row's old name (renumbering, or swapping two
  samples); those are renamed through temporary names so nothing is
  overwritten.
```

## fetchReads.py

```text
usage: fetchReads.py [-h] [-p NAME [NAME ...]] [-o DIR] [-e EXT]
                     [--undetermined] [--merge-reads] [--keep-parts]
                     [--dry-run]

Download pooled Illumina reads from BaseSpace, merge lane parts without decompressing, and flatten to one file per pool.

options:
  -h, --help            show this help message and exit
  -p, --project NAME [NAME ...]
                        BaseSpace project name(s) to download, each into its
                        own subdirectory; omit to work on files already here
  -o, --outdir DIR      directory to download into and flatten into (default:
                        cwd)
  -e, --ext EXT         read file extension; also the only files downloaded
                        (default: .fastq.gz)
  --undetermined        also download the Undetermined reads, which are left
                        out by default
  --merge-reads         concatenate lane parts into one file per pool; without
                        it, pools with more than one part are left unprocessed
                        and reported
  --keep-parts          keep the per-lane files after merging instead of
                        deleting them once the merged size checks out
  --dry-run             report what would happen and change nothing

examples:
  # download two projects and flatten to one file per pool
  fetchReads.py -p JA23123 JA23124

  # same, and combine lane parts where a pool has them
  fetchReads.py -p JA23123 JA23124 --merge-reads

  # also download the Undetermined reads
  fetchReads.py -p JA26157 --undetermined

  # submit it as a single launcher job
  launcher_creator.py -b 'srun fetchReads.py -p JA23123 JA23124' -n fetchReads \
      -q shortq7 -t 06:00:00 -e you@fau.edu
  sbatch fetchReads.slurm

  # already downloaded; just merge and flatten what is there
  fetchReads.py --merge-reads

  # show the plan and change nothing
  fetchReads.py --dry-run

  # read counts afterwards, straight off the gz
  readCounts.py -f fastq.gz -o gom -j 16

what is downloaded:
  Only files ending in -e (default .fastq.gz). A BaseSpace project also
  holds the demultiplexer's reports, logs, InterOp and sample sheet; those are
  not fetched. Undetermined reads are left out unless --undetermined is given.
  The filter is passed to bs as --exclude/--include, which needs BaseSpace CLI
  1.4.0 or later.

pools, not samples:
  A file here is a pool. Samples inside it are separated later by in-line
  barcode, so identity is pool plus barcode and two pools must never be
  combined. That is why identically named files in different projects get the
  project appended instead of being merged: they are either the same file
  twice or two different pools, and merging the second case would put two
  specimens under one barcode.

failures:
  A problem with one pool prints a reason, leaves that pool's files at the top
  level under their original Illumina names, and the run carries on. Long
  names mark them as unprocessed. The exit status is non-zero if anything was
  skipped. Only a bad --outdir, a failed download, or no read files at all
  stop the run.
```

## readCounts.py

```text
usage: readCounts.py [-h] [-f F] [-o O] [-j J] [-l L]

Prints a tab-separated file of Illumina read counts for files with a specified extension.
Counts newlines in binary chunks and reads files in parallel, and handles gzipped input.

Output is 'file<TAB>reads' under a header line, so it reads straight into R
or awk:
  read.delim('trimCounts')
  awk 'NR>1 {t+=$2} END {print t}' trimCounts

usage examples:
  readCounts.py                                 # *fq in the working directory -> readCounts
  readCounts.py -f fastq -o raw                 # *fastq -> rawCounts
  readCounts.py -f fq.gz -j 8                   # gzipped reads, 8 files at a time
  readCounts.py -f trim -l 1 -o trim            # one line per read instead of four

options:
  -h, --help  show this help message and exit
  -f F        File extension for files to add, don't need to include '.' [e.g.
              fastq]
  -o O        Base name for Counts file [e.g. read -> readCounts]
  -j J        Number of files to count in parallel [default = 8]
  -l L        Lines per read [default = 4 for fastq; use 1 for one-line-per-
              read formats]
```

## bowtie2Launcher.py

```text
usage: bowtie2Launcher.py [-h] [-f F] [-g G] [-s S] [-p THREADS] [--args ARGS]
                          [--plain] [--samdir SAMDIR] [--bam] [-l L]
                          [--logdir LOGDIR] [-n N] [--launcher] [-N N] [-w W]
                          [-t T] [-e E] [-q Q]

Prints a list of bowtie2 commands, one per read file (for comrad-jobs or launcher_creator.py). Default alignment settings (--no-unal --score-min L,16,1 --local -L 16) are tuned for 36 bp 2bRAD tags.

options:
  -h, --help            show this help message and exit
  -f F                  File type to map, don't include [e.g. trim]
  -g G                  Path to reference genome
  -s S                  Optional tag inserted into the SAM and log names, e.g.
                        -s sym gives SAMPLE.trim.sym.bt2.sam and
                        SAMPLE.trim.sym.bt2log. Use a different -s and/or
                        --samdir per reference so two mappings of the same
                        reads don't overwrite each other
  -p, --threads THREADS
                        bowtie2 threads per sample (bowtie2 -p) [default: 1]
  --args ARGS           bowtie2 alignment options, replacing the default '--
                        no-unal --score-min L,16,1 --local -L 16'
  --plain               Name outputs <sample>.<tag>.bam (or .sam) and logs
                        <sample>.<tag>.bowtie2.log, where <sample> is the file
                        name without -f; needs -s
  --samdir SAMDIR       Optional directory for alignment output [default: cwd]
  --bam                 If used, pipe bowtie2 through samtools and write .bam
                        instead of .sam. Same size as a gzip -9 SAM but needs
                        no decompressing, and splitHostSym.py reads it
                        directly
  -l L                  BAM compression level for --bam [default: 9]
  --logdir LOGDIR       Directory for per-sample bowtie2 stderr logs [default:
                        logs], '' to disable
  -n N                  Name of output file to write the commands to
  --launcher            If used, create .slurm file with launcher_creator.py
  -N N                  With --launcher: the number of nodes to split job
                        across in slurm
  -w W                  With --launcher: wayness, the number of commands to
                        run per node
  -t T                  With --launcher: walltime [default: 6:00:00]
  -e E                  With --launcher: email to send SLURM info to
  -q Q                  With --launcher: partition to submit job to

examples:
  # map every .trim in the working directory
  bowtie2Launcher.py -g $GENOME -f .trim -n maps --logdir logs

  # write BAM instead of SAM: same size as a gzip -9 SAM, but usable as-is
  bowtie2Launcher.py -g $GENOME -f .trim -n maps --bam --logdir logs

  # gzipped reads work too; bowtie2 decompresses them itself
  bowtie2Launcher.py -g $GENOME -f .trim.gz -n maps --bam --logdir logs

  # 4 bowtie2 threads per sample (ask the scheduler for 4 CPUs per task as well)
  bowtie2Launcher.py -g $GENOME -f .trim.gz -n maps --bam -p 4
  comrad-jobs maps -n maps --cpus 4 --mem 8GB

  # comRAD names: <sample>.host.bam and logs/<sample>.host.bowtie2.log, where the
  # sample is the file name without -f
  bowtie2Launcher.py -g $GENOME -f .filtered.fastq.gz -n hostMaps -s host --plain

  # different alignment settings (replaces the 2bRAD defaults entirely)
  bowtie2Launcher.py -g $GENOME -f .trim.gz -n maps --bam --args '--no-unal --local --very-sensitive-local'

  # two independent mappings of the SAME reads, for splitHostSym.py.
  # -s tags each pass so the two don't overwrite each other; --samdir is
  # optional on top of that.
  bowtie2Launcher.py -g $SYMGENOME  -f .trim -n symMaps  -s sym --bam
  bowtie2Launcher.py -g $HOSTGENOME -f .trim -n hostMaps -s host --bam

  # pull the alignment rates out of the logs afterwards
  grep -H 'overall alignment rate' logs/*.bt2log | sed 's|logs/||; s|\\.bt2log:| |; s|% overall alignment rate||'

  # rerun only the jobs whose log has no alignment rate (TAG is the -s value,
  # or '' if -s was not used)
  awk -v tag=.sym '{if(!match($0, /-U [^ ]+/)) next;
        s=substr($0, RSTART+3, RLENGTH-3); sub(/\\.gz$/, "", s);
        f="logs/" s tag ".bt2log"; ok=0;
        while((getline l < f)>0) if(l ~ /overall alignment rate/) ok=1;
        close(f); if(!ok) print}' symMaps > symMapsRerun

index:
  -g is the bowtie2 index basename, which for these genomes is the .fasta
  path itself: bowtie2-build writes <basename>.1.bt2 .. .rev.2.bt2 beside it
  (.bt2l for a large index). The job list is not written unless those exist,
  since bowtie2 would otherwise exit before emitting a single SAM line and
  every task would fail with samtools reporting 'fail to read the header'.

threads:
  -p sets bowtie2's own thread count. The scheduler's CPU request only reserves
  cores; bowtie2 uses one thread unless -p says otherwise.

exit status:
  With --bam the command is run through 'bash -o pipefail' so that a bowtie2
  that dies mid-run fails the task. Without it the pipeline reports samtools'
  status, and samtools exits 0 on a truncated stream -- leaving a valid but
  incomplete BAM, which splitHostSym.py cannot distinguish from reads that
  genuinely did not map.
```

## splitHostSym.py

```text
usage: splitHostSym.py [-h] COMMAND ...

Split host and symbiont 2bRAD alignments of the same trimmed reads into two disjoint sets, dropping every read that aligned to both references.

positional arguments:
  COMMAND
    scan      find the cross-mapping loci and write them as BEDs
    split     write the two disjoint bams, optionally masked

options:
  -h, --help  show this help message and exit

the two subcommands:
  scan    find the cross-mapping loci, all samples, no alignments written
  split   write the partitioned bams, optionally masking what scan found

  They cannot be one command: masking is a reduction over ALL samples, so
  split for one sample needs every other sample's BEDs, and launcher tasks are
  independent processes that cannot wait for each other.

examples:
  # partition only, every sample, sorted and indexed
  splitHostSym.py split -H '*.trim.host.bt2.sam' -S 'symbiont/*.trim.sym.bt2.sam' \
      --outdir final --sort --launcher -e you@fau.edu

  # the default path: find cross-mapping loci, then split with them masked
  splitHostSym.py scan  -H '*.trim.host.bt2.sam' -S 'symbiont/*.trim.sym.bt2.sam' \
      --outdir bed --launcher -e you@fau.edu
  splitHostSym.py split -H '*.trim.host.bt2.sam' -S 'symbiont/*.trim.sym.bt2.sam' \
      --bed bed --outdir final --sort --launcher -e you@fau.edu

  # one sample, right now
  splitHostSym.py split -H '*.trim.host.bt2.bam' -S '*.trim.sym.bt2.bam' \
      --sample XXX13.1

  # collect the summaries
  cat logs/*.scanlog  | sort -k1,1
  cat logs/*.splitlog | sort -k1,1

sample names:
  The sample name is whatever the '*' in the pattern covers: the file name
  with the literal text after the last '*' removed. With -H
  '*.trim.host.bt2.bam',
    OGM001.trim.host.bt2.bam  -> OGM001
    XXX13.1.trim.host.bt2.bam -> XXX13.1
  The patterns must be quoted. If the shell expands them the script never sees
  the suffix, and it stops with an error rather than guess.

output names:
  OGM001 -> OGM001.bam (+ .bai with --sort) and OGM001.sym.bam
  scan writes OGM001.host.bed and OGM001.sym.bed. The host bams keep the bare
  sample name because they are the primary analysis product; only the
  symbiont side is tagged. List them with
    ls *.bam | grep -v '\.sym\.' > hostBams
    ls *.sym.bam > symBams

pairing:
  -H and -S are paired by sample name. A name matched by only one pattern is
  listed and skipped, the rest still run, and the exit status is non-zero. Two
  files giving the same name is an error, since their outputs would overwrite
  each other. A single pair, or --sample, runs immediately; more than one pair
  writes a job list.

before running either:
  Both mappings must have finished for every sample. A bowtie2 task that died
  mid-write leaves reads missing from its SAM, and a read missing from the
  symbiont SAM is indistinguishable from a host-only read. Check with
    ls *.trim | wc -l
    grep -lc 'overall alignment rate' logs/*.sym.bt2log  | wc -l
    grep -lc 'overall alignment rate' logs/*.host.bt2log | wc -l
```

## packMetagenome.py

```text
usage: packMetagenome.py [-h] -i INPUTS [INPUTS ...] -o OUTPUT -m MAP
                         [--prefixes PREFIXES [PREFIXES ...]]
                         --chrs-per-genome CHRSPERGENOME [--pad PAD]
                         [--plain-names] [--gap-seq GAPSEQ]

Pack genomes into a metagenome of balanced, per-genome named chromosomes.

options:
  -h, --help            show this help message and exit
  -i, --inputs INPUTS [INPUTS ...]
                        input FASTA files, one per genome, .gz ok (order =
                        chromosome order)
  -o, --output OUTPUT   output metagenome FASTA
  -m, --map MAP         output TSV mapping contigs to chromosomes
  --prefixes PREFIXES [PREFIXES ...]
                        chromosome name prefix per input, same order as -i
                        [e.g. Sym Brev Clad Dur]
  --chrs-per-genome CHRSPERGENOME
                        number of chromosomes to split each genome into (>=1)
  --pad PAD             zero-padding width for the chromosome number
  --plain-names         name chromosomes chr1, chr2, ... instead of by prefix
  --gap-seq GAPSEQ      gap sequence inserted between contigs

chromosome naming:
  Chromosomes are named <prefix><number>, e.g. Sym001, Sym002, Brev001. The
  prefix comes from --prefixes, one per input genome in the same order as -i.
  Two inputs may share a prefix (two Cladocopium genomes, say); numbering then
  continues across them rather than restarting. Without --prefixes the input
  file stem is used, and with --plain-names the old chr1, chr2, ... naming is
  kept.

  The point of the naming is that samtools idxstats output can be summed per
  genome by name alone. The map TSV carries the same grouping in its genome and
  prefix columns, so either can be used downstream.

outputs:
  <output>   metagenome FASTA, one record per chromosome
  <map>      TSV: chr, prefix, genome, contig_id, contig_length,
             contig_start_on_chr, contig_end_on_chr

examples:
  # four Symbiodiniaceae genera, five chromosomes each
  packMetagenome.py -i Symbiodinium.fasta Breviolum.fasta Cladocopium.fasta Durusdinium.fasta \
      --prefixes Sym Brev Clad Dur --chrs-per-genome 5 \
      -o symMetaGenome.fasta -m symMetaGenome.map.tsv

  # two Cladocopium assemblies share a prefix and number consecutively
  packMetagenome.py -i CladA.fa.gz CladB.fa.gz --prefixes Clad Clad --chrs-per-genome 2 \
      -o clad.fasta -m clad.map.tsv

  # then index it
  bowtie2-build symMetaGenome.fasta symMetaGenome.fasta
  samtools faidx symMetaGenome.fasta
```

## plotQC.py

```text
usage: plotQC.py [-h] -p PREFIX -b BAMS [-q D [D ...]] [-x XL] [-c D [D ...]]
                 [--chunk CHUNK] [--no-counts]

Plot the ANGSD depth/quality QC (a Python port of ngsTools' plotQC.R) and write each sample's proportion of covered sites at depth >= D, for filtering samples by coverage.

options:
  -h, --help            show this help message and exit
  -p, --prefix PREFIX   prefix of the angsd -doDepth output
  -b, --bams BAMS       bam list, in the order angsd saw them
  -q, --quality D [D ...]
                        minimum depths, inclusive; one quality<D>X.txt per
                        value (default 5)
  -x, --xl XL           max depth to plot (default 20)
  -c, --coverage D [D ...]
                        depth cutoffs to mark on the depth panel; one or more
                        (default 10)
  --chunk CHUNK         rows per chunk when streaming counts
  --no-counts           skip the .counts pass (panels 3 and 4)

examples:
  plotQC.py -p sintDepth -b sintAllBams                # sintDepth.pdf, quality5X.txt, 10x marked
  plotQC.py -p sintDepth -b sintAllBams -q 3 5 10      # quality3X.txt, quality5X.txt, quality10X.txt
  plotQC.py -p sintDepth -b sintAllBams -c 5 10 20     # mark three depths on the per-sample panel
  plotQC.py -p sintDepth -b sintAllBams --no-counts    # skip the slow .counts panels

inputs:
  ANGSD run with -doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2 -out <prefix>, which writes
  <prefix>.qs, <prefix>.depthSample, <prefix>.depthGlobal and <prefix>.counts.gz.

thresholds:
  both -q and -c are inclusive: 5 means depth >= 5. quality<D>X.txt lists each sample's
  proportion of its COVERED sites at >= D, worst first; filterBams.py -q reads it.
```

## filterBams.py

```text
usage: filterBams.py [-h] -b BAMS -q QUALITY [-p PROP] [-c CLONES] [-o OUT]

Drop low-quality bams and/or prune clone groups to their best member.

options:
  -h, --help            show this help message and exit
  -b, --bams BAMS       bam list, one file/path per line
  -q, --quality QUALITY
                        two columns: sample, value
  -p, --prop PROP       minimum value to keep; -p 0.25 -> >=0.25
  -c, --clones CLONES   clone groups file, one group per line
  -o, --out OUT         output list (default: <bams>.q)
```

## paralogFilter.py

```text
usage: paralogFilter.py [-h] (-i LR [LR ...] | --collect DIR) -s FILE
                        -o PREFIX [--fdr Q]

FDR-filter ngsParalog calcLR output and write the sites to keep for ANGSD -sites.

options:
  -h, --help            show this help message and exit
  -i, --input LR [LR ...]
                        calcLR output file(s)
  --collect DIR         directory of chunk_NNNN site lists with .lr and .log
                        beside them; checks all finished and combines them
                        into <o>.lr first
  -s, --sites FILE      site list calcLR was run on, scaffold <TAB> position
  -o, --out PREFIX      output prefix
  --fdr Q               Benjamini-Hochberg threshold (default .05)

examples:
  # one calcLR run
  paralogFilter.py -i mcav.lr -s mcavddSites -o mcav

  # chunked run: check every chunk finished, combine into mcav.lr, then filter
  paralogFilter.py --collect paralog -s mcavddSites -o mcav

  # stricter
  paralogFilter.py -i mcav.lr -s mcavddSites -o mcav --fdr 0.01

  # then, for ANGSD
  angsd sites index mcavKeepSites

outputs:
  <o>.lr             combined calcLR output (--collect only)
  <o>ParalogSites    flagged sites, scaffold <TAB> position
  <o>KeepSites       -s sites minus the flagged ones, original order
  <o>KeepRegions     scaffolds in <o>KeepSites, for ANGSD -rf

p-values:
  p = 0.5 * P(chi2_1 > LR), Benjamini-Hochberg adjusted across every site
  together. calcLR fits each site on its own reads, so chunking does not change
  the LR; the correction is the only step that needs all sites at once.
```

## pruneSnps.py

```text
usage: pruneSnps.py [-h] -m FILE -b FILE [-x FILE] [--tag BP] [--window BP]
                    -o PREFIX

One SNP per 2bRAD tag (highest MAF), outlier tags removed, and the matching cut-down beagle.

options:
  -h, --help          show this help message and exit
  -m, --mafs FILE     ANGSD .mafs(.gz) from the same run as the beagle
  -b, --beagle FILE   ANGSD .beagle.gz
  -x, --exclude FILE  chr <TAB> pos of SNPs whose tags are dropped (e.g.
                      pcangsdOutliers.py output)
  --tag BP            tag length; SNPs closer than this are on one tag
                      (default 36)
  --window BP         drop SNPs within this many bp of an excluded SNP
                      (default 35)
  -o, --out PREFIX    output prefix

examples:
  # all individuals: drop outlier tags, one SNP per tag
  pruneSnps.py -m mcav.mafs.gz -b mcav.beagle.gz -x mcavOutliers -o mcavUnlinked

  # one lineage, same outliers
  pruneSnps.py -m lineages/mcavCluster1.mafs.gz -b lineages/mcavCluster1.beagle.gz \
      -x mcavOutliers -o lineages/mcavCluster1Unlinked

  # linkage only, no outlier removal
  pruneSnps.py -m mcav.mafs.gz -b mcav.beagle.gz -o mcavThin

output:
  <out>.beagle.gz   the beagle, header plus the kept SNPs, in beagle order
  <out>.sites       chr <TAB> pos of the kept SNPs
```

## pcangsdFisToText.py

```text
usage: pcangsdFisToText.py [-h] -o FILE [--lineages FILE]
                           [--keep-suffix KEEP_SUFFIX]
                           NAME [NAME ...]

Collect PCAngsd per-individual inbreeding (--inbreed-samples) into one TSV, paired with bam names from NAME.keep.

positional arguments:
  NAME                  PCAngsd -o prefix; NAME.keep (bam names, one per line,
                        in beagle order) must exist alongside

options:
  -h, --help            show this help message and exit
  -o, --out FILE        output TSV
  --lineages FILE       bam in column 1, lineage in column 2 (e.g.
                        mcavBamsClusters); adds a lineage column
  --keep-suffix KEEP_SUFFIX
                        suffix of the bam-name file (default .keep)

examples:
  # one run: bam, Fis
  pcangsdFisToText.py mcavPcangsd -o mcavFis.txt

  # add each individual's lineage from the NGSadmix assignment
  pcangsdFisToText.py mcavPcangsd --lineages ../angsd/mcavBamsClusters -o mcavFis.txt

  # several runs: bam, run, Fis
  pcangsdFisToText.py mcavCluster1Pcangsd mcavCluster2Pcangsd -o mcavLineageFis.txt

output columns:
  bam, then run (only with more than one run), then lineage (only with
  --lineages), then Fis.
```

## ngsAdmixLauncher.py

```text
usage: ngsAdmixLauncher.py [-h] [-f F] [--minK MINK] [--maxK MAXK] [-n N]
                           [-P THREADS] [--minMaf MINMAF] [--launcher] [-e E]
                           [-r R]

Prints a list of NGSadmix commands for a range of K (for comrad-jobs or launcher_creator.py)

options:
  -h, --help            show this help message and exit
  -f F                  Input beagle file name
  --minK MINK           Minium value of K
  --maxK MAXK           Maximum value of K
  -n N                  Common name of outputs
  -P, --threads THREADS
                        NGSadmix threads per run (NGSadmix -P) [default: 10]
  --minMaf MINMAF       NGSadmix -minMaf (default 0: keep all beagle sites so
                        .fopt matches the beagle for evalAdmix/NGSremix)
  --launcher            If used, create .slurm file with launcher_creator.py
  -e E                  If used, email to send SLURM info to
  -r R                  Number of replicate simulations to run per K

examples:
  ngsAdmixLauncher.py -f ../angsd/apal.beagle.gz --maxK 7 -r 20 -n apal
      beagle already MAF-filtered by ANGSD; NGSadmix keeps every site (-minMaf 0)
  ngsAdmixLauncher.py -f ../angsd/apal.beagle.gz --maxK 7 -r 20 -n apal --minMaf 0.05
      apply NGSadmix's own MAF filter (its built-in default is 0.05)
  ngsAdmixLauncher.py -f ../angsd/apal.beagle.gz --maxK 7 -r 20 -n apal -P 4
  comrad-jobs apalNgsAdmix -n apalNgsAdmix --cpus 4 --mem 8GB
      4 NGSadmix threads per run; ask the scheduler for the same number of CPUs

threads:
  -P sets NGSadmix's own thread count. The scheduler's CPU request only reserves
  cores; NGSadmix uses what -P says.
```

## pickQ.py

```text
usage: pickQ.py [-h] [--min MIN] [--max MAX] [--thr THR] [--glob GLOB]
                [--out OUT] [--eval BEAGLE] [--minMaf MINMAF] [--fit]
                [--threads THREADS] [--no-copy]

pickQ.py -- pick the best ngsAdmix run within the dominant mode, per K,
and score model fit against K from the evalAdmix residuals.

For each K, all replicate runs are compared pairwise after optimally matching
their clusters (Hungarian assignment, exact at any K). Runs are grouped into
modes by average-linkage clustering on that similarity, and the highest-
likelihood run in the largest mode is chosen.

Similarity is 1 - ||Q1 - Q2||^2 / (2n): 1 for identical solutions, 0 for
maximally different ones.

  pickQ.py --min 1 --max 15 --glob 2bradNoClones
  pickQ.py --min 1 --max 15 --glob 2bradNoClones --eval 2bradNoClones.beagle.gz
  pickQ.py --glob 2bradNoClones --fit

options:
  -h, --help         show this help message and exit
  --min MIN          minimum K (default 1)
  --max MAX          maximum K (default 10)
  --thr THR          similarity above which two runs are the same mode
                     (default 0.9)
  --glob GLOB        prefix glob for Q files (default '*')
  --out OUT          table to write (default <glob>_bestQ.txt, or
                     <glob>_fit.txt with --fit)
  --eval BEAGLE      also write an evalAdmix job list using this beagle file
  --minMaf MINMAF    -minMaf for the evalAdmix jobs; must match the NGSadmix
                     runs (default 0, as ngsAdmixLauncher.py)
  --fit              second pass: score model fit against K from the
                     *_best.corres residual matrices and mark the elbow
  --threads THREADS  -P for the evalAdmix jobs
  --no-copy          don't copy winners to *_best.qopt

examples:
  # pass 1: modes, best run per K, and the evalAdmix job list
  pickQ.py --min 1 --max 10 --glob gommcav --eval mcavNoClones.beagle.gz

  # same, for NGSadmix runs made with ngsAdmixLauncher.py --minMaf 0.05
  pickQ.py --min 1 --max 10 --glob gommcav --eval mcavNoClones.beagle.gz --minMaf 0.05

  # submit that job list and wait for it
  launcher_creator.py -j gommcav_evalAdmixJobs -n evalAdmix -q shortq7 -t 6:00:00 -e $EMAIL
  sbatch evalAdmix.slurm

  # pass 2: model fit against K, from the residual matrices those jobs wrote
  pickQ.py --glob gommcav --fit

why two passes:
  --fit needs the *_best.corres files evalAdmix produces, which do not exist
  when the job list is written, so it cannot run in the same invocation. It
  skips the mode analysis entirely -- that work is already on disk -- and
  reads inMode back out of <glob>_bestQ.txt for the support column.

--minMaf:
  evalAdmix filters the beagle with its own -minMaf (built-in default 0.05)
  and stops if the .fopt has a different number of sites. It must match the
  -minMaf the NGSadmix runs used; ngsAdmixLauncher.py defaults to 0, and so
  does this.

the elbow:
  Kneedle on the RMSD of the off-diagonal residual correlations: the K
  furthest from the straight line joining the first and last points. Same
  quantity and same formula as admix_knee() in admixfit.R, so the two cannot
  disagree on the same input. Residual RMSD is used rather than the
  likelihood because it plateaus once K captures the structure, while the
  likelihood keeps climbing with every added cluster.

  It is a geometric heuristic, not a test. It always returns something, it
  moves if the range of K changes, and on a shallow curve it is close to
  arbitrary -- read it next to the inMode column, and next to
  CLUMPAK/StructureSelector.
```

## pcangsdOutliers.py

```text
usage: pcangsdOutliers.py [-h] -b FILE -s PREFIX [--fdr Q] -o PREFIX

Per-SNP p-values and FDR outliers from a PCAngsd --selection scan.

options:
  -h, --help         show this help message and exit
  -b, --beagle FILE  the beagle PCAngsd was run on
  -s, --scan PREFIX  PCAngsd -o prefix; needs PREFIX.selection and
                     PREFIX.sites (--sites-save)
  --fdr Q            q-value cutoff (default .05)
  -o, --out PREFIX   output prefix

examples:
  # PCAngsd run with --selection --sites-save -o mcavSelPcangsd
  pcangsdOutliers.py -b ../angsd/mcav.beagle.gz -s mcavSelPcangsd -o mcav

  # stricter cutoff
  pcangsdOutliers.py -b ../angsd/mcav.beagle.gz -s mcavSelPcangsd --fdr .01 -o mcav

output:
  <out>Selection.tsv   chr, pos, then stat / p / q for each PC, every SNP
                       PCAngsd tested
  <out>Outliers        chr <TAB> pos of SNPs with q < --fdr on any PC

p = P(chi-square, 1 df > stat), as in the PCAngsd README. q is
Benjamini-Hochberg within each PC.
```

## subsetBeagle.py

```text
usage: subsetBeagle.py [-h] -b BEAGLE -s SAMPLES -k KEEP -o OUT

Subset an ANGSD beagle genotype-likelihood file (-doGlf 2) to some individuals.

options:
  -h, --help            show this help message and exit
  -b, --beagle BEAGLE   input beagle (.gz or plain)
  -s, --samples SAMPLES
                        every sample in the beagle, one per line, in bamlist
                        order
  -k, --keep KEEP       the samples to keep, one per line, in output order
  -o, --out OUT         output beagle; gzipped if it ends in .gz

examples:
  subsetBeagle.py -b sintNeutral.beagle.gz -s ../angsd/sintBams -k ../angsd/sintNeutralBams -o sintNeutralNoKin.beagle.gz
      the neutral beagle with kin left out (both lists are bam paths: names must match exactly)
  subsetBeagle.py -b xesto.beagle.gz -s samples.txt -k subsample01.keep -o subsample01.beagle.gz

notes:
  --samples lists every individual in the beagle, in bamlist order; --keep lists the ones to
  keep, in the order to write them. The order written is echoed to stderr and <out>.order.
  The output is gzipped if its name ends in .gz.
```

## momentsModels.py

```text
usage: momentsModels.py [-h] -i FILE -n N N [-p NAME NAME] -o PREFIX
                        [-m M [M ...]] [--ghost] [-r N] [--maxiter N]
                        [--fold F] [--mu RATE] [-L N] [--project N N]
                        [--boot FILE] [--eps F] [--seed SEED] [--plot]

Fit two-population divergence models to a folded realSFS spectrum with moments; optional Godambe CLAIC and SEs from block-bootstrap spectra.

options:
  -h, --help            show this help message and exit
  -i, --sfs FILE        realSFS 2d spectrum for the two populations
  -n, --sizes N N       individuals per population, in the order given to
                        realSFS
  -p, --pops NAME NAME  population names
  -o, --out PREFIX      output prefix
  -m, --models M [M ...]
                        models to fit (default: all non-ghost models)
  --ghost               also fit the ghost-population models
  -r, --runs N          perturbed starting points per model (default 20)
  --maxiter N           optimizer iterations per run (default 100)
  --fold F              how far to perturb each start, in factors of two
                        (default 1)
  --mu RATE             mutation rate per site per generation; without it the
                        scaled columns are left empty
  -L, --sites N         sites the spectrum represents (default: its total,
                        monomorphic included)
  --project N N         down-project to this many individuals per population
                        before folding; applied to --boot spectra too
  --boot FILE           realSFS -bootstrap output (one spectrum per line);
                        adds CLAIC and Godambe SEs for each model's best run
  --eps F               fractional step for Godambe derivatives (default .01)
  --seed SEED           rng seed (default 1)
  --plot                write data/model/residual panels per model

models:
  SI SI2N IM IM2N IM2m IM2m2N IMa IMa2N AM AM2N AM2m AM2m2N SC SC2N SC2m SC2m2N
  ghost (with --ghost): IMghost SCghost

examples:
  # every non-ghost model, 20 restarts each, coalescent units only
  momentsModels.py -i mcavC1.C2.sfs -n 51 15 -p C1 C2 -o C1.C2.all

  # one model per job, with Godambe CLAIC/SEs, scaled units and fit plots
  momentsModels.py -i mcavC1.C2.sfs -n 51 15 -p C1 C2 -m SC2m -r 50 --mu 2e-8 \
      --boot mcavC1.C2.boot.sfs --plot -o C1.C2.SC2m

  # ghost models only, down-projected to keep the 3-population integration fast
  momentsModels.py -i mcavC1.C2.sfs -n 51 15 -p C1 C2 -m IMghost SCghost \
      --project 20 10 --boot mcavC1.C2.boot.sfs -o C1.C2.ghost

output:
  <out>.tsv           one row per model per run
  <out>.best.tsv      best run of each model, sorted by CLAIC if --boot, else AIC
  <out>.godambe.tsv   with --boot: value and Godambe SE of every parameter
                      (coalescent units; theta included)
  <out>.<M>.pdf       with --plot: data / model / residual panels

notes:
  -n is INDIVIDUALS, in the order given to realSFS.
  deltaAIC / deltaCLAIC are within this invocation. When models are fitted in
  separate jobs, recompute them after combining the .best.tsv files.
  Convergence is the spread of ll across runs: if the top runs of a model
  disagree, raise -r before reading its parameters.
```

## zipper.py

```text
usage: zipper.py [-h] [-f F [F ...]] [--gunzip] [-a A] [-n N] [--launcher]
                 [-e E] [-w W]

Prints a list of gzip/gunzip commands to be used with launcher_creator.py

options:
  -h, --help    show this help message and exit
  -f F [F ...]  File type(s) to gzip, don't include '.' [e.g. fastq tr0 sam]
  --gunzip      If used, gunzip all .gz files
  -a A          Optional arguments to pass to g/unzip
  -n N          Name of output file to write launcher commands to
  --launcher    If used, create .slurm file with launcher_creator.py
  -e E          If used, email to send SLURM info to
  -w W          Wayness for launcher

examples:
  # gzip every .fastq in the current directory, maximum compression
  zipper.py -f fastq -a -9 --launcher -e you@fau.edu

  # several file types at once
  zipper.py -f fastq tr0 trim sam -a -9 --launcher -e you@fau.edu

  # gunzip everything ending in .gz
  zipper.py --gunzip --launcher -e you@fau.edu

  # write the commands file only, then submit it yourself
  zipper.py -f sam -a -9 -n zipSams
```

## mind

```text
mind <bamlist> [proportion]  -- prints ceiling(nlines * proportion)
```

## cloneGroups.py

```text
usage: cloneGroups.py [-h] -i IBS [-g GENO] -b BAMS [-r TECHREPS]
                      [--metric {auto,geno,ibs,both}] [--cut CUT]
                      [--conf CONF [CONF ...]] [--min-shared MIN_SHARED]
                      [--min-well MIN_WELL] [--max-tight MAX_TIGHT]
                      [--window WINDOW] [--buffer BUFFER] [--maxdrop MAXDROP]
                      [--check CUT] [--quality QUALITY] [--meta META]
                      [--meta-sample META_SAMPLE]
                      [--meta-groups [META_GROUPS ...]] [--color COLOR]
                      [--shape SHAPE] [--levels LEVELS] [--colors COLORS]
                      [--shapes SHAPES] -o OUT [--dir DIR]

Builds an average-linkage (UPGMA) tree of the samples from one of two distances, cuts it into clone groups, checks the technical replicates and draws the tree. The distances: IBS (1 - IBS from ANGSD -doIBS 1: one random read per sample and site) and genotype distance (the share of sites where two samples' genotype calls differ, among sites both called with probability above a confidence level, from ANGSD -doGeno 8 -doPost 1). With --metric auto the genotype distance is used when the data support it (enough shared confident sites for most pairs, and replicates tight compared with typical distances), at 0.99 confidence, else 0.95; otherwise IBS.

options:
  -h, --help            show this help message and exit
  -i, --ibs IBS         ANGSD .ibsMat (1 - IBS), one row per bam
  -g, --geno GENO       ANGSD .geno.gz with genotype probabilities (-doGeno 8
                        or 15, -doPost 1); without it, IBS only
  -b, --bams BAMS       the bam list ANGSD read, in the same order
  -r, --techreps TECHREPS
                        technical replicates, space-separated <sample>.<n>
                        names
  --metric {auto,geno,ibs,both}
                        distance (default auto)
  --cut CUT             auto, gap, replicates or a number (default auto)
  --conf CONF [CONF ...]
                        genotype call confidence levels, tried in order
                        (default 0.99 0.95)
  --min-shared MIN_SHARED
                        shared confident sites for a pair to count as well
                        measured (default 1000)
  --min-well MIN_WELL   share of pairs that must be well measured (default
                        0.90)
  --max-tight MAX_TIGHT
                        largest replicate distance / median distance allowed
                        (default 0.15)
  --window WINDOW       gap cut: gaps must start within this x the replicates'
                        largest distance (default 3)
  --buffer BUFFER       treeTrimmer cut: the replicates' deepest join x this
                        (default 1.10)
  --maxdrop MAXDROP     exit 4 if the groups would drop more than this
                        proportion of the samples
  --check CUT           check mode: tree of a clone-free set with this cut;
                        warn for pairs below it
  --quality QUALITY     <SP>Quality<N>X.txt (bam, proportion of sites at >= N
                        x), for coverage in the borderline report
  --meta META           sample metadata CSV, for the dendrogram's tip colours
                        and shapes
  --meta-sample META_SAMPLE
                        metadata column with the sample names (default sample)
  --meta-groups [META_GROUPS ...]
                        metadata columns: the first fills the tips, the second
                        shapes them
  --color COLOR         metadata column that colours the tips (default: first
                        --meta-groups)
  --shape SHAPE         metadata column that shapes the tips (default: second
                        --meta-groups)
  --levels LEVELS       level order per column, as PLOT_LEVELS: 'bank: WFGB
                        EFGB; depthZone: Shallow Mesophotic'
  --colors COLORS       colours per level, as PLOT_COLORS: 'bank: WFGB=D53E4F
                        EFGB=F88D51' (hex without #, or names)
  --shapes SHAPES       tip shapes per level, as PLOT_SHAPES: 'depthZone:
                        Shallow=^ Mesophotic=v'
  -o, --out OUT         output prefix, e.g. the species code
  --dir DIR             folder for the outputs (default here)

examples:
  cloneGroups.py -i ofavClones.ibsMat -g ofavClones.geno.gz -b ofavCloneBams -r "OGM024.1 OGM024.2 OGM066.1 OGM066.2" -o ofav
      automatic choice of distance and cut; writes ofavCloneGroups, ofavCloneReplicates.txt,
      ofavCloneMetric.txt, ofavCloneBorderline.txt and ofavCloneDend.pdf (plus the other distance's tree)
  cloneGroups.py -i ofavClones.ibsMat -b ofavCloneBams -r "OGM024.1 OGM024.2" --metric ibs --buffer 1.10 -o ofav
      IBS only, cut at the replicates' deepest join x 1.10 (treeTrimmer)
  cloneGroups.py -i pastClones.ibsMat -g pastClones.geno.gz -b pastCloneBams -r "$TECHREPS" --cut 0.075 -o past
      a cut you chose, in the units of the distance used (here the automatic one)
  cloneGroups.py -i ofavClones.ibsMat -g ofavClones.geno.gz -b ofavCloneBams -r "OGM024.1 OGM024.2" -o ofav \
      --meta metadata.csv --meta-sample sample --meta-groups bank depthZone --quality ofavQuality5X.txt
      tips coloured by bank and shaped by depthZone; coverage reported for borderline joins
  cloneGroups.py -i ofav.ibsMat -g ofav.geno.gz -b ofavBams --metric geno --conf 0.99 --check 0.0676 -o ofav
      a clone-free set (the genotype step's output): tree with the clone cut, and any pair below it

distance and cut (--metric, --cut):
  auto    genotype distance at the first --conf level that passes both checks
          (--min-well of pairs with at least --min-shared shared confident sites; replicate
          tightness, largest replicate distance / median distance, at most --max-tight), else IBS
  geno    genotype distance at the first --conf level that passes, else the last level, with a warning
  ibs     IBS only
  both    as geno, and the IBS tree and grouping are written too, with the pairs they disagree on
  cut auto: for genotype distance, the widest gap between joins whose lower edge is within --window x
  the replicates' largest distance; for IBS, treeTrimmer (the replicates' deepest join x --buffer).
  gap and replicates force one of those; a number sets the cut directly.

replicates:
  -r names technical replicates as <sample>.<n>; a set is every name with the same <sample>.
  Names with no row are left out with a warning; a set needs two or more to count.

exit status:
  0 done; 3 no replicate set has two or more samples and the cut is not a number;
  4 the groups would drop more than --maxdrop of the samples (reason in <out>CloneStop.txt); 1 otherwise.
```

## admixResiduals.py

```text
usage: admixResiduals.py [-h] -s SP -k MAXK -b BAMS [--cut CUT] [--meta META]
                         [--meta-sample META_SAMPLE]
                         [--meta-groups [META_GROUPS ...]]

evalAdmix leaves positive residual correlations among individuals that share something the admixture model at that K does not capture. This puts each individual in its main cluster at --maxK (the finest split NGSadmix found) and averages the residuals within and between those blocks for every K; a K that lumps two real lineages leaves positive blocks within them. The suggested K is the smallest whose worst cluster-block mean |residual| is under --cut. Sampling groups from the metadata are summarised too, as an independent check.

options:
  -h, --help            show this help message and exit
  -s, --sp SP           prefix of the NGSadmix/evalAdmix files
                        (<sp>_k<K>_best.*)
  -k, --maxK MAXK       largest K; its best .qopt sets the blocks
  -b, --bams BAMS       bam list, in the order the beagle had the samples
  --cut CUT             largest worst-block mean |residual| that counts as
                        fitted (default 0.01)
  --meta META           sample metadata CSV, for the sampling-group check
  --meta-sample META_SAMPLE
                        metadata column with the sample names (default sample)
  --meta-groups [META_GROUPS ...]
                        metadata columns whose combination is a sampling group

examples:
  admixResiduals.py -s ofav -k 7 -b ../angsd/ofavBams
      reads ofav_k1_best.corres .. ofav_k7_best.corres and ofav_k7_best.qopt here; writes
      ofav_resid.txt, ofav_residBlocks.txt and ofavResiduals.pdf
  admixResiduals.py -s ofav -k 7 -b ../angsd/ofavBams --cut 0.02 \
      --meta ../metadata.csv --meta-sample sample --meta-groups bank depthZone
      a looser cut-off, and sampling groups (bank x depthZone) as the second check

outputs:
  <sp>_resid.txt       K, worstCluster, worstClusterBlock, worstGroup, residK ('<-' marks the suggested K)
  <sp>_residBlocks.txt K, a, b, within, mean, n: every cluster block with >= 3 samples on each side
  <sp>Residuals.pdf    one residual heatmap per K, then the block means for all K
```

## cloneTestRun.py

```text
usage: cloneTestRun.py [-h] -s SP -b BAMS [--depth-global DEPTH_GLOBAL]
                       [--depth-mult DEPTH_MULT] [--maxdepth-from POS.GZ]
                       [--maxdepth MAXDEPTH]
                       [-v {base,hetbias,maxdepth,post2,mind50,maf01,dep3,dep5} [{base,hetbias,maxdepth,post2,mind50,maf01,dep3,dep5} ...]]
                       [-o OUTDIR] [--min-ind MIN_IND] [--args ARGS]
                       [-p PARTITION] [-c CPUS] [--mem MEM] [--dry-run]

Submits one ANGSD job per variant of the clone step, each in its own folder under --outdir, with the IBS matrix and genotype probabilities that cloneTestScore.py compares. Run it from the species' angsd folder. Nothing in that folder is changed.

options:
  -h, --help            show this help message and exit
  -s, --sp SP           species prefix
  -b, --bams BAMS       bam list for the clone step (e.g.
                        clones/<sp>CloneBams); paths as ANGSD reads them from
                        here
  --depth-global DEPTH_GLOBAL
                        ANGSD .depthGlobal from the depth pass (default
                        <sp>Depth.depthGlobal here)
  --depth-mult DEPTH_MULT
                        maxdepth variant: limit = this x the median total
                        depth (default 2)
  --maxdepth-from POS.GZ
                        maxdepth variant: take the median total depth from
                        this ANGSD .pos.gz (e.g. cloneTests/base/<sp>.pos.gz
                        once base has run)
  --maxdepth MAXDEPTH   maxdepth variant: use this limit directly
  -v, --variants {base,hetbias,maxdepth,post2,mind50,maf01,dep3,dep5} [{base,hetbias,maxdepth,post2,mind50,maf01,dep3,dep5} ...]
                        variants to run (default all but maxdepth, which needs
                        a depth source: see --maxdepth-from)
  -o, --outdir OUTDIR   folder for the runs (default cloneTests)
  --min-ind MIN_IND     -minInd as a proportion of samples, for every variant
                        but mind50 (default 0.75)
  --args ARGS           the base ANGSD filters (default: comRAD's [clones]
                        defaults)
  -p, --partition PARTITION
                        SLURM partition (default shortq7)
  -c, --cpus CPUS       ANGSD threads and SLURM CPUs per job (default 8)
  --mem MEM             memory per job (default 100G)
  --dry-run             print the commands; submit nothing

variants:
  base      the current [clones] options (control)
  hetbias   + -hetbias_pval 1e-5 (drop sites with skewed read ratios at heterozygous calls)
  maxdepth  + -setMaxDepth at --depth-mult x the median total depth (collapsed duplicates collect about twice the reads); the median from --maxdepth-from, else <sp>Depth.depthGlobal; or --maxdepth
  post2     -doPost 2 (flat genotype prior) instead of 1 (allele-frequency prior)
  mind50    -minInd at 50%% of samples instead of --min-ind
  maf01     -minMaf 0.01 instead of 0.05
  dep3      + -setMinDepthInd 3 (a sample's data at a site counts only with 3 or more reads)
  dep5      + -setMinDepthInd 5

examples:
  cd $koko/2bRAD/nwgom/past/angsd
  cloneTestRun.py -s past -b clonesIBSnofilter/pastCloneBams --dry-run
      print the commands without submitting
  cloneTestRun.py -s past -b clonesIBSnofilter/pastCloneBams
      submit every variant but maxdepth to shortq7 (8 CPUs, 100G each); outputs in cloneTests/<variant>/
  cloneTestRun.py -s ofav -b clones/ofavCloneBams -v base post2 -p longq7
      two variants only, on another partition
  cloneTestRun.py -s past -b clonesIBSnofilter/pastCloneBams -v maxdepth --maxdepth-from cloneTests/base/past.pos.gz
      when the depth pass's histogram is capped (deeply sequenced samples): the limit from base's site depths
  cloneTestScore.py -s past -r "PGM016.1 PGM016.2 PGM016.3 PGM076.1 PGM076.2 PGM076.3 PGM119.1 PGM119.2 PGM119.3"
      once they finish: score them
```

## cloneTestScore.py

```text
usage: cloneTestScore.py [-h] -s SP -r TECHREPS [-d DIR]
                         [--variants [VARIANTS ...]] [-b BAMS]
                         [--conf CONF [CONF ...]] [--min-shared MIN_SHARED]
                         [-o OUT] [--ratio RATIO] [--min-width MIN_WIDTH]
                         [--max-mult MAX_MULT] [--buffer BUFFER] [--truth TSV]
                         [--plot PDF]

For each variant folder made by cloneTestRun.py, builds two distances between samples: ANGSD's IBS (1 - IBS) and genotype discordance (the share of sites where two samples' genotype calls differ, among sites both called with probability above --conf), and scores how well each separates clones. tight = largest distance within a replicate set / median distance between samples (lower is better); clarity = widest empty stretch in the average-linkage tree above the replicates / median distance (higher is better). The cut is at the bottom of that stretch.

options:
  -h, --help            show this help message and exit
  -s, --sp SP           species prefix used by cloneTestRun.py
  -r, --techreps TECHREPS
                        technical replicates, space-separated <sample>.<n>
  -d, --dir DIR         folder of variant runs (default cloneTests)
  --variants [VARIANTS ...]
                        variants to score (default every folder in --dir with
                        outputs)
  -b, --bams BAMS       bam list the runs used (default: the -b from the first
                        variant's command.txt)
  --conf CONF [CONF ...]
                        genotype call confidence levels (default 0.90 0.95
                        0.99)
  --min-shared MIN_SHARED
                        shared confident sites a pair needs to count as well
                        measured (default 1000)
  -o, --out OUT         also write the table here (tab-separated)
  --ratio RATIO         first-gap cut: the next join must be at least this
                        many times higher (default 1.5)
  --min-width MIN_WIDTH
                        first-gap cut: and the gap at least this share of the
                        median distance (default 0.02)
  --max-mult MAX_MULT   first-gap cut: look no higher than this many times the
                        replicates' largest distance (default 3)
  --buffer BUFFER       treeTrimmer cut (cut4): the replicates' largest
                        distance x this (default 1.10)
  --truth TSV           simulated truth (simCloneKin.py): sample1, sample2,
                        relation for every related pair; adds, per cut, the
                        true clone pairs grouped and the related or unrelated
                        pairs wrongly grouped
  --plot PDF            also draw every tree (one page per variant and
                        distance) into this PDF

examples:
  cd $koko/2bRAD/nwgom/past/angsd
  cloneTestScore.py -s past -r "PGM016.1 PGM016.2 PGM016.3 PGM076.1 PGM076.2 PGM076.3 PGM119.1 PGM119.2 PGM119.3"
      every variant in cloneTests/, IBS and genotypes at call confidence 0.90, 0.95 and 0.99
  cloneTestScore.py -s ofav -r "$(awk -F' *= *' '$1=="TECHREPS" {print $2}' ../popgen.config)" --conf 0.95
      replicates taken from the config; one confidence level
  cloneTestScore.py -s past -r "..." -d cloneTests --variants base post2 -o pastCloneTests.tsv
      two variants, and the table also written to a file
  cloneTestScore.py -s past -r "..." --conf 0.99 --plot pastCloneTests.pdf
      the table, plus one tree per variant and distance with its cut and the replicates coloured

cuts:
  cut  (red in the plots): the bottom of the widest empty stretch above the replicates.
  cut2 (blue): walking up from the replicates' largest distance, the bottom of the first gap
  where the next join is at least --ratio times higher and the gap at least --min-width of the
  median distance; none if no such gap within --max-mult x the replicates' distance.
  cut3 (green): the widest gap whose lower edge is within --max-mult x the replicates' largest
  distance (the widest gap near the replicates; gaps between relatives start above that).
  cut4: treeTrimmer, the replicates' largest distance x --buffer (comRAD's cut for IBS).

reading it:
  compare rows for the same species. A good distance has low tight, high clarity,
  every replicate set in one group, and (for genotypes) most pairs with at least
  --min-shared shared sites. Clone counts that change between rows show what each
  choice would do to the data.
```

## simCloneKin.py

```text
usage: simCloneKin.py [-h] -o OUTDIR [--variant VARIANT] [-s SP]
                      [--sites SITES] [--fst FST] [--cov COV]
                      [--cov-sd COV_SD] [--weak WEAK] [--err ERR]
                      [--somatic SOMATIC] [--big-clone BIG_CLONE]
                      [--fixed-het FIXED_HET] [--scale SCALE] [--seed SEED]

Simulates samples from two lineages, each with unrelated founders, full-sib families whose parents are also sampled (parent-offspring pairs), half-sibs, selfed sibs, clones (extra ramets of a genet) and technical replicates (the same colony sequenced again). Reads are drawn per sample and site with the sample's own coverage and a sequencing error rate; genotype probabilities use a population-frequency prior (as ANGSD -doPost 1) and IBS uses one random read per sample and site (as -doIBS 1). Writes <out>/<variant>/<sp>.geno.gz, .ibsMat and command.txt, the bam list, and a truth table, so cloneTestScore.py can score it with --truth.

options:
  -h, --help            show this help message and exit
  -o, --outdir OUTDIR   folder for the simulation (created)
  --variant VARIANT     name of this run under <outdir> (default base)
  -s, --sp SP           file prefix (default sim)
  --sites SITES         SNPs (default 8000)
  --fst FST             differentiation between the two lineages (default 0.1)
  --cov COV             median per-sample coverage (default 7)
  --cov-sd COV_SD       spread of coverage between samples, log scale (default
                        0.5)
  --weak WEAK           extra weak samples per lineage, at about 2x, all
                        clones or replicates of others (default 6)
  --err ERR             sequencing error per base (default 0.005)
  --somatic SOMATIC     somatic mutation rate between clone ramets, per site
                        (default 2e-4)
  --big-clone BIG_CLONE
                        lineage B: one genet with this many ramets (default 0)
  --fixed-het FIXED_HET
                        with --big-clone: share of sites forced heterozygous
                        in that genet (default 0)
  --scale SCALE         multiply the design counts (default 1)
  --seed SEED           random seed (default 1)

examples:
  simCloneKin.py -o simTest
      the default design (about 125 samples, 8000 SNPs) in simTest/cloneTests/base; then:
  cd simTest && cloneTestScore.py -s sim -r "$(cat techreps.txt)" --truth truth.tsv --conf 0.95 0.99 --plot sim.pdf
  simCloneKin.py -o simPast --big-clone 25 --fixed-het 0.3
      lineage B gets one genet with 25 ramets and extra fixed heterozygosity (past's comb)
  simCloneKin.py -o simLow --cov 2.5
      everything at about 2.5x: where genotype calls run short and IBS should do better
  simCloneKin.py -o simTest --variant noisy --err 0.02
      a second run of the same design (same samples) with more sequencing error, scored next to base;
      a different design (counts, --big-clone, --seed) needs its own --outdir

design (per lineage; --scale multiplies the counts):
  12 unrelated founders; 4 full-sib families (2 parents sampled + 3 offspring);
  2 half-sib groups (1 parent, 2 mates, 2 offspring each); 1 selfed family (parent + 3 selfed offspring);
  4 genets with 1-3 extra clone ramets; 2 colonies sequenced 3 times as technical replicates
```

## plotStep.py

```text
usage: plotStep.py [-h] [-s SP] [--dir DIR] [--meta META]
                   [--meta-sample META_SAMPLE]
                   [--meta-groups [META_GROUPS ...]] [--color COLOR]
                   [--shape SHAPE] [--order [ORDER ...]] [--levels LEVELS]
                   [--colors COLORS] [--shapes SHAPES] [--kcolor KCOLOR]
                   [--prefix PREFIX] [-k K] [--kin KIN] [--qdepth QDEPTH]
                   [--pipeline PIPELINE]
                   {reads,pools,symbionts,pca,admix,kin,outliers,neutral,fst,het,theta}

Draws comRAD's quick-look plots for one step from the files that step wrote, as one PDF per step (several panels or pages where there is more than one thing to show), plus PNGs of each page in run/<pipeline>/plots/ for comrad report. Meant for a first look; final figures belong in R.

positional arguments:
  {reads,pools,symbionts,pca,admix,kin,outliers,neutral,fst,het,theta}
                        which step's plots

options:
  -h, --help            show this help message and exit
  -s, --sp SP           species prefix (SP)
  --dir DIR             the pipeline's folder (default here)
  --meta META           sample metadata CSV
  --meta-sample META_SAMPLE
                        metadata column with the sample names (default sample)
  --meta-groups [META_GROUPS ...]
                        metadata columns: sampling groups (colour, then shape)
  --color COLOR         factor (metadata column, or lineage) that colours
                        points (default: first --meta-groups)
  --shape SHAPE         factor that sets point shapes (default: second --meta-
                        groups)
  --order [ORDER ...]   factors that group and order samples along bar plots
                        (default: --meta-groups)
  --levels LEVELS       level order per factor: 'bank: WFGB EFGB Stetson;
                        depthZone: Shallow Meso'
  --colors COLORS       colours per level: 'bank: WFGB=1b9e77 EFGB=d95f02;
                        lineage: Mc1=red' (hex without #, or names)
  --shapes SHAPES       point shapes per level: 'depthZone: Shallow=o Meso=^'
                        (matplotlib markers)
  --kcolor KCOLOR       admixture bars: this colour for the first cluster, the
                        rest lighter toward white in MAXK steps (hex without
                        #, or a name)
  --prefix PREFIX       admix: draw the run whose files start <SP><prefix>
                        (Neutral: the kin-free, neutral K run, written to
                        ngsAdmix/<SP>NeutralAdmix.pdf)
  -k K                  admix: K that orders the samples (default 2); neutral:
                        the chosen K
  --kin KIN             kin: the KIN cut (default 0.125)
  --qdepth QDEPTH       het: QDEPTH, for the Quality<N>X.txt coverage file
                        (default 5)
  --pipeline PIPELINE   pipeline name, for run/<pipeline>/plots (default
                        popgen)

kinds:
  reads      filteredReads/<SP>Reads.pdf: reads per sample; host (and symbiont) alignment rate per sample
  pools      RawCounts.pdf (reads pipeline, project folder): raw and kept reads per pool
  symbionts  bams/<SP>Symbionts.pdf: share of symbiont reads per genus, per sample
  pca        pcangsd/<SP>Pca.pdf: PC1-2 and PC3-4 of the PCAngsd covariance, % variance on the axes
  admix      ngsAdmix/<SP>Admix.pdf: NGSadmix bars for K = 2..MAXK (best run, one sample order, colours matched across K), and PCAngsd's own admixture at its K; evalAdmix residual RMSD and log-likelihood against K
  kin        relatedness/<SP>Kin.pdf: pairwise kinship (theta) from NGSremix, the KIN cut marked
  outliers   angsd/lineages/<SP>Outliers.pdf: PCAngsd selection statistic along the genome per lineage (outliers in red); within-lineage PCAs
  neutral    ngsAdmix/<SP>Neutral.pdf: neutral, kin-free bars next to the full run's at K; neutral PCA
  fst        sfs/<SP>Fst.pdf: weighted FST between lineages (neutral and all sites) and between groups within lineages
  het        het/<SP>Het.pdf: heterozygosity per sample by lineage and group; against coverage
  theta      theta/<SP>Theta.pdf: pi, Watterson's theta and Tajima's D per lineage (and the THETA_N subsample)

examples (from the species folder):
  plotStep.py pca -s sint --meta meta.csv --meta-groups bank depthZone
      PCA coloured by bank and shaped by depthZone
  plotStep.py pca -s mcav --meta meta.csv --meta-groups bank depthZone --color source
      coloured by another metadata column (PCA_COLOR in popgen.config), e.g. published vs new reads
  plotStep.py admix -s ofav --meta meta.csv --meta-groups bank depthZone -k 2
      bars for every K, samples ordered by group, then by their main cluster at K = 2
  plotStep.py kin -s sint --kin 0.125
  plotStep.py pools --dir /path/to/project
      the reads pipeline's per-pool yield

sample order and colours:
  Samples are grouped by the --order factors (default --meta-groups, in that order), with a divider and label
  per group, and ordered within each group by main cluster at -k (default 2), then by its share. Colours of
  clusters are matched across K (each K's clusters to the previous K's, by overlap), so a lineage keeps its colour.
  A plot that already exists is kept: a rerun writes <name>_v2.pdf, then _v3, and so on (the report shows the newest).
  --kcolor colours the admixture clusters from one colour toward white in MAXK steps (PLOT_KCOLOR). --levels, --colors and --shapes set level order, colours and shapes for named factors; levels not listed
  come after the listed ones, with default colours and shapes. In popgen.config these are PLOT_COLOR,
  PLOT_SHAPE, PLOT_ORDER, PLOT_LEVELS, PLOT_COLORS and PLOT_SHAPES.

  plotStep.py pca -s ofav --meta meta.csv --meta-groups bank depthZone \
      --levels "bank: WFGB EFGB Bright McGrail" --colors "depthZone: Shallow=gold Meso=1f78b4" --color depthZone --shape bank
```

## kmeansBic.py

```text
usage: kmeansBic.py [-h] -c COV -k MAXK [--var VAR] [--pcs PCS]
                    [--starts STARTS] [--seed SEED] -o OUT

Clusters samples by k-means on their principal components and scores each K by BIC, as adegenet's find.clusters does before DAPC (Jombart et al. 2010): BIC = n ln(WSS / n) + K ln(n), with WSS the within-cluster sum of squares. The suggested K is the lowest BIC (adegenet's 'min' criterion). The PCs come from a PCAngsd covariance matrix, so genotype uncertainty is accounted for as in the PCA. An independent read on K: like PCAngsd's own count, it responds to any discrete grouping, families and batches included.

options:
  -h, --help       show this help message and exit
  -c, --cov COV    PCAngsd covariance matrix (.cov)
  -k, --maxK MAXK  largest K to test (from 1)
  --var VAR        keep the PCs that explain this share of the variance
                   (default 0.8)
  --pcs PCS        keep this many PCs instead
  --starts STARTS  k-means random starts per K; the best is kept (default 50)
  --seed SEED      random seed (default 1)
  -o, --out OUT    output table

examples:
  kmeansBic.py -c ../pcangsd/ofavPcangsd.cov -k 7 -o ofav_bic.txt
      K = 1..7 on the PCs that explain 80% of the variance (at most n/3 of them)
  kmeansBic.py -c ../pcangsd/ofavNeutralPcangsd.cov -k 7 -o ofavNeutral_bic.txt --var 0.9
      the kin-free, neutral run, on PCs explaining 90%
  kmeansBic.py -c ../pcangsd/sintPcangsd.cov -k 10 --pcs 20 -o sint_bic.txt
      a fixed number of PCs

output:
  K, BIC, WSS, and '<-' on the lowest BIC; a comment line gives the PCs used.

tested on simulated genotypes with noise (50 per population, 3000 SNPs): one population -> 1; two at FST 0.1
or 0.02 -> 2; three at FST 0.05 -> 3; three at FST 0.01 -> 1 (too weak). Clone groups and close families
also form clusters, so run it on clone-free data and read it beside the kin-free rerun.
```
