# Settings

Every setting, its default and what it does, for each pipeline. This page is written from the
pipeline templates by `tools/makeDocs.py`; `comrad defaults <pipeline>` prints the same thing in
the terminal, and `comrad init <pipeline>` writes it as a config to fill in.

A config is `KEY = value` lines, then one `[section]` per job list. In a section, `cpus` and `mem`
are per job (per array task), and `cpus` is also the tool's own thread count. `args` lists
**changes** to the tool's default options: an option you give replaces the default's value or is
added, and `-option none` removes it. `queue` (`long` = the cluster's `long_queue`) and `slurm`
(extra sbatch flags as `key=value`) can go in any section. Settings are read from the template's
defaults, then `cluster.config`, then the pipeline's own config; later lines win.

## reads

Run from the project folder. **Required** settings must be filled in; `comrad check reads` says which are missing.

| setting | default | what it does |
|---|---|---|
| `RAW` **required** | `rawReads` | folder for the raw .fastq.gz (downloads land here) |
| `BS_PROJECTS` | *(empty)* | BaseSpace projects to download into RAW, space-separated; empty = reads already in RAW |
| `SAMPLE_SHEET` **required** | *(empty)* | CSV for sampleRename.py: old name (pool_barcode), new sample name |
| `RENAME_ARGS` | *(empty)* | extra sampleRename.py options, e.g. --header |
| `SPECIES` **required** | *(empty)* | species=prefix pairs, e.g. mcav=MGM ofav=OGM; samples starting with the prefix go to &lt;species&gt;/demuxReads |
| `MIN_KEPT` | `25` | stop if any pool keeps less than this percent of its reads (kept = after dedup) |

### Job sections

| section | runs | cpus | mem | default args |
|---|---|---|---|---|
| `[download]` | fetchReads.py, one job | 1 | 200GB |  |
| `[demux]` | trim2bRAD.py: find the tag, deduplicate, demultiplex; one job per pool | 2 | 16GB | `--gzip` |

## popgen

Run from the species folder. **Required** settings must be filled in; `comrad check popgen` says which are missing.

| setting | default | what it does |
|---|---|---|
| `SP` **required** | *(empty)* | prefix for every output file, e.g. sint |
| `GENOME_FASTA` **required** | *(empty)* | host reference, already bowtie2-indexed and samtools-faidx'd |
| `SYMGENOME` | *(empty)* | symbiont metagenome, already indexed; empty = host only |
| `SYMGENUS` | *(empty)* | with SYMGENOME: chromosome&lt;TAB&gt;genus table |
| `META` **required** | *(empty)* | sample metadata CSV, one row per sample |
| `META_SAMPLE` | `sample` | column with the sample names (bam names without .bam) |
| `META_GROUPS` | *(empty)* | columns whose combination is a sampling group, e.g. site depth; empty = no group FST |
| `EXTRA_READS` | *(empty)* | folder of extra &lt;sample&gt;.filtered.fastq.gz, already quality-filtered (e.g. published reads) |
| `MIN_IND` | `0.75` | proportion of samples an ANGSD site must have data in, for every analysis run ([depth] has its own) |
| `SKIP` | *(empty)* | steps to leave out: het theta groups fstall symbionts neutral kin outliers fst (check says what each changes) |

### Clones

| setting | default | what it does |
|---|---|---|
| `TECHREPS` | *(empty)* | technical replicates as &lt;sample&gt;.1 &lt;sample&gt;.2 ...; empty = set CLONE_CUT to a number |
| `CLONE_METRIC` | `auto` | clone distance: auto (genotype distance when coverage supports it, else IBS), geno, ibs or both (geno, with IBS compared) |
| `CLONE_CUT` | `auto` | auto (geno: widest gap just above the replicates; ibs: replicates' deepest join x TREE_BUFFER), gap, replicates, or a number |
| `TREE_BUFFER` | `1.10` | IBS cut: deepest join within any replicate set x this |
| `CLONE_PAUSE` | `yes` | pause after the dendrogram; to go on set it to no (or set CLONE_CUT to a number) and run again |
| `QDEPTH` | `5` | depth for the sample-quality filters below (proportion of sites at &gt;= QDEPTH x, from the depth step) |
| `CLONE_MINQ` | `0.20` | leave samples below this proportion out of the clone step (near-empty samples look close to everything; their genotype distances are unmeasurable) |
| `CLONE_MAXDROP` | `0.5` | stop if clone grouping would drop more than this proportion of the samples (a threshold or data problem) |
| `MINQ` | `0.25` | after clones: drop samples below this proportion |

### Structure and lineages

| setting | default | what it does |
|---|---|---|
| `MAXK` | `10` | NGSadmix K = 1..MAXK |
| `ADMIX_REPS` | `20` | NGSadmix replicates per K |
| `K` | *(empty)* | empty: pause at the K step to choose; a number: use that K |
| `RESID_MAX` | `0.01` | residual-block K: smallest K whose worst cluster-block mean \|evalAdmix residual\| is below this (a judgement call) |
| `ASSIGN` | `0.90` | ancestry needed to assign an individual to a lineage |
| `MINPURE` | `2` | fewest pure individuals a lineage may have |
| `LINEAGE_MAXHET` | *(empty)* | within each lineage, remove sites heterozygous in more than this share of samples (e.g. 0.5: likely paralogs); empty = report only |
| `KIN` | `0.125` | NGSremix kinship at or above which one of a pair is dropped |
| `MIN_GROUP` | `5` | smallest lineage x sampling group kept for group FST (smaller ones are named in the diary) |

### Diversity

| setting | default | what it does |
|---|---|---|
| `THETA_N` | *(empty)* | theta always from all pure individuals per lineage; a number or 'min' (smallest lineage) also gives theta from that many per lineage (Sub columns) |
| `PLOT_COLOR` | *(empty)* | factor that colours points; empty = the first META_GROUPS column |
| `PLOT_SHAPE` | *(empty)* | factor that sets point shapes; empty = the second META_GROUPS column |
| `PLOT_ORDER` | *(empty)* | factors that group and order samples along bar plots; empty = META_GROUPS |
| `PLOT_LEVELS` | *(empty)* | level order, e.g.  bank: WFGB EFGB Bright; depthZone: Shallow Meso; lineage: Mc1 Mc2 |
| `PLOT_COLORS` | *(empty)* | colours, e.g.  bank: WFGB=1b9e77 EFGB=d95f02; lineage: Mc1=red  (hex without the #, or names) |
| `PLOT_SHAPES` | *(empty)* | point shapes, e.g.  depthZone: Shallow=o Meso=^  (matplotlib markers: o ^ s D v P X *) |
| `PLOT_KCOLOR` | *(empty)* | admixture bar colours: the first cluster in this colour, the rest lighter toward white in MAXK steps, e.g. 1B9E77; empty = default palette |

### Job sections

| section | runs | cpus | mem | default args |
|---|---|---|---|---|
| `[cutadapt]` | quality trim, one job per sample | 1 | 4GB | `-q 15,15 -m 36` |
| `[map]` | bowtie2, one job per sample (host and symbiont) | 1 | *(see note)* | `--no-unal --score-min L,16,1 --local -L 16` |
| `[split]` | host/symbiont split, scans and sorting | 1 | 8GB |  |
| `[depth]` | ANGSD depth pass and plotQC.py, one job each<br>`min_ind = 0.5` | 1 | 200GB |  |
| `[clones]` | ANGSD IBS and genotype calls on every sample above CLONE_MINQ, one job; no per-sample depth floor (low-coverage samples still count, as in Matz's 2bRAD protocol) | 1 | 200GB | `-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -skipTriallelic 1 -snp_pval 1e-5 -minMaf 0.05` |
| `[paralog]` | ngsParalog, one job per chunk of sites<br>`min_ind = 0.5` | 1 | 8GB |  |
| `[genotype]` | ANGSD on the kept samples, one job | 1 | 200GB | `-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -snp_pval 1e-5 -minMaf 0.05` |
| `[structure]` | PCAngsd (PCA and Fis) | 8 | 64GB |  |
| `[admix]` | NGSadmix, one job per K and replicate | 10 | 8GB |  |
| `[evaladmix]` | evalAdmix, one job per K | 4 | 16GB |  |
| `[kin]` | NGSremix, one job | 8 | 100GB |  |
| `[lineagesnps]` | ANGSD within each lineage | 1 | 32GB | `-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -snp_pval 1e-5 -minMaf 0.05` |
| `[selection]` | PCAngsd selection scan within each lineage | 8 | 64GB |  |
| `[outliers]` | PCAngsd within each lineage, outliers removed | 8 | 64GB |  |
| `[neutral]` | PCAngsd on neutral, kin-free SNPs (NGSadmix and evalAdmix there use [admix] and [evaladmix]) | 8 | 64GB |  |
| `[sfssites]` | ANGSD site discovery for the SFS, per lineage | 1 | 64GB | `-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -skipTriallelic 1 -sb_pval 1e-3 -hetbias_pval 1e-3` |
| `[saf]` | ANGSD site allele frequencies (SFS, theta, FST): per lineage, full panel and sampling group | 1 | 32GB | `-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30` |
| `[fst]` | realSFS 2D spectra and FST, per pair | 20 | 64GB |  |
| `[het]` | ANGSD + realSFS per sample on the SFS sites, for heterozygosity; 10x floor: one sample alone has no pooling to correct low-depth bias | 1 | 8GB | `-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -setMinDepthInd 10` |
| `[theta]` | realSFS and thetaStat per lineage | 4 | 32GB |  |

`[map]` with `mem` empty asks for the bowtie2 index size × 1.5 + 2 GB per task.

## demography

Run from the species folder, after popgen. **Required** settings must be filled in; `comrad check demography` says which are missing.

| setting | default | what it does |
|---|---|---|
| `SP` **required** | *(empty)* | same prefix as popgen.config |
| `MU` **required** | `2e-8` | mutation rate per site per generation, for scaled units |
| `MODELS` | `SI SI2N IM IM2N IM2m IM2m2N IMa IMa2N AM AM2N AM2m AM2m2N SC SC2N SC2m SC2m2N IMghost SCghost` |  |
| `BOOT_JOBS` | `25` | bootstrap jobs per pair, each with its own seed |
| `BOOT_REPS` | `4` | replicates per job (total = BOOT_JOBS x BOOT_REPS) |
| `FIT_SEEDS` | `5` | moments jobs per model, each with its own seed |
| `FIT_RUNS` | `10` | optimisations per job (total per model = FIT_SEEDS x FIT_RUNS) |
| `MAXITER` | `1000` | optimiser iterations per run |
| `PROJECT` | *(empty)* | down-project every fit to N1 N2 haploid samples; empty = none |

### Job sections

| section | runs | cpus | mem | default args |
|---|---|---|---|---|
| `[boot]` | realSFS contig bootstraps, BOOT_JOBS per pair | 20 | 8GB |  |
| `[moments]` | momentsModels.py, one job per model and seed; long runs<br>`queue = long` | 4 | 32GB |  |

## cluster.config

One per project, written by `comrad init cluster`. comrad reads it from the folder a pipeline runs in
or up to three folders above, so one in the project folder covers every species folder.

| setting | default | what it does |
|---|---|---|
| `queue` | *(empty)* | partition for jobs; empty = the cluster's default partition |
| `long_queue` | *(empty)* | partition for jobs that hit a time limit, and for sections with queue = long; empty = none (they stop) |
| `max_mem` | `250GB` | out-of-memory retries double a job's memory up to this; set it to your largest node's memory |
| `max_running` | `64` | array tasks running at once, per job list |
| `email` | *(empty)* | mail on failed jobs and when the pipeline ends; empty = no mail |
| `slurm` | *(empty)* | extra sbatch flags for every job, as key=value or bare flags: e.g. account=mylab constraint=amd exclusive |
| `driver_queue` | *(empty)* | partition for the driver job, which runs for the whole pipeline; empty = long_queue, else queue |
| `driver_time` | *(empty)* | driver time limit, e.g. 2-00:00:00 (helps it start sooner on a busy partition); empty = the partition maximum |
| `driver_cpus` | `8` | driver CPUs, used by the steps' own short work (sorting, R plots) |
| `driver_mem` | `32GB` | driver memory |
