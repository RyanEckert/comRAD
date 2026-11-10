# Tutorial: from raw 2bRAD reads to lineages, FST and demography

This walks through one project from start to finish: raw pool fastqs from one
sequencing run holding two coral species, through lineages, diversity and FST, to
moments demographic models. Every command is run on the cluster's login node;
`comrad` submits its own jobs.

The example files are in the repository's
[`examples/`](https://github.com/RyanEckert/comRAD/tree/main/examples) folder: `cluster.config`,
`reads.config`, `popgen.config`, `demography.config`, `sampleID.csv`, `metadata.csv`. The
[guide](guide/setup.md) explains each stage in more depth; this page is the path through it.

## 1. What you need

**Software** (`comrad check <pipeline>` lists anything missing):

- SLURM (job arrays; nothing else needed for submitting)
- ANGSD (with realSFS, thetaStat), bowtie2, samtools
- NGSadmix, evalAdmix, NGSremix, ngsParalog (PCAngsd and cutadapt install with comRAD)
- Python 3.9+ (comRAD installs the Python packages it needs; R is not needed)

The examples keep everything under one storage folder exported as `$koko`; any
exported variable works in a config path, or none.

**Install comRAD** (straight from GitHub; pip puts the commands in `~/.local/bin`):

```bash
python3 -m pip install --user git+https://github.com/RyanEckert/comRAD.git
comrad -h
```

Load the Python module you will run with before installing: pip writes that
interpreter into every script it installs.

**Install the other tools:** `comrad setup` checks for each one and, for any
that is missing, looks for an environment module, then uses conda/mamba (with no
conda on the cluster, it installs Miniforge into the folder first; `--no-conda`
skips that), then downloads or builds it from source (pip for the Python tools). It installs into the folder you
give it, never your home folder, and prints the manual steps for anything it could
not install (`docs/install.md` has the same list).

```bash
comrad setup --prefix $koko/bin/comradDeps
echo "source $koko/bin/comradDeps/comrad.env" >> ~/.bashrc    # PATH for what it installed
```

A module it finds is reported, not loaded: add the `module load` lines it prints
to `~/.bashrc` (they are also in `comrad.env`).

**Cluster settings**, once per project: partitions, the largest node's memory, email
and any extra sbatch flags your cluster needs (an account, a constraint).

```bash
cd $koko/2bRAD/myProject
comrad init cluster && nano cluster.config      # see examples/cluster.config
```

Every pipeline run in that folder, or up to three folders below it, reads it.
On a busy partition a driver with no time limit can wait behind other jobs;
`driver_time = 2-00:00:00` (or `comrad run popgen -t 2-00:00:00`) can let it
start sooner. If it runs out, `comrad run popgen` picks up at the step it was in.

**Genomes:** indexed once, beforehand.

```bash
bowtie2-build sintGenome.fa sintGenome.fa && samtools faidx sintGenome.fa
```

## 2. Inputs

| file | what it holds |
|---|---|
| `rawReads/<pool>.fastq.gz` | one fastq per pool (i7/i5 index), named by pool; trim2bRAD.py names each library `<pool>_<barcode>` from it |
| `sampleID.csv` | `pool_barcode,sample` per line, no header. Technical replicates are `<sample>.1`, `.2`, `.3`; the original library is always `.1` |
| `metadata.csv` | one row per sample, replicates included; a `sample` column plus whatever you group by (`bank`, `depthZone`, ...) |

The sample-name prefix decides the species (`SGM` = sint, `PGM` = past, set in
`SPECIES`). Every sample that ends up with reads must have a metadata row:
`comrad check popgen` stops and lists any that don't.

```
myProject/                    cluster.config  reads.config  sampleID.csv  metadata.csv  rawReads/
  sint/                       popgen.config  demography.config
    demuxReads/ filteredReads/ bams/ angsd/ pcangsd/ ngsAdmix/ relatedness/ sfs/ het/ theta/ demo/
    run/popgen/               one .log and .done per step, popgenSummary.txt (the diary), versions.txt, settings copies
  past/                       ...
```

## 3. Reads: trim, demultiplex, rename, split by species

```bash
cd $koko/2bRAD/myProject
comrad init reads           # or: cp $koko/bin/comRAD/examples/reads.config .
nano reads.config
comrad check reads
comrad run reads
comrad status reads         # last events, and whether the driver is still queued
```

Steps: `download` (skipped without `BS_PROJECTS`), `demux`, `demuxcheck`, `rename`,
`species`. When it finishes, each species folder has `demuxReads/*.demux.fastq.gz`,
one per sample, and `rawReads/RawCounts.pdf` plots the yield per pool. Check:

```bash
cat rawReads/RawCounts                            # raw and kept reads per pool
ls sint/demuxReads/*.demux.fastq.gz | wc -l       # samples per species
ls unassigned/demuxReads/                         # libraries with no sample name: empty wells, barcode errors
```

A barcode with fewer than 100,000 reads (trim2bRAD.py's default) gets no file, so
failed libraries are simply absent. If a sample you expect is missing, find its
barcode in `rawReads/logs/<pool>.demux.log`.

**File names along the way:** `demuxReads/<s>.demux.fastq.gz` (demultiplexed,
deduplicated) → `filteredReads/<s>.filtered.fastq.gz` (cutadapt) →
`filteredReads/<s>.host.bam`, `<s>.symbiont.bam` (bowtie2; removed by `clean`) →
`bams/<s>.bam`, `bams/<s>.sym.bam` (after the host/symbiont split). Logs are in
each folder's `logs/`; job lists are `<SP><Step>Jobs`.

## 4. Popgen: mapping to lineages, diversity and FST

```bash
cd sint
comrad init popgen          # or copy examples/popgen.config and edit it
nano popgen.config
comrad check popgen
comrad run popgen
```

The steps, in order (`comrad list popgen` shows which are done):

| stage | steps |
|---|---|
| reads and mapping | cutadapt, readcounts, map (host and symbiont), alignrates, split, clean |
| QC and clones | depth, clones, **clonegroups** (pause 1), filterbams |
| genotypes and structure | paralog, genotype, structure, admix, symbionts, **K** (pause 2) |
| lineages | assign, kin, lineagebams, lineagesnps, outliers, neutral |
| statistics | sfssites, saf, fst, fstall, groups, het, theta, tarball |

### Leaving analyses out

`SKIP` lists steps not to run, e.g. `SKIP = het theta groups`. `check` says what
each one changes, as a note or a WARNING, and stops on a step later ones need (only
`het theta groups fstall symbionts neutral kin outliers fst` can be skipped):

```
note: SKIP het: no per-individual heterozygosity
WARNING: SKIP kin: related individuals stay in: SFS, theta, Tajima's D and FST will include close kin, ...
WARNING: SKIP outliers: SNPs under selection stay in the SFS and FST: demography assumes neutral sites, ...
WARNING: SKIP fst: no lineage FST, fstall is skipped as well, and the demography pipeline cannot run ...
```

The same lines go in the diary when the run starts, a skipped step's files are left
out of the results tarball, and the demography `check` repeats the kin and outlier
warnings. To run a step you skipped earlier, take it out of `SKIP` and use
`comrad from popgen <step>`.

### Pause 1: clones (with `CLONE_PAUSE = yes`)

`CLONE_PAUSE = yes` is the default. Samples with almost no data look close to
everything and pull the clone groups together, so samples below `CLONE_MINQ`
(default 0.20) of sites at ≥ `QDEPTH`x (default 5) are left out of the clone step;
the diary names them and `angsd/clones/sintCloneLowCoverage.txt` lists them. They
would fall below `MINQ` (0.25) afterwards in any case. Check that no technical
replicate is among them: a replicate left out cannot help set the cut.

With `CLONE_METRIC = auto`, the clone step measures **genotype distance** (the share
of sites where two samples' confident genotype calls differ) when coverage supports
it, and **IBS** otherwise; the diary and `sintCloneMetric.txt` say which and why.
The cut (`CLONE_CUT = auto`) comes from the technical replicates: for genotype
distance, the bottom of the widest gap just above them; for IBS, their deepest join
× `TREE_BUFFER` (the treeTrimmer method). [Clones and K](guide/decisions.md) explains
both.

The run stops after drawing `angsd/sintCloneDend.pdf`: the tree, technical
replicates in bold, and the cut as a red line. Checks that come with it:
`sintCloneReplicates.txt` gives each replicate's nearest sample, and a replicate
whose nearest sample is outside its own clone group gives a WARNING (a mislabelled
or failed library); `sintCloneBorderline.txt` lists joins near the cut with each
side's coverage (a weak library on one side points to a clone, good coverage on both
to a close relative); and if the clone groups would drop more than `CLONE_MAXDROP`
(0.5) of the samples, the run stops: that is a cut or data problem, not a result.

Look for a gap between the clone-level pairs (replicates and true clones, near the
bottom) and everything else. If the red line sits in that gap, accept it; if not,
set your own.

```bash
nano popgen.config     # accept:  CLONE_PAUSE = no      or set your own:  CLONE_CUT = 0.08
comrad run popgen
```

If you already trust the clone call (for example on a rerun), set `CLONE_PAUSE = no` before
running and it never stops here.

`filterbams` then keeps the best-covered sample of each clone group and drops
samples below `MINQ`. `genotype` warns about samples covering under half the median
proportion of SNPs; check where they fall in the PCA.

### Pause 2: K (with `K` empty)

Look at `pcangsd/sintPca.pdf` (the PCA, coloured by your sampling groups) and
`ngsAdmix/sintAdmix.pdf` (admixture bars for every K, one sample order). The diary, and
`run/popgen/K.log`, give five reads on K:

```
note  K: evalAdmix elbow K = 2 (residual RMSD against K)
note  K: PCAngsd K = 2 (eigenvectors kept by its MAP test + 1; ...)
note  K: residual-block K = 2 (smallest K whose worst cluster-block mean |residual| < 0.01; plots: ngsAdmix/sintResiduals.pdf)
note  K: k-means BIC K = 2 (lowest BIC of k-means on the PCA scores, as adegenet find.clusters; ...)
note  K: runs in the dominant mode per K: 1:20/20 2:20/20 3:18/20 4:9/20
```

- **evalAdmix elbow:** where the residual curve flattens. A geometric heuristic: it
  always returns something, and is never K = 1 or MAXK.
- **PCAngsd:** the number of eigenvectors its MAP test kept, plus one; its own admixture at
  that K is the bottom row of the admixture plot.
- **k-means BIC:** k-means clustering of the PCA scores, the K with the lowest BIC (DAPC's
  `find.clusters`). Independent of NGSadmix; like PCAngsd, it also responds to families and
  batches.
- **Residual-block K:** read from the evalAdmix residual heatmaps in
  `ngsAdmix/sintResiduals.pdf`. Each sample is placed in its main cluster at MAXK;
  a K that lumps two real lineages leaves red blocks within them and blue between.
  It's the smallest K with no block mean above `RESID_MAX`. It cannot tell you K
  is too high, and the cut-off is a judgement call: look at the plots.
- **Mode support:** how many of the 20 NGSadmix runs agree at each K. Low support
  means that K is not a stable solution.

Then:

```bash
nano popgen.config     # K = 2
comrad run popgen
```

If you already know K, set it before the run and it never pauses here.

`K = 1` (no structure) is allowed: every sample is one lineage, `fst` and `fstall`
are skipped (FST among sampling groups still runs in `groups`), and the demography
pipeline's `check` stops, since its models need two lineages. At K ≥ 2, fewer than
two lineages with pure individuals still stops at `assign`.

### When it finishes

Each step has drawn its quick-look plots along the way (reads, symbionts, PCA, admixture,
kinship, outliers, the neutral rerun, F<sub>ST</sub>, heterozygosity, θ), and
`comrad report popgen` gathers them, with the warnings and key tables, into
`run/popgen/popgenReport.html`. [Results, plots and the report](guide/results.md) describes
each, and how to set their colours and order.

`sintResults.tar.gz` holds everything the R side needs:

| folder | contents |
|---|---|
| `.` | bam lists, lineage calls (`BamsClusters`), kin table, alignment and retention rates |
| `qc` | depth distributions, per-sample coverage at 3/5/10x |
| `clones` | clone groups, distance used and cut, dendrogram(s), replicate check, borderline joins, samples left out of the clone step |
| `structure` | beagle, PCAngsd covariance and F_IS, best NGSadmix Q per K, evalAdmix fit and residuals, neutral rerun and its comparison with the lineages |
| `plots` | every step's quick-look PDF |
| `relatedness` | NGSremix table |
| `diversity` | heterozygosity (with `HetDepth.txt`), π/θ/Tajima's D (and the `THETA_N` subsample lists), FST between lineages and groups, 2D spectra |
| `run` | the diary, `popgenReport.html`, `versions.txt` and every copy of the settings |

**Heterozygosity and coverage.** Heterozygosity is estimated per sample on the SFS
sites with at least 10 reads in that sample. A sample with less coverage reaches
10x at fewer, and different, sites, so `het` also tests whether heterozygosity
follows coverage (Spearman's ρ against `nSites` and against the proportion of
sites ≥ `QDEPTH`x, permutation p). The result is `het/<SP>HetDepth.txt` and a diary note;
p < 0.05 gives a WARNING. If it warns, compare groups with coverage as a
covariate, or downsample the bams to one depth and rerun `het`.

**θ and sample size.** `theta/<SP>Theta.out` gives π, Watterson's θ, Tajima's D and
the number of sites per lineage, from all its pure, kin-free individuals (`n`). π is
the comparison to make between lineages of different size; θ_W and Tajima's D depend
on n. With `THETA_N` set (a number, or `min` for the smallest lineage) the same
statistics are also estimated from that many individuals per lineage, in the `Sub`
columns; the individuals used are listed in `theta/sub/` and are the same on every
run. Compare θ_W and Tajima's D in those columns when lineage sizes differ much.

**Sampling groups.** `groups` keeps lineage × sampling-group combinations with at
least `MIN_GROUP` (default 5) individuals for group FST; the diary names each one left out
and its n. Heterozygosity per group is the per-sample `het/<SP>.het` joined with
`angsd/<SP>Groups`, which the R side does.

## 5. Demography: moments models between lineages

Run in the same species folder, after popgen has finished:

```bash
comrad init demography      # or copy examples/demography.config
nano demography.config
comrad check demography
comrad run demography
```

For each pair of lineages it builds 100 contig-bootstrap spectra, fits every model
in `MODELS` 50 times (`FIT_SEEDS` × `FIT_RUNS`), keeps the best run of each, and
ranks them by CLAIC: AIC corrected for linked sites using the bootstrap spectra.

Two checks run automatically:

- **convergence:** do the best runs of a model agree?
- **nesting:** does every larger model fit at least as well as the smaller model it
  contains? If not, its optimiser stopped early (`demo/nesting.txt`).

The result is `sintDemography.tar.gz`. The model names follow the SI / IM / AM / SC
convention with `2m` (heterogeneous migration) and `2N` (heterogeneous Ne)
suffixes; `momentsModels.py -h` describes each.

## 6. Rerunning

| you want to | command |
|---|---|
| carry on after a stop or pause | `comrad run popgen` |
| redo from a step because its inputs or settings changed | `comrad from popgen <step>` |
| start at a step when earlier ones were done by hand | `comrad from popgen <step> --mark-earlier` |
| run one step now, in this shell | `comrad step popgen <step>` |
| see which steps are done | `comrad list popgen` |

`run` skips steps already done, and inside a step skips jobs that already
finished. It does not notice files you changed by hand, or a setting used by a
step that is already done: that's what `from` is for. Every job in a rewound step
reruns, even where an old output file exists.

Moving a species folder while its driver is running breaks the run (paths are
absolute); move it between runs.

## 7. What gets recorded

The diary (`run/popgen/popgenSummary.txt`) has one line per event and a `note` with the
numbers for your methods at every step; `comrad summary popgen` prints them step by step.
`run/popgen/` also keeps the tool versions and a copy of the settings each time they change.
[Running and rerunning](guide/running.md) has the details, along with how jobs, logs and memory
are handled.

## 8. When something stops

```bash
comrad status popgen              # the reason is the last line
comrad summary popgen             # every step's notes and warnings, and any open pause or stop
less run/popgen/<step>.log        # the step's full output
```

comrad-jobs has already retried failed jobs (on `long_queue` after a time limit, with more
memory after running out); a stop means that wasn't enough. [Troubleshooting](troubleshooting.md)
lists the usual messages and what to do.
