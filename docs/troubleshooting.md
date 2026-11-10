# Troubleshooting

## Finding out why it stopped

```bash
comrad status popgen                 # the last events; the reason for a pause or stop is the last line
comrad summary popgen                # each step's notes and warnings from its last run, and any open pause or stop
less run/popgen/<step>.log           # the step's full output
ls -t <folder>/<jobname>*.e* | head  # the jobs' own error logs, next to their outputs
```

From a project folder, `comrad summary popgen --dir sp1 sp2 sp3` does several species at once.

A **pause** is a decision for you (the clone dendrogram, K); a **stop** is something a rerun on
its own will not fix. In both cases the message says what to look at. Fix it, then
`comrad run popgen` carries on from the step that stopped.

comrad-jobs has already retried failed jobs before a step stops: on `long_queue` after a time
limit, with double the memory (up to `max_mem`) after running out, on another node after a node
failure, and once more for anything else. A stop means that wasn't enough.

## Common messages

### Before the run

| message | what to do |
|---|---|
| `not on PATH: ...` | the external tools are missing: `comrad setup --prefix DIR` installs them, then `source DIR/comrad.env` |
| `numpy, scipy or matplotlib missing for ...` | comRAD was installed with a different Python: load the Python you run with and reinstall |
| `python cannot import moments and demes` | the same, for the demography pipeline |
| `... has no bowtie2 index` / `has no .fai` | index the genome: `bowtie2-build genome.fa genome.fa && samtools faidx genome.fa` |
| `<n> samples without a row ...` | add the samples to the metadata CSV (`META`), with every `META_GROUPS` column filled in |
| `CLONE_CUT = auto needs TECHREPS` | list your technical replicates, or set `CLONE_CUT` to a number |
| `PLOT_... names <x>, which is not a column of ...` | a plot setting names a factor that is not a metadata column (or `lineage`); check the spelling and case |
| `WARNING: PLOT_LEVELS: <factor> has no level <x>` | a level is spelled differently from the metadata (case counts): `Shallow` is not `shallow` |
| `<KEY> is from an older config: read as ...` | a renamed setting; it still works, but rename it as the message says |

### Reads and mapping

| message | what to do |
|---|---|
| `pools below <n>% kept` | a pool lost most of its reads to deduplication or demultiplexing; check its barcodes in `rawReads/logs/<pool>.demux.log`, or lower `MIN_KEPT` if that is expected |
| `<n> .demux.fastq.gz match no prefix in SPECIES` | sample names whose prefix is not in `SPECIES`; they stay in the project folder |
| `host bams vs samples: <a> vs <b>` | some mapping jobs did not finish; look at the `map` step's `.e` logs |
| `rows with a chromosome missing from SYMGENUS` | the symbiont genus table does not cover every sequence in `SYMGENOME` |

### Clones

| message | what to do |
|---|---|
| `clone grouping at <x> would drop <n> of <m> samples` | more than `CLONE_MAXDROP` would go as clones. Look at `<SP>CloneDend.pdf`, `<SP>CloneReplicates.txt` and `<SP>CloneBorderline.txt`. In a highly clonal species this can be real: raise `CLONE_MAXDROP`. Otherwise set `CLONE_CUT` yourself |
| `no technical replicate set has two or more bams` | replicates are missing, or below `CLONE_MINQ`; set `CLONE_CUT` to a number |
| `WARNING: ... TECHREPS are not in the clone step` | those replicates had no reads, or too little coverage; the rest still set the cut |
| `WARNING: replicates whose nearest sample is outside their clone group` | a mislabelled or failed library, or a cut that is too low; check those samples in the dendrogram |
| `genotype distance did not pass its checks ...; using IBS` | coverage is too low for confident genotype calls; IBS is used, and close relatives can be grouped with clones: check `<SP>CloneBorderline.txt`. See [Clones and K](guide/decisions.md) |
| `WARNING: <n> pairs of clone-free samples are below the clone cut` | after genotyping, some kept samples still look like clones. The sites differ from the clone run, so check those pairs before removing anything |

### Structure and lineages

| message | what to do |
|---|---|
| `a lineage has fewer than <n> pure individuals` | K is too high, or few samples reach `ASSIGN`; choose a lower K, or lower `MINPURE` on purpose |
| `fewer than 2 lineages with pure individuals at K = <k>` | the same, at K ≥ 2; consider `K = 1` |
| `WARNING: the neutral, kin-free run does not match the lineages` | structure may be driven by outlier loci or kin; compare `ngsAdmix/<SP>Neutral.pdf` with the full run |
| `<lineage> selection scan did not converge` | the PCAngsd selection scan hit its iteration cap; look at `angsd/lineages/<SP><lineage>SelPcangsd.log` |

### Spectra, FST and demography

| message | what to do |
|---|---|
| `a spectrum did not converge` | a realSFS fit hit its iteration cap; the CHECK lines in the step log say which. Raise `-maxIter` in that section's `args` |
| `THETA_N = <n> is more than the smallest lineage` | set `THETA_N = min` or a smaller number |
| `one lineage ...: the demographic models here need two lineages` | demography needs K ≥ 2 |
| `moments jobs failed` | see the `.e` logs in `demo/`; a model that cannot be fitted to your spectrum can be left out of `MODELS` |

### Warnings that do not stop the run

| warning | meaning |
|---|---|
| `samples covering under half the median proportion of SNPs` | low-coverage samples can sit apart in the PCA for technical reasons; check where they fall |
| `heterozygosity vs <x>: rho ... p < 0.05` | heterozygosity follows coverage; compare groups with coverage as a covariate, or downsample the bams and rerun `het` |
| `<kind> plots failed` | a quick-look plot could not be drawn; the step's log has the error. The results are not affected |

## Rerunning after a fix

| you want to | command |
|---|---|
| carry on after a stop or pause | `comrad run popgen` |
| redo from a step because its inputs or settings changed | `comrad from popgen <step>` |
| start at a step when earlier ones were done by hand | `comrad from popgen <step> --mark-earlier` |
| run one step now, in this shell | `comrad step popgen <step>` |

`run` does not notice a setting you changed for a step that is already done: use `from`.
Deleting an output folder by hand (for example `angsd/`) also needs a `from` at the step that
makes it; the depth step remakes the bam lists if `angsd/` is gone.

## Cluster problems

- **The driver waits in the queue for a long time.** A driver with no time limit can sit
  behind other jobs; set `driver_time` in `cluster.config` (or `comrad run popgen -t 2-00:00:00`).
  If it runs out, `comrad run popgen` picks up at the step it was in.
- **`sacct did not answer`** in the job notes: peak memory is not recorded, nothing else is
  affected.
- **bgzf errors from ANGSD or realSFS** on parallel file systems (BeeGFS), e.g.
  `Failed to read BGZF header` then `safreader.cpp ... Assertion failed` (exit 134), on files
  that pass `gzip -t`. The genotype step copies the bams to node-local disk first. Every realSFS
  and thetaStat call goes through `localSaf.py`, which fails (exit 75, so the job is retried) on
  any bgzf message instead of keeping a partial result. If the errors persist, set
  `local_copy = yes` in `cluster.config`: realSFS and thetaStat then read copies of their SAF,
  theta and FST files on node-local disk (`$TMPDIR`), which needs space for the largest SAF.
- **Moving a species folder** while its driver runs breaks the run (paths are absolute); move it
  between runs.
