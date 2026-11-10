# Popgen

The `popgen` pipeline takes one species from reads to lineages, kinship, diversity and
F<sub>ST</sub>. Run it from the species folder.

```bash
cd sint
comrad init popgen
nano popgen.config             # SP, GENOME_FASTA, META, TECHREPS, ... (the required lines are marked)
comrad check popgen
comrad run popgen
```

`check` confirms every program, file and metadata row is there, and prints the full options
each step will run with. `comrad list popgen` shows which steps are done.

## Steps

| stage | step | does |
|---|---|---|
| **Reads and mapping** | `cutadapt` | quality-trims each sample |
| | `readcounts` | counts reads per sample |
| | `map` | bowtie2 to the host and, if given, the symbiont metagenome |
| | `alignrates` | alignment-rate tables |
| | `split` | keeps each read in host or symbiont, dropping reads that map to both; sorts and indexes `bams/<s>.bam` |
| | `clean` | removes the raw mapping files; writes the bam lists |
| **QC and clones** | `depth` | ANGSD depth pass and `plotQC.py`: each sample's proportion of sites at 3, 5, 10 and `QDEPTH`x |
| | `clones` | ANGSD IBS and genotype calls on all samples, near-empty ones left out |
| | `clonegroups` | clone distance (genotype or IBS) and cut from technical replicates, clone groups, dendrogram. **Pause 1** |
| | `filterbams` | keeps the best-covered sample of each clone group; drops samples below `MINQ` |
| **Genotypes and structure** | `paralog` | ngsParalog: removes sites that look like collapsed paralogs |
| | `genotype` | genotype likelihoods at SNPs in all kept samples; coverage at those SNPs per sample; clone-free check |
| | `structure` | one SNP per tag; PCAngsd PCA and per-sample inbreeding |
| | `admix` | NGSadmix at K = 1 to `MAXK`, `ADMIX_REPS` runs each; evalAdmix fit and residuals |
| | `symbionts` | reads per symbiont genus (host-only runs skip it) |
| | `K` | lays out the evidence for K. **Pause 2** |
| **Lineages** | `assign` | assigns each sample to a lineage if its ancestry is at least `ASSIGN`, otherwise *Admixed* |
| | `kin` | NGSremix; drops one of each pair related at or above `KIN` |
| | `lineagebams` | pure, kin-free samples per lineage |
| | `lineagesnps` | SNPs within each lineage, with heterozygosity per site (`LINEAGE_MAXHET`); PCAngsd selection scan |
| | `outliers` | outlier SNPs (FDR < 0.05 on any PC, per lineage); within-lineage PCA without them |
| | `neutral` | reruns structure without outliers or kin and compares the lineage calls |
| **Statistics** | `sfssites` | sites that pass the filters in every lineage, outlier tags removed |
| | `saf` | site allele frequencies per lineage |
| | `fst` | 2D spectra and F<sub>ST</sub> between lineages, on those neutral sites |
| | `fstall` | the same on the full panel, outliers in, for comparison |
| | `groups` | F<sub>ST</sub> among sampling groups within each lineage |
| | `het` | heterozygosity per sample |
| | `theta` | π, Watterson's θ and Tajima's D per lineage |
| | `tarball` | the results tarball and the HTML report |

The two pauses are explained in [Clones and K](decisions.md), and the outputs in [Results](results.md).

## Filters

`MIN_IND` (0.75) is the proportion of samples a site must have data in, for every ANGSD
analysis; a section's `min_ind` overrides it. The depth pass uses 0.5.

No step has a per-sample depth floor except heterozygosity (10 reads). The others pool
genotype likelihoods across samples, which handles low depth; one sample on its own cannot.
To add a floor for a sensitivity run, it is one line in that step's section:

```ini
[genotype]
args = -setMinDepthInd 3
```

Each section's `args` lists **changes** to the defaults, which `init` writes above it as a
comment. `comrad defaults popgen` prints them all, and the
[settings reference](../reference/settings.md#popgen) lists every step's defaults.

## Leaving analyses out

`SKIP` lists steps not to run, e.g. `SKIP = het theta groups`. Only
`het theta groups fstall symbionts neutral kin outliers fst` can be skipped. `check` says what
each one changes, as a note or a WARNING, and stops on any step a later one needs:

```text
note: SKIP het: no per-individual heterozygosity
WARNING: SKIP kin: related individuals stay in: SFS, theta, Tajima's D and FST will include close kin, ...
WARNING: SKIP outliers: SNPs under selection stay in the SFS and FST: demography assumes neutral sites, ...
```

A skipped step's files stay out of the results tarball. To run it later, take it out of `SKIP`
and run `comrad from popgen <step>`.

## File names along the way

`demuxReads/<s>.demux.fastq.gz` → `filteredReads/<s>.filtered.fastq.gz` (cutadapt) →
`filteredReads/<s>.host.bam` and `<s>.symbiont.bam` (bowtie2, removed by `clean`) →
`bams/<s>.bam` and `bams/<s>.sym.bam`. Analysis files start with `SP` (`sintClones.ibsMat`),
logs are in each folder's `logs/`, and job lists are `<SP><Step>Jobs`.

## Plots

Each step draws quick-look plots of what it made; [Results, plots and the report](results.md#plots)
lists them and how to set their colours and order.
