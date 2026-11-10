# Results, plots and the report

When `popgen` finishes, the species folder holds `<SP>Results.tar.gz`: everything an
analysis in R or Python needs, without the bams and intermediate files.

| folder | holds |
|---|---|
| `.` | bam lists, lineage calls (`<SP>BamsClusters`), sampling groups, kin table, alignment and retention rates, read counts |
| `qc` | depth distributions, each sample's proportion of sites at 3, 5, 10 and `QDEPTH`x |
| `clones` | clone groups, distance used and cut, dendrogram(s), replicate check, borderline joins, samples left out of the clone step |
| `structure` | beagle, PCAngsd covariance and inbreeding, best NGSadmix Q per K, evalAdmix fit and residuals, neutral rerun and its comparison with the lineages, outlier scans |
| `plots` | quick-look PDFs from each step (reads, symbionts, PCA, admixture, kinship, outliers, neutral run, F<sub>ST</sub>, heterozygosity, θ); also inline in the report |
| `relatedness` | NGSremix table |
| `diversity` | heterozygosity, π/θ/Tajima's D, F<sub>ST</sub> between lineages and among groups, 2D spectra |
| `run` | the diary, the report, tool versions and every copy of the settings |

## Plots


Each step draws quick-look plots of what it made (`plotStep.py`), one PDF per step. They are
for a first look; final figures belong in R. A plot that fails gives a WARNING and never stops
the run.

| step | file | shows |
|---|---|---|
| `alignrates` | `filteredReads/<SP>Reads.pdf` | reads per sample; host and symbiont alignment rates |
| `symbionts` | `bams/<SP>Symbionts.pdf` | share of symbiont reads per genus, per sample |
| `structure` | `pcangsd/<SP>Pca.pdf` | PCA, PC1–2 and PC3–4, % variance on the axes |
| `admix` | `ngsAdmix/<SP>Admix.pdf` | admixture bars for K = 2…`MAXK`, one sample order, cluster colours matched across K; evalAdmix fit and likelihood against K |
| `kin` | `relatedness/<SP>Kin.pdf` | pairwise kinship, the `KIN` cut marked |
| `outliers` | `angsd/lineages/<SP>Outliers.pdf` | selection statistic along the genome per lineage, outliers in red; within-lineage PCAs |
| `neutral` | `ngsAdmix/<SP>Neutral.pdf` | neutral, kin-free bars beside the full run's at K; neutral PCA coloured by lineage |
| `fst`, `fstall`, `groups` | `sfs/<SP>Fst.pdf` | F<sub>ST</sub> heatmaps on one colour scale: between lineages (neutral and all sites), between groups within lineages |
| `het` | `het/<SP>Het.pdf` | heterozygosity per sample by lineage and by group; against coverage |
| `theta` | `theta/<SP>Theta.pdf` | π, θ<sub>W</sub> and Tajima's D per lineage (and the `THETA_N` subsample) |
| reads pipeline | `rawReads/RawCounts.pdf` | raw and kept reads per pool |

The clone dendrograms (`angsd/clones/<SP>CloneDend.pdf`, and `angsd/<SP>CloneFreeDend.pdf` from
`genotype`) and the depth plots (`angsd/<SP>Depth.pdf`) come from their own steps.

A rerun never overwrites a plot: it writes `<name>_v2.pdf`, then `_v3`, and so on. The report
shows the newest; the tarball's `plots/` folder has them all.

### Colours, shapes and order

Samples are grouped along bar plots by `META_GROUPS` (or `PLOT_ORDER`); points are coloured by
the first group column and shaped by the second (or `PLOT_COLOR`, `PLOT_SHAPE`). `PLOT_LEVELS`,
`PLOT_COLORS` and `PLOT_SHAPES` set level order, colours and shapes per factor; a factor is a
metadata column named in the config, or `lineage`:

```
PLOT_COLOR = depthZone
PLOT_SHAPE = bank
PLOT_LEVELS = bank: WFGB EFGB Bright McGrail; depthZone: Shallow Meso; lineage: Of1 Of2
PLOT_COLORS = depthZone: Shallow=f4a582 Meso=0571b0; lineage: Of1=red Of2=1b9e77
PLOT_SHAPES = bank: WFGB=o EFGB=^ Bright=s McGrail=D
PLOT_KCOLOR = 7570B3
```

`PLOT_KCOLOR` colours the admixture clusters from one colour toward white in `MAXK` steps,
so each cluster keeps its colour in every K panel.

Colours are hex without the `#` (the config treats `#` as a comment) or names. `comrad check`
stops on a factor that is not a metadata column and warns about levels it cannot find.

## The report

`run/popgen/popgenReport.html` (also in the tarball) is one page per run:

- **Read first:** every warning, and any pause or stop still open.
- **Samples:** lineage sizes, samples left out of the clone step, the replicate check, clone
  groups, kin dropped, coverage per sample.
- **Structure:** the K evidence, NGSadmix fit, the neutral comparison, inbreeding.
- **Differentiation and diversity:** F<sub>ST</sub>, θ, heterozygosity.
- **Plots:** every step's quick-look plots, inline.
- **Settings, tool versions** and the full diary.

`comrad report popgen` writes it at any point, including on a paused or stopped run.

## Diversity

**Heterozygosity** is estimated per sample, on the SFS sites with at least 10 reads in that
sample (`het/<SP>.het`; means per lineage in `<SP>Ho.txt`). A sample with less coverage
reaches 10 reads at fewer, and different, sites, so `het` also tests whether heterozygosity
follows coverage: Spearman's ρ against the number of sites and against the proportion of
sites at `QDEPTH`x, with a permutation p (`het/<SP>HetDepth.txt`). A p below 0.05 is a
WARNING. If you see it, compare groups with coverage as a covariate, or downsample the bams to
one depth and rerun `het`.

**θ.** `theta/<SP>Theta.out` gives π, Watterson's θ, Tajima's D and the number of sites per
lineage, from all its pure, kin-free samples (`n`). π compares well between lineages of
different size; θ<sub>W</sub> and Tajima's D depend on n. With `THETA_N` set (a number, or
`min` for the smallest lineage's n), the same statistics are also estimated from that many
samples per lineage, in the `Sub` columns. The samples used are listed in `theta/sub/` and
are the same on every run. Compare θ<sub>W</sub> and Tajima's D in those columns when lineage
sizes differ much.

**Sampling groups.** `groups` estimates F<sub>ST</sub> among sampling groups (the combination
of `META_GROUPS` columns) within each lineage, for groups with at least `MIN_GROUP` (5) samples;
the diary names each group left out and its n. Heterozygosity per group is `het/<SP>.het`
joined with `angsd/<SP>Groups`.

## F<sub>ST</sub>

`sfs/<SP>LineageFst.out` has unweighted and weighted F<sub>ST</sub> for each pair of lineages, on
sites outside the outlier tags (realSFS `-whichFst 1`: Hudson's estimator as in Bhatia et al. 2013). `sfs/<SP>LineageFstAll.out` has the
same on every site, outliers in, so you can see how much the outliers add.
