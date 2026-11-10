# Clones in a highly heterozygous coral

Four coral species from the northwestern Gulf of Mexico (*Montastraea cavernosa*, *Orbicella
faveolata*, *Stephanocoenia intersepta* and *Porites astreoides*), 2bRAD, mapped to each
species' genome. Each had technical replicates: the same colony, sequenced as two or three
libraries.

## What went wrong with IBS

The first clone step used IBS: ANGSD's `-doIBS 1`, which takes one random read per sample at each
site, and the treeTrimmer cut (the replicates' deepest join × `TREE_BUFFER`). It worked for three
species. For *P. astreoides* (128 samples) it did not:

- only 4 of 9 replicates had one of their own replicates as their nearest sample;
- replicates of one colony sat further apart than some unrelated-looking pairs, so no cut could
  group the replicates without grouping other samples too;
- the dendrogram had no clean gap between clone-level and other pairs.

**Why.** With one random read per site, two libraries of the same colony draw different alleles
at about half of its heterozygous sites. So the IBS distance between clones is not near zero: it
grows with the colony's heterozygosity. In a species with high heterozygosity, clones sit as
far apart as close relatives. Filters did not fix it: removing sites with skewed read ratios at
heterozygous calls (`-hetbias_pval`), depth floors, a maximum depth, a flat genotype prior and a
lower minor-allele cut-off each left the IBS tree as confused, and the read-ratio filter alone
removed 32–58% of sites.

**Was the heterozygosity real?** Collapsed paralogs (two gene copies mapped as one) also look
heterozygous. Here, sites heterozygous in most samples had only 1.14–1.20 times the usual depth
(paralogs would show about twice), balanced reference and alternative read counts, and were spread
across all 14 chromosomes of the reference (jaPorAstr1). The heterozygosity was real.

## Genotype distance

Genotype distance asks a different question: at sites where both samples' genotype calls are
confident (probability above 0.99), how often do the calls differ? Clones share genotypes, so
their distance sits near zero however heterozygous they are; relatives differ at a share of
sites set by how related they are. It needs enough coverage for confident calls, which these
data had.

With genotype distance at 0.99 and the cut at the widest gap just above the replicates:

| species | clone groups | samples dropped as clones | pairs with ≥ 1000 shared sites | replicate tightness |
|---|---|---|---|---|
| *P. astreoides* | 19 | 82 | 100% | 0.064 |
| *O. faveolata* | 13 | 18 | 100% | 0.086 |
| *S. intersepta* | 9 | 16 | 100% | 0.122 |
| *M. cavernosa* | 4 | 5 | 100% | 0.080 |

Replicate tightness is the replicates' largest distance divided by the median distance between
all samples: small means replicates sit far below everything else. The groups barely moved across the
filter variants: only a depth floor of 5 reads per sample changed one *P. astreoides* call
(18 groups, 83 dropped). For the three species where IBS worked, genotype distance agreed
with it.

Two details from this data set:

- **A weak library joins its clone a little higher.** In *O. faveolata*, OGM009 joined two clone
  mates at 0.068, above most clone pairs. It had 0.56 of its sites at ≥ 5x against 0.87 and
  0.89 for the other two, and shared about half as many confident sites with them. A weak
  library, not a relative. This is why `<SP>CloneBorderline.txt` lists each side's coverage.
- **The IBS cut needs its buffer.** For *O. faveolata*, IBS with treeTrimmer grouped correctly
  only with `TREE_BUFFER` between 1.05 and 1.15 (1.00 split a replicate pair, 1.20 merged
  another), hence the default of 1.10.

## What comRAD does now

- `CLONE_METRIC = auto` (the default) uses genotype distance when the data support it: at least
  90% of pairs share 1000 or more confident sites, and replicate tightness is at most 0.15. It
  tries 0.99, then 0.95, and otherwise falls back to IBS.
- `CLONE_CUT = auto` uses the widest gap above the replicates for genotype distance, and
  treeTrimmer for IBS.
- The clone ANGSD run already writes genotype probabilities, so this costs no extra job.
- `<SP>CloneMetric.txt` and the diary record the choice and the numbers behind it.

[Clone calls on a simulated pedigree](clone-simulation.md) tests both choices against a known
truth.
