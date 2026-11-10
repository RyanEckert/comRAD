# Clones and K

`popgen` pauses twice, where the right answer depends on your data: the clone threshold and
the number of genetic clusters, K. In both cases comRAD lays out the evidence, and you answer
in the config, so the decision is recorded with the run. If you already know the answer, put
it in the config before running and the pause is skipped.

## Pause 1: clones

Clonal organisms put several samples of one genotype in a collection. Keeping them all
inflates relatedness and skews allele frequencies, so the `clonegroups` step finds clone
groups and `filterbams` keeps the best-covered sample of each.

**Which samples go in.** A sample with almost no data looks close to everything in an IBS
tree and pulls unrelated samples together. Samples covering less than `CLONE_MINQ` (0.20) of
sites at `QDEPTH`x (5x) are left out of the clone step. They would fall below `MINQ` (0.25)
afterwards in any case. The diary names them, and `angsd/clones/<SP>CloneLowCoverage.txt`
lists them. Check that none of your technical replicates are among them.

**The distance.** Two ways to measure how different two samples are:

- **Genotype distance**: the share of sites where the two samples' genotype calls differ, among
  sites both are called at with high confidence (probability above 0.99, or 0.95). Clones and
  replicates sit near 0, whatever the species' heterozygosity.
- **IBS** (1 - IBS): one random read per sample and site. It works at very low coverage, but two
  libraries of one heterozygous genotype differ at about half its heterozygous sites, so clones
  of a highly heterozygous species sit far from 0 and close to their relatives.

With `CLONE_METRIC = auto` (the default) comRAD uses genotype distance when the data support
it: at least 90% of pairs share 1000 or more confident sites, and the replicates are tight
(their largest distance at most 0.15 of the median distance). It tries 0.99, then 0.95, and
otherwise uses IBS. The diary and `<SP>CloneMetric.txt` give the numbers behind the choice.
`geno` or `ibs` forces one; `both` uses genotypes and also writes the IBS tree, with the pairs
the two group differently.

In simulations with clones, replicates, parents and offspring, full, half and selfed sibs at
about 7x, genotype distance found every clone with no relative merged, including a large clone
of a highly heterozygous genet; IBS merged parent-offspring and selfed-sib pairs and missed
that clone. At about 2.5x there are too few confident calls, `auto` falls back to IBS, and IBS
still finds the clones but can merge close relatives: check the borderline joins.

**The cut.** Technical replicates are known clones, so they show how different two libraries of
one genotype can look. `CLONE_CUT = auto` places the cut:

- for genotype distance, at the bottom of the widest gap between joins that starts within three
  times the replicates' largest distance (clones sit just above the replicates; relatives are
  much further up);
- for IBS, at the deepest join within any replicate set times `TREE_BUFFER` (1.10), the method
  of [treeTrimmer](https://github.com/RyanEckert/treeTrimmer).

`gap` or `replicates` forces one rule; a number sets the cut yourself (in the units of the
distance used, given in `<SP>CloneMetric.txt`).

The run pauses after drawing `angsd/<SP>CloneDend.pdf`: the tree, technical replicates in bold,
and the cut as a red line (with genotype distance, `<SP>CloneDendIBS.pdf` shows the IBS tree for
comparison). Look for a gap between the clone-level pairs (replicates and true clones, near the
bottom) and everything else. If the red line sits in that gap, accept it. If not, set your own.

```bash
nano popgen.config      # accept:  CLONE_PAUSE = no       or your own:  CLONE_CUT = 0.08
comrad run popgen
```

Checks that come with it:

- **Replicate check.** `<SP>CloneReplicates.txt` gives each replicate's nearest sample and the
  nearest sample outside its set. A replicate whose nearest sample is outside its own clone
  group gives a WARNING: a mislabelled or failed library.
- **Borderline joins.** `<SP>CloneBorderline.txt` lists joins grouped as clones but more than
  1.5 times the replicates' distance, and the nearest join left ungrouped, with each side's
  coverage. A weak library on one side points to a clone; good coverage on both, to a close
  relative.
- **Too many clones.** If the clone groups would drop more than `CLONE_MAXDROP` (half) of the
  samples, the run stops. That is a cut or data problem, not a result.
- **Clone-free check.** After genotyping, the `genotype` step draws `angsd/<SP>CloneFreeDend.pdf`
  for the kept samples, with the same distance and cut, and warns about any pair still below it.

Replicates with no reads, or left out for low coverage, are warned about and left out; a
replicate set needs two samples to count. If no set has two, set `CLONE_CUT` to a number.

## Pause 2: K

`admix` runs NGSadmix at every K from 1 to `MAXK`, `ADMIX_REPS` times each, and evalAdmix on
the best run of each K. The `K` step then gives five reads on K, in the diary and
in `run/popgen/K.log`:

```text
note  K: evalAdmix elbow K = 2 (residual RMSD against K)
note  K: PCAngsd K = 2 (eigenvectors kept by its MAP test + 1; ...)
note  K: residual-block K = 2 (smallest K whose worst cluster-block mean |residual| < 0.01; plots: ngsAdmix/sintResiduals.pdf)
note  K: k-means BIC K = 2 (lowest BIC of k-means on the PCA scores, as adegenet find.clusters; table in ngsAdmix/sint_bic.txt)
note  K: runs in the dominant mode per K: 1:20/20 2:20/20 3:18/20 4:9/20
```

| read | what it is |
|---|---|
| evalAdmix elbow | where the residual curve flattens. A geometric rule: it always returns something, and never K = 1 or `MAXK` |
| PCAngsd | the number of eigenvectors PCAngsd's MAP test kept, plus one (Meisner & Albrechtsen 2018); its own admixture at that K is the bottom row of `ngsAdmix/<SP>Admix.pdf` |
| k-means BIC | k-means clustering of the PCA scores, K with the lowest BIC (the `find.clusters` step of DAPC; Jombart et al. 2010). Independent of NGSadmix. Like PCAngsd, it also answers to families and batches |
| residual-block K | from the evalAdmix residuals in `ngsAdmix/<SP>Residuals.pdf`. Each sample is placed in its main cluster at `MAXK`; a K that lumps two real lineages leaves positive residuals within them. This is the smallest K with no block mean above `RESID_MAX`. It cannot tell you K is too high, and the cut-off is a judgement call: look at the plots |
| mode support | how many of the NGSadmix runs agree at each K. Low support means that K is not a stable solution |

```bash
nano popgen.config      # K = 2
comrad run popgen
```

`K = 1` (no structure) is allowed: every sample is one lineage, F<sub>ST</sub> between lineages
is skipped (F<sub>ST</sub> among sampling groups still runs), and demography does not apply.
At K ≥ 2, the run stops at `assign` if fewer than two lineages have pure samples, or any has
fewer than `MINPURE`.

## After K

- **Lineages.** Each sample with at least `ASSIGN` (0.90) ancestry from one cluster belongs to
  that lineage; the rest are *Admixed* and stay out of the per-lineage statistics.
- **Kin.** NGSremix estimates relatedness; for each pair at or above `KIN` (0.125), one sample is
  dropped: the one in most remaining pairs, then the one with less coverage.
- **Outliers.** A PCAngsd selection scan within each lineage; SNPs with FDR < 0.05 on any PC are
  outliers, and their tags are removed from the SFS.
- **Neutral check.** Structure is rerun without outliers or kin: PCAngsd, and NGSadmix with
  evalAdmix at every K from 1 to `MAXK`, with the same reads on K as at the pause
  (`ngsAdmix/<SP>NeutralAdmix.pdf`, `<SP>NeutralResiduals.pdf`). Close kin can pull PCs and
  clusters of their own and leave residual blocks, so this is K without them. If none of the
  elbow, PCAngsd, residual-block and k-means BIC K gives the K you chose, that is a WARNING. The lineage
  calls at your K are then compared with the full run's; a mismatch is a WARNING too, with the
  table in `ngsAdmix/<SP>NeutralCompare.txt`. Either way, change `K` and
  `comrad from popgen assign` if the kin-free run says otherwise.
