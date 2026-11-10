# Worked examples

Two examples of the clone step, the decision that shapes every analysis after it.

- [Clones in a highly heterozygous coral](clones-past.md): real data from four Gulf of Mexico
  corals, where IBS could not separate clones from relatives in *Porites astreoides* and
  genotype distance could. How comRAD came to use genotype distance by default.
- [Clone calls on a simulated pedigree](clone-simulation.md): a simulated population with known
  clones, replicates and relatives, at good and at low coverage. What each distance and cut gets
  right and wrong, and how `CLONE_METRIC = auto` chooses between them.

To run either on your own data, `cloneTestRun.py`, `cloneTestScore.py` and `simCloneKin.py` are
installed with comRAD ([scripts](../reference/scripts.md)).
