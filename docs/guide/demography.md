# Demography

The `demography` pipeline fits two-population divergence models with
[moments](https://momentsld.github.io/moments/) to the 2D spectrum of every pair of lineages.
It is an extra: run it in the species folder after `popgen` has finished, with K of at least 2.

```bash
comrad init demography
nano demography.config         # SP (as in popgen.config), MU
comrad check demography
comrad run demography
```

## Steps

| step | does |
|---|---|
| `setup` | samples per lineage (`demo/<SP>PopN.txt`) and the lineage pairs |
| `boot` | contig-bootstrap 2D spectra for each pair: `BOOT_JOBS` jobs × `BOOT_REPS` replicates (100 by default) |
| `bootjoin` | joins the bootstrap spectra per pair |
| `moments` | every model in `MODELS`, `FIT_SEEDS` jobs × `FIT_RUNS` optimisations each (50 by default) |
| `collect` | the best run of each model, the model table, and the checks below |
| `tarball` | `<SP>Demography.tar.gz` and the HTML report |

Models are ranked by CLAIC: AIC corrected for linked sites, from the bootstrap spectra, with
Godambe standard errors for every parameter. `MU` (mutation rate per site per generation)
converts the scaled parameters into individuals and generations.

## Models

The names follow the usual convention: **SI** strict isolation, **IM** isolation with
continuous migration, **AM** ancient migration (migration, then isolation), **SC** secondary
contact (isolation, then migration). The suffix **2m** adds heterogeneous migration across the
genome, and **2N** heterogeneous effective size; the **ghost** models add an unsampled
population. `momentsModels.py` (see the [scripts reference](../reference/scripts.md#momentsmodelspy))
fits them, and its source defines each one.

## Checks

- **Convergence:** do a model's best runs agree, across seeds? A model whose best fit came from
  only one seed is flagged in `demo/<SP>Convergence.txt`; give it more seeds before reading its
  parameters.
- **Nesting:** every larger model should fit at least as well as the smaller model it contains.
  If one does not, its optimiser stopped early (`demo/<SP>Nesting.txt`).
- **Bounds and Godambe warnings** from each kept run's log go in `demo/<SP>Warnings.txt`.

`popgen` decisions carry over: if kin or outliers were skipped there, the demography `check`
warns, because the models assume unrelated samples and neutral sites.
