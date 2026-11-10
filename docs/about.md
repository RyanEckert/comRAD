# About

## Citing comRAD

If you use comRAD, cite it (the repository's `CITATION.cff` has the details; GitHub's
"Cite this repository" button formats it) and the tools it runs:

- **ANGSD**, realSFS and thetaStat: Korneliussen, Albrechtsen & Nielsen 2014, *BMC Bioinformatics*
- **PCAngsd**: Meisner & Albrechtsen 2018, *Genetics*
- **NGSadmix**: Skotte, Korneliussen & Albrechtsen 2013, *Genetics*
- **evalAdmix**: Garcia-Erill & Albrechtsen 2020, *Molecular Ecology Resources*
- **NGSremix**: Nøhr, Hanghøj, Garcia-Erill et al. 2021, *G3*
- **ngsParalog**: Linderoth, ngsParalog (GitHub)
- **bowtie2**: Langmead & Salzberg 2012, *Nature Methods*
- **samtools**: Danecek et al. 2021, *GigaScience*
- **cutadapt**: Martin 2011, *EMBnet.journal*
- **moments**: Jouganous et al. 2017, *Genetics*
- **2bRAD**: Wang, Meyer, McKay & Matz 2012, *Nature Methods*

`run/<pipeline>/versions.txt` records the version of each tool your run used.

## Where the scripts come from

Most scripts in `scripts/` are new. Some are ports or rewrites of earlier work:

| script | from | author | license |
|---|---|---|---|
| `trim2bRAD.py` | `trim2bRAD_2barcodes_dedup.pl` and `2bRAD_trim_launch_dedup.pl` in [2bRAD_denovo](https://github.com/z0on/2bRAD_denovo); fixes the Perl's misaligned quality strings (see the script's help) | Mikhail Matz | none declared |
| `bowtie2Launcher.py` | `2bRAD_bowtie2_launch.pl` in 2bRAD_denovo | Mikhail Matz | none declared |
| `ngsAdmixLauncher.py` | the NGSadmix job-list approach in 2bRAD_denovo | Mikhail Matz | none declared |
| `plotQC.py` | `plotQC.R` in [ngsTools](https://github.com/mfumagalli/ngsTools) | Matteo Fumagalli | GPL-3.0 |
| `cloneGroups.py` | the IBS cut is the method of [treeTrimmer](https://github.com/RyanEckert/treeTrimmer) | Ryan J. Eckert | |

## License

GPL-3.0-or-later, as required by the ngsTools port. The full text is in the repository's
`LICENSE`.

## Renamed from radTools

comRAD was radTools before its first release. The commands are now `comrad` (was
`radpipe.py`) and `comrad-jobs` (was `rtjobs.py`); the old names are still installed and do the
same thing. Configs, `.rtjobs/` stamps and run folders are unchanged, so a run started under
radTools carries on under comRAD. Older settings are read under their new names, and
`comrad check` says which to rename.

## Tests and development

```bash
git clone https://github.com/RyanEckert/comRAD.git && cd comRAD
python3 -m pip install --user .        # after editing a script in scripts/, install again
bash tests/run_mock.sh                 # reads -> popgen -> demography against stand-ins for every tool
bash tests/sim_clones.sh               # the clone call on simulated pedigrees, against the truth
bash tests/check_scripts.sh            # every script is listed, installed and executable
python3 tools/makeDocs.py --check      # the generated reference pages are up to date
```

`tests/run_mock.sh` checks the wiring between steps (`tests/fake.py` stands in for each
tool), including the pauses, kin drop, outliers, the host-only path and K = 1. It does not test
the tools themselves. Variants: `MOCK_HOSTONLY=1`, `MOCK_K=1`, `MOCK_METRIC=geno`,
`MOCK_OLDKEYS=1`. `tests/sim_clones.sh` simulates clones, replicates and relatives at about 7x
and 2.5x and checks which are grouped.

The documentation site is built with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/):
`pip install mkdocs-material` then `mkdocs serve` in the repository folder.
