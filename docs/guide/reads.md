# Reads

The `reads` pipeline turns raw pool fastqs into one file of clean, deduplicated reads per
sample, sorted into a folder per species. Run it from the project folder.

```bash
comrad init reads
nano reads.config              # SAMPLE_SHEET, SPECIES (and BS_PROJECTS to download)
comrad check reads
comrad run reads
```

## Steps

| step | does |
|---|---|
| `download` | downloads the pools from BaseSpace into `rawReads/` (skipped when `BS_PROJECTS` is empty) |
| `demux` | `trim2bRAD.py`, one job per pool: finds the 2bRAD tag, removes PCR duplicates, and writes one file per in-line barcode, `<pool>_<barcode>.demux.fastq.gz` |
| `demuxcheck` | writes `rawReads/RawCounts` (raw and kept reads per pool), and stops if a pool kept less than `MIN_KEPT` percent |
| `rename` | `sampleRename.py` renames each library from the sample sheet, after a dry run that changes nothing if anything is wrong |
| `species` | moves each species' files into `<species>/demuxReads/` |

A barcode with fewer than 100,000 reads (the `trim2bRAD.py` default) gets no file, so failed
libraries are simply absent. To find a missing sample, look up its barcode in
`rawReads/logs/<pool>.demux.log`.

!!! note "Libraries with no sample name"
    A library missing from the sample sheet keeps its `<pool>_<barcode>` name, and the
    `species` step stops if any file matches no prefix in `SPECIES`. Add a catch-all for your
    pool names, e.g. `unassigned=FG` for pools named `FG...`, to collect them in their own folder.

## Checking the result

```bash
cat rawReads/RawCounts                            # raw and kept reads per pool
ls sint/demuxReads/*.demux.fastq.gz | wc -l       # samples per species
ls unassigned/demuxReads/                         # empty wells, barcode errors, unnamed libraries
```

## Starting from your own reads

`popgen` only needs `<species>/demuxReads/<sample>.demux.fastq.gz`, one per sample. If your
reads were demultiplexed elsewhere, put them there under those names and skip `reads`.
Reads that are already quality-filtered can go in as `filteredReads/<sample>.filtered.fastq.gz`;
see [Running and rerunning](running.md#starting-part-way).
