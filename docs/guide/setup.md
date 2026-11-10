# Before you start

## What you need

- A Linux cluster with **SLURM**. comRAD submits job arrays; nothing else is needed for submitting.
- **comRAD** and the programs it runs: see [Install](../install.md). `comrad check <pipeline>`
  lists anything missing.
- **Reference genomes**, indexed once:

```bash
bowtie2-build sintGenome.fa sintGenome.fa && samtools faidx sintGenome.fa
```

For a host/symbiont system you also need a symbiont metagenome (indexed the same way) and a
two-column table of its chromosomes and their genus (`chromosome<TAB>genus`), for the
symbiont read counts. `packMetagenome.py` packs several symbiont genomes into one metagenome
with each chromosome named after its genome (`Sym001`, `Brev001`, ...), which makes that table
a one-line job.

## Inputs

| file | what it holds |
|---|---|
| `rawReads/<pool>.fastq.gz` | one fastq per pool (one Illumina index), named by pool. Downloaded for you from BaseSpace if you give `BS_PROJECTS` |
| `sampleID.csv` | one line per library, no header: `<pool>_<barcode>,<sample>` |
| `metadata.csv` | one row per sample, replicates included: a `sample` column plus whatever you group by (`bank`, `depthZone`, ...) |

Filled-in examples of each are in [`examples/`](https://github.com/RyanEckert/comRAD/tree/main/examples).

## Naming conventions

comRAD relies on a few naming rules:

- **Sample names** are the bam file names without `.bam`, and the values in the metadata's
  sample column (`META_SAMPLE`). Every sample with reads needs a metadata row;
  `comrad check popgen` lists any that don't have one.
- **Technical replicates** are `<sample>.1`, `<sample>.2`, ..., with the original library as
  `.1`. They set the clone threshold, so name them this way and list them in `TECHREPS`.
- **Species** are told apart by sample-name prefix (`SPECIES = mcav=MGM ofav=OGM`): every
  sample starting `MGM` goes to the `mcav` folder.

## Folder layout

```text
myProject/                    cluster.config  reads.config  sampleID.csv  metadata.csv  rawReads/
  sint/                       popgen.config  demography.config
    demuxReads/ filteredReads/ bams/ angsd/ pcangsd/ ngsAdmix/ relatedness/ sfs/ het/ theta/ demo/
    run/popgen/               a .log and .done per step, the diary, versions, settings copies
  past/                       ...
```

Each pipeline runs from the folder it works in: `reads` from the project folder, `popgen` and
`demography` from a species folder. Paths in a config can be absolute, relative to that
folder, or use any exported variable (`$STORAGE/genomes/sint.fa`).

## Cluster settings

Once per project, in the project folder:

```bash
comrad init cluster
nano cluster.config
```

| setting | set it to |
|---|---|
| `queue` | the partition jobs start on |
| `long_queue` | the partition for jobs that run out of time (empty: they stop instead) |
| `max_mem` | your largest node's memory; out-of-memory retries double up to this |
| `email` | where SLURM mails failures and the end of a pipeline |
| `slurm` | extra sbatch flags every job needs, e.g. `account=mylab constraint=amd` |

`sinfo -o "%P %m %l"` lists the partitions, their node memory (MB) and time limits. Every
pipeline run in this folder, or up to three folders below it, reads this file. All settings
are listed in the [settings reference](../reference/settings.md#clusterconfig).
