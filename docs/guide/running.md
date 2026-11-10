# Running and rerunning

## How a run works

`comrad run <pipeline>` submits one small **driver** job that runs the pipeline's steps in
order. Each step builds a list of jobs and hands it to `comrad-jobs`, which:

- **skips** every job whose output is already there, so a rerun only does what is missing;
- **submits** the rest as one SLURM job array, at most `max_running` tasks at once;
- **reruns only what fails**: on `long_queue` after a time limit, with double the memory (up to
  `max_mem`) after running out, on another node after a node failure, and once more for
  anything else;
- **records** what each job list used, in the diary.

A step ends in one of four ways:

| outcome | means | what to do |
|---|---|---|
| done | finished | nothing |
| skip | does not apply (no symbiont genome, in `SKIP`, ...) | nothing |
| **pause** | a decision is needed: the clone threshold or K | answer it in the config, then `comrad run` again |
| **stop** | something a rerun will not fix | read the reason, fix it, then `comrad run` again |

Problems that do not need to stop the run (replicates that don't group, samples with little
coverage at the SNPs, a mismatch in the neutral check, heterozygosity that follows coverage)
are WARNING lines in the diary and the report.

## Where it is, and why it stopped

```bash
comrad status popgen                 # the last events; the reason for a pause or stop is the last line
comrad list popgen                   # every step, done or pending
less run/popgen/<step>.log           # one step's full output
ls -t <folder>/<jobname>*.e* | head  # the jobs' own error files
comrad report popgen                 # run/popgen/popgenReport.html: warnings first, then the key tables
comrad summary popgen --dir mcav ofav sint   # from the project folder: each step's notes and warnings, and what waits on you
```

| message | usually means |
|---|---|
| `missing or empty: <file>` | an earlier step's output is gone; `comrad from` that step |
| `<n> samples ... no row in META` | add the rows to the metadata |
| `a spectrum did not converge` | a realSFS fit hit its iteration cap; see the CHECK lines in the step log |
| `clone grouping would drop ...` | the clone threshold is too high, or the samples are near-empty; set `CLONE_CUT` |
| `a lineage has fewer than MINPURE pure individuals` | K is too high, or `ASSIGN` too strict |
| `not on PATH: ...` | a program is missing; see [Install](../install.md) |

## Rerunning

| you want to | run |
|---|---|
| carry on after a pause or a stop | `comrad run popgen` |
| redo from a step, because its inputs or settings changed | `comrad from popgen <step>` |
| start at a step, because the earlier ones were done another way | `comrad from popgen <step> --mark-earlier` |
| run one step now, in this shell | `comrad step popgen <step>` |
| run the whole thing in this shell instead of a driver job | `comrad run popgen --here` |

`run` skips steps that are done, and within a step, jobs that are done. It does not notice
files you changed by hand, or a new setting for a step that has already run: that is what
`from` is for. Every job in a rewound step reruns, even where an old output exists.

A driver stopped by its time limit loses nothing; `comrad run` picks up at the step it was in.
Don't move a species folder while its driver is running (paths are absolute); move it between runs.

### Starting part way

If the reads were already quality-filtered, put them in `filteredReads/<sample>.filtered.fastq.gz`
with a read-count table `filteredReads/<SP>ReadCounts` (`readCounts.py -f filtered.fastq.gz -o <SP>Read`),
then start at mapping:

```bash
comrad check popgen && comrad from popgen map --mark-earlier
```

## Jobs, logs and memory

Every job list has its own section in the config, which `init` writes with its defaults:

```ini
[map]                    # bowtie2, one job per sample (host and symbiont)
cpus = 1
mem =                    # empty = the bowtie2 index size x 1.5 + 2 GB
# default args: --no-unal --score-min L,16,1 --local -L 16
args =                   # changes only

[admix]                  # NGSadmix, one job per K and replicate
cpus = 10
mem = 8GB
```

- `cpus` is both what SLURM reserves and the tool's own thread count (`-P`, `--threads`,
  `-p`): a reservation alone does not make a program use the cores.
- `mem` is per task. A task that runs out is resubmitted with double, up to `max_mem`, so a
  low figure costs a retry, not the run. If tasks keep running out, raise it.
- `queue = long` sends that job list straight to `long_queue`, and `slurm = constraint=amd`
  adds `--constraint=amd` to its jobs. In `cluster.config`, `slurm` applies to every job.
- No time limit is asked for unless you set one, so each task gets its partition's own.

Each job list's logs are combined into one pair next to its outputs when the array ends: a
summary table with failed tasks first, then a section per task with its node, exit code and run
time. The diary notes what each job list used:

```text
note  map: sintHostMap: 128 jobs, peak 3.1 GB of 6GB, longest 0.68 h, 61.2 CPU-hours
```

## What gets recorded

| file | holds |
|---|---|
| `run/<pipeline>/<pipeline>Summary.txt` | the diary: every submit, step start and finish, `note` lines with the numbers for your methods (samples kept, SNPs, clone threshold, K evidence, kin dropped, F<sub>ST</sub>, best model), warnings, pauses and stops |
| `run/<pipeline>/versions.txt` | comRAD, Python and its packages, every program's version (or path and checksum, when it has no version flag), every script's checksum, the loaded modules |
| `run/<pipeline>/<pipeline>.config.<date>` | every setting in effect, saved whenever one changes |
| `run/<pipeline>/<step>.log` | each step's full output |
| `run/<pipeline>/<pipeline>Report.html` | the report |

All of these go in the results tarball.
