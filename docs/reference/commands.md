# Commands

Two commands: `comrad`, which runs the pipelines, and `comrad-jobs`, which submits and watches job
lists (the pipelines call it; it also works on its own). Below is each one's `--help`, as
`tools/makeDocs.py` captured it.

## comrad

```text
usage: comrad [-h] command ...

Run a comRAD pipeline from the folder it works in.

Pipelines:
  reads        project folder: download, trim and demultiplex, rename, split by species
  popgen       species folder: QC, mapping, clones, genotyping, structure, K, kin,
               lineages, outliers, SFS, FST, diversity, results tarball
  demography   species folder, after popgen: bootstrap SFS, moments, results tarball

Each pipeline reads <pipeline>.config in the current folder (`init` writes one to
fill in), after cluster.config (partitions, memory limit, email, extra sbatch
flags) from this folder or up to three above it (`init cluster` writes one).
Every job list has its own [section] in the config: cpus, mem, queue, slurm and
the tool's own options. State is kept in run/<pipeline>/: a .done marker and a
.log per step, and <pipeline>Summary.txt with one line per event. `run` submits
one small driver job (no time limit given) that runs the steps in order.
Steps submit their own jobs through comrad-jobs, which reruns failed jobs by
itself; a step stops the run only when something is wrong that a rerun would
not fix, and pauses it where a decision is needed (K). Submitting again picks
up after the last finished step.

positional arguments:
  command
    setup     find or install the tools the pipelines run
    defaults  print every default setting of a pipeline
    init      write <pipeline>.config (or cluster.config) to fill in
    check     check <pipeline>.config and its inputs
    run       submit the driver job
    status    last events, and whether the driver is queued
    list      every step, done or pending
    from      clear the done markers from STEP on, then submit
    step      run one step now, in this shell
    report    write run/<pipeline>/<pipeline>Report.html: warnings, pauses and stops, then the key
              tables
    summary   print open pauses and stops, then each step's notes and warnings from its last run

options:
  -h, --help  show this help message and exit

examples:
  cd /path/to/myProject
  comrad init cluster               # writes cluster.config (partitions, memory, email); once per project
  comrad init reads                 # writes reads.config; fill it in
  comrad check reads                # is everything it names there?
  comrad run reads                  # submit; species folders appear when it finishes

  cd /path/to/myProject/sint
  comrad init popgen && nano popgen.config
  comrad run popgen                 # pauses at K
  comrad status popgen              # where it is, or why it stopped
  nano popgen.config && comrad run popgen   # K = 2 in the config
  comrad from popgen clonegroups    # forget that step and everything after, resubmit
  comrad from popgen kin --mark-earlier --no-submit   # earlier steps were run by hand
  comrad step popgen fst            # run one step now, in this shell
  comrad list demography            # the steps, done or pending
  comrad report popgen              # run/popgen/popgenReport.html: warnings first, then the key tables
  comrad summary popgen --dir mcav ofav sint past   # from the project folder: what each step said, and what is waiting on you
  comrad setup --prefix /path/to/bin/comradDeps   # find or install the external tools
```

### comrad init

```text
usage: comrad init [-h] [--dir DIR] [--force] pipeline

positional arguments:
  pipeline    reads, popgen or demography (or cluster: write cluster.config)

options:
  -h, --help  show this help message and exit
  --dir DIR   folder the pipeline works in (default: here)
  --force     overwrite an existing .config
```

### comrad check

```text
usage: comrad check [-h] [--dir DIR] pipeline

positional arguments:
  pipeline    reads, popgen or demography

options:
  -h, --help  show this help message and exit
  --dir DIR   folder the pipeline works in (default: here)
```

### comrad run

```text
usage: comrad run [-h] [--dir DIR] [-p PARTITION] [-c CPUS] [--mem MEM] [-t TIME] [--here]
                  pipeline

positional arguments:
  pipeline              reads, popgen or demography

options:
  -h, --help            show this help message and exit
  --dir DIR             folder the pipeline works in (default: here)
  -p, --partition PARTITION
                        driver partition (default: [cluster] driver_queue, else long_queue, else
                        queue)
  -c, --cpus CPUS       driver CPUs (default: [cluster] driver_cpus)
  --mem MEM             driver memory (default: [cluster] driver_mem)
  -t, --time TIME       driver time limit, e.g. 2-00:00:00 (default: [cluster] driver_time, else
                        none)
  --here                run the steps in this shell instead of submitting
```

### comrad status

```text
usage: comrad status [-h] [--dir DIR] pipeline

positional arguments:
  pipeline    reads, popgen or demography

options:
  -h, --help  show this help message and exit
  --dir DIR   folder the pipeline works in (default: here)
```

### comrad summary

```text
usage: comrad summary [-h] [--dir DIR [DIR ...]] pipeline

positional arguments:
  pipeline             reads, popgen or demography

options:
  -h, --help           show this help message and exit
  --dir DIR [DIR ...]  folder(s) the pipeline works in, e.g. one per species (default: here)
```

### comrad list

```text
usage: comrad list [-h] [--dir DIR] pipeline

positional arguments:
  pipeline    reads, popgen or demography

options:
  -h, --help  show this help message and exit
  --dir DIR   folder the pipeline works in (default: here)
```

### comrad from

```text
usage: comrad from [-h] [--dir DIR] [--mark-earlier] [--no-submit] [-p PARTITION] [-c CPUS]
                   [--mem MEM] [-t TIME] [--here]
                   pipeline STEP

positional arguments:
  pipeline              reads, popgen or demography
  STEP                  step name, as `list` prints it

options:
  -h, --help            show this help message and exit
  --dir DIR             folder the pipeline works in (default: here)
  --mark-earlier        mark every step before STEP done (they were run by hand)
  --no-submit           clear the markers only
  -p, --partition PARTITION
                        driver partition (default: [cluster] driver_queue, else long_queue, else
                        queue)
  -c, --cpus CPUS       driver CPUs (default: [cluster] driver_cpus)
  --mem MEM             driver memory (default: [cluster] driver_mem)
  -t, --time TIME       driver time limit, e.g. 2-00:00:00 (default: [cluster] driver_time, else
                        none)
  --here                run the steps in this shell instead of submitting
```

### comrad step

```text
usage: comrad step [-h] [--dir DIR] pipeline STEP

positional arguments:
  pipeline    reads, popgen or demography
  STEP        step name, as `list` prints it

options:
  -h, --help  show this help message and exit
  --dir DIR   folder the pipeline works in (default: here)
```

### comrad report

```text
usage: comrad report [-h] [--dir DIR] pipeline

positional arguments:
  pipeline    reads, popgen or demography

options:
  -h, --help  show this help message and exit
  --dir DIR   folder the pipeline works in (default: here)
```

### comrad defaults

```text
usage: comrad defaults [-h] pipeline

positional arguments:
  pipeline    reads, popgen, demography or cluster

options:
  -h, --help  show this help message and exit

examples:
  comrad defaults popgen
  comrad defaults popgen | grep -A4 '\[clones\]'
  comrad defaults cluster
```

### comrad setup

```text
usage: comrad setup [-h] [--prefix PREFIX] [--only TOOL [TOOL ...]] [--no-modules] [--no-conda]
                    [--optional]

comrad setup: find or install the external tools the pipelines run.

For each tool that is not on PATH, in order:
  1. an environment module (module avail): reported, and loaded by comrad.env;
  2. conda / mamba (bioconda, conda-forge), into <prefix>/conda. With no conda or
     mamba on PATH, Miniforge (conda-forge's installer, no admin rights needed) is
     installed into <prefix>/miniforge first; --no-conda skips both;
  3. a download or build from source (GitHub releases, git, make), in <prefix>/src;
  4. pip, for the Python tools.
Every program found or installed is linked into <prefix>/bin, and
<prefix>/comrad.env puts it on PATH.
When nothing works, the manual instructions for that tool are printed. The same
list writes docs/install.md (`comrad setup --write-docs docs/install.md`).

options:
  -h, --help            show this help message and exit
  --prefix PREFIX       folder to install into (bin/, src/, Rlib/, conda/ go there)
  --only TOOL [TOOL ...]
                        just these tools (names as in docs/install.md)
  --no-modules          do not look for environment modules
  --no-conda            do not use conda/mamba, and do not install Miniforge
  --optional            also try optional tools (BaseSpace CLI)

examples:
  comrad setup --prefix /path/to/bin/comradDeps      # everything missing
  comrad setup --prefix /path/to/bin/comradDeps --only evalAdmix NGSremix
  comrad setup --prefix /path/to/bin/comradDeps --no-conda
  echo 'source /path/to/bin/comradDeps/comrad.env' >> ~/.bashrc
```

## comrad-jobs

```text
usage: comrad-jobs [-h] -n NAME [-c TEMPLATE] [-x REGEX] [--res SECTION] [-q QUEUE] [-t TIME]
                   [--long LONG] [--long-time LONG_TIME] [--cpus CPUS] [--mem MEM]
                   [--max-mem MAX_MEM] [--max-running MAX_RUNNING] [--slurm FLAGS] [--email EMAIL]
                   [--constraint CONSTRAINT] [--exclusive] [--retries RETRIES]
                   jobs

Submit a list of commands to SLURM, wait, check each one, rerun the failures.

Each line of JOBS is one command. A line's check decides whether it is done:

  check<TAB>command      explicit: the check is the text before the first tab
  command                with -c TEMPLATE: the check is TEMPLATE, with {1}, {2}
                         ... filled from the -x regexes matched against the line

A check is a bash test (`test -s out.saf.idx`, `grep -q Finished! x.log`). A job
is done when its check passes and its command exited 0 (without a check, exit 0
alone). Lines already done are not submitted again, so a rerun only does what is
missing.

Under comrad, "already done" also needs the job to have succeeded after
RT_SINCE (when the steps before it last finished, or it was rewound): every
success leaves a stamp in .rtjobs/<name>/, named after the command, with its SLURM
job and log. .rtjobs/<name>/logs lists the log behind every job in the list
(`rtlogs NAME` in the pipelines).

How jobs are submitted:
  one line    a plain sbatch job, logs <name>.o<jobid> / <name>.e<jobid>
  more lines  a SLURM job array, one task per line, at most --max-running at once
              (default: cluster max_running, else 64), split into arrays of <= 1000. Each task
              writes its own header (node, start, end, exit code) into its log, in
              .rtjobs/<name>/tasks.<...>/; when the array ends they are combined
              into <name>.e<arrayid> / .o<arrayid> (a summary table, failed tasks
              first, then one section per task) and the task folder is removed.
              It is kept if combining fails.
--cpus and --mem are per job (per task for arrays). No time limit is asked for
unless -t is given: the partition's own default applies. No partition is asked
for unless -q is given (or the config sets one): the cluster's default applies.

Under comrad, --res SECTION takes cpus, mem, queue and slurm from that
section of the pipeline config, for any of them not given on the command line,
and [cluster] gives the defaults: queue, long_queue, max_mem, max_running,
email and slurm (extra sbatch flags for every job).

After each round, only the failed jobs are resubmitted, by what went wrong:
  time limit     -> on --long (cluster long_queue), if one is set; otherwise stop
  out of memory  -> that job's --mem doubled, up to --max-mem (cluster max_mem, default 250GB)
  node failure   -> again, avoiding the node it ran on (does not use up a retry)
  anything else  -> once more (--retries), then stop
Exits 1, listing the lines still failing, if any are left. Before submitting, it
waits for jobs of the same name still in the queue (from an interrupted run).

Under comrad, a diary note per job list gives the tasks' peak memory (from sacct, if it answers), longest run and
CPU-hours.

positional arguments:
  jobs                  job list, one command per line (optionally check<TAB>command)

options:
  -h, --help            show this help message and exit
  -n, --name NAME       job name; logs are named after it
  -c, --check TEMPLATE  bash test per line, with {1} {2} ... from -x
  -x, --extract REGEX   regex matched against each line; its groups fill {1} {2} ... in order
                        (repeatable)
  --res SECTION         take cpus, mem, queue and slurm from this config section (under comrad)
  -q, --queue QUEUE     partition for the first attempt (default: cluster queue, else the
                        cluster's default)
  -t, --time TIME       time limit on --queue; default none (the partition's own)
  --long LONG           partition after a time-limit failure (default: cluster long_queue; none =
                        stop)
  --long-time LONG_TIME
                        time limit on --long; default none
  --cpus CPUS           CPUs per job (per array task) (default: the section's cpus, else 1)
  --mem MEM             memory per job (per array task), e.g. 8GB (default: the section's mem)
  --max-mem MAX_MEM     out-of-memory retries stop above this (default: cluster max_mem, else
                        250GB)
  --max-running MAX_RUNNING
                        array tasks running at once (default: cluster max_running, else 64)
  --slurm FLAGS         extra sbatch flags as 'key=value key ...' (added to cluster and section
                        slurm)
  --email EMAIL         mail on failure to this address (default: cluster email)
  --constraint CONSTRAINT
                        sbatch --constraint, e.g. amd (same as --slurm constraint=amd)
  --exclusive           a single job takes the whole node (only for jobs that need it)
  --retries RETRIES     attempts after the first for anything but node failures (default 3)

examples:
  # angsd per lineage; done when the .saf.idx exists
  comrad-jobs sfsJobs -n sintSFS -c 'test -s {1}.saf.idx' -x='-out (\S+)' --mem 32GB

  # NGSadmix job list from ngsAdmixLauncher.py (-P 10 each), at most 40 at once
  comrad-jobs sintNgsAdmix -n sintNgsAdmix -c 'test -s {1}.qopt' -x='-o (\S+)' --cpus 10 --mem 8GB --max-running 40

  # explicit checks: one 'check<TAB>command' per line
  printf 'test -s a.out\tsleep 5; echo hi >a.out\n' >jobs
  comrad-jobs jobs -n demo

  # one NGSremix job: a plain sbatch job with 8 CPUs and 100 GB
  comrad-jobs remixJob -n sintRemix --cpus 8 --mem 100GB

  # bootstrap spectra on AMD nodes, 20 CPUs each, on a named partition
  comrad-jobs bootJobs -n sintBoot -c 'test -s {1}' -x='>(\S+\.sfs)' --cpus 20 --mem 8GB -q shortq7 --slurm 'constraint=amd'

  # under comrad: cpus, mem, queue and slurm from the config's [boot] section
  comrad-jobs bootJobs -n sintBoot --res boot

  # a regex starting with '-' must be written -x='...' so it is not read as an option
```
