# Quick start

One project, from raw reads to results. Every command runs on the cluster's login node;
comRAD submits its own jobs. Each step is explained in the [guide](guide/setup.md).

## Once per cluster

```bash
python3 -m pip install --user git+https://github.com/RyanEckert/comRAD.git
comrad setup --prefix /path/to/storage/comradDeps          # finds or installs ANGSD, bowtie2, ...
echo 'source /path/to/storage/comradDeps/comrad.env' >> ~/.bashrc && source ~/.bashrc
```

Index each reference genome once:

```bash
bowtie2-build sintGenome.fa sintGenome.fa && samtools faidx sintGenome.fa
```

## Once per project

```bash
cd /path/to/myProject                # rawReads/, sampleID.csv and metadata.csv go here
comrad init cluster                  # partitions, largest node memory, email
nano cluster.config
```

## Reads: raw pools to one file per sample

```bash
comrad init reads && nano reads.config      # SAMPLE_SHEET, SPECIES
comrad check reads
comrad run reads
comrad status reads                         # where it is, or why it stopped
```

## Popgen: one species

```bash
cd sint
comrad init popgen && nano popgen.config    # SP, GENOME_FASTA, META, TECHREPS, ...
comrad check popgen
comrad run popgen
comrad summary popgen                       # what each step found so far, and what waits on you
```

It pauses twice. (Set `CLONE_PAUSE = no` and `K` in the config beforehand and it never does.)
First at the **clone dendrogram**:

```bash
comrad status popgen                        # the last line says what to look at
nano popgen.config                          # CLONE_PAUSE = no  (accept)  or  CLONE_CUT = 0.08  (your own)
comrad run popgen
```

Then at **K**:

```bash
less run/popgen/K.log                       # the reads on K, side by side
# and the plots: pcangsd/sintPca.pdf, ngsAdmix/sintAdmix.pdf, ngsAdmix/sintResiduals.pdf
nano popgen.config                          # K = 2
comrad run popgen
```

When it finishes, `sintResults.tar.gz` holds the results and `run/popgen/popgenReport.html`
summarises the run, warnings first.

## Demography (optional)

```bash
comrad init demography && nano demography.config   # SP, MU
comrad run demography
```
