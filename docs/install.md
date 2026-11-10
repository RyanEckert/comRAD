# Install

comRAD needs a Linux cluster with SLURM and Python 3.9 or newer. Installing it is two steps:
comRAD itself with pip, then the external programs it runs with `comrad setup`.

## 1. comRAD

```bash
module load python                 # whichever Python you will run with (if your cluster uses modules)
python3 -m pip install --user git+https://github.com/RyanEckert/comRAD.git
comrad -h
```

This installs `comrad`, `comrad-jobs`, every script the pipelines call, and the Python packages
they use (numpy, scipy, matplotlib, moments, demes, PCAngsd and cutadapt). R is not needed.

!!! tip "Small home folder?"
    `pip --user` installs into `~/.local`. To put it on bigger storage instead, set
    `PYTHONUSERBASE` before installing, and put its `bin` on your `PATH`:

    ```bash
    echo 'export PYTHONUSERBASE=/path/to/storage/.local' >> ~/.bashrc
    echo 'export PATH=$PYTHONUSERBASE/bin:$PATH' >> ~/.bashrc
    source ~/.bashrc
    ```

Load the same Python before installing that you will run comRAD with: pip writes that
interpreter into each script's first line. If pip itself is missing, `python3 -m ensurepip
--upgrade --user` installs it. To update comRAD later, run the same `pip install` with
`--upgrade`; to install a tagged version, add `@v0.1.0` after `.git`.

## 2. The programs comRAD runs

```bash
comrad setup --prefix /path/to/storage/comradDeps
echo 'source /path/to/storage/comradDeps/comrad.env' >> ~/.bashrc
source ~/.bashrc
```

Below, `DIR` is that `--prefix` folder; put it on storage with space, not a small home folder.
For each program that is not already on your `PATH`, `comrad setup` tries, in order:

1. an **environment module** (`module avail`); `comrad.env` loads the ones it finds;
2. **conda/mamba** from bioconda, into `DIR/conda`. With no conda or mamba on the cluster, it
   first installs [Miniforge](https://github.com/conda-forge/miniforge) into `DIR/miniforge`
   (no admin rights needed); `--no-conda` skips conda altogether;
3. a **download or build from source** in `DIR/src` (needs git, make and a C/C++ compiler);
4. **pip**, for the Python tools.

Everything it installs is linked into `DIR/bin`, which `comrad.env` puts on your `PATH`.
Anything it could not install, it prints the steps for; they are also below.
`comrad check popgen` (in a species folder) lists anything still missing.

This page is written by `comrad setup --write-docs docs/install.md`, from the same list setup uses.

## The programs

| tool | programs | needed by | source |
|---|---|---|---|
| ANGSD | `angsd`, `realSFS`, `thetaStat` | popgen, demography | https://github.com/ANGSD/angsd |
| NGSadmix | `NGSadmix` | popgen | https://github.com/ANGSD/angsd (built in misc/) |
| bowtie2 | `bowtie2`, `bowtie2-build` | popgen | https://github.com/BenLangmead/bowtie2 (prebuilt release binaries) |
| samtools | `samtools` | popgen | https://www.htslib.org/download/ |
| cutadapt | `cutadapt` | popgen | https://cutadapt.readthedocs.io |
| PCAngsd | `pcangsd` | popgen | https://github.com/Rosemeis/pcangsd |
| evalAdmix | `evalAdmix` | popgen | https://github.com/GenisGE/evalAdmix |
| NGSremix | `NGSremix` | popgen | https://github.com/KHanghoj/NGSremix |
| ngsParalog | `ngsParalog` | popgen | https://github.com/tplinderoth/ngsParalog |
| BaseSpace CLI | `bs` | reads, only to download from BaseSpace (BS_PROJECTS) | https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview |

### ANGSD

Module (most clusters): `module avail angsd`, then add `module load <name>` to `~/.bashrc`.

Conda:

```bash
conda install -c conda-forge -c bioconda angsd
```

From source (needs git, make, a C/C++ compiler and zlib, bzip2, xz (lzma) and libcurl development files):

```bash
cd DIR/src
git clone --recurse-submodules https://github.com/ANGSD/angsd.git && cd angsd && make
```

Then put `angsd/angsd`, `angsd/misc/realSFS`, `angsd/misc/thetaStat` on PATH (for example, link them into `DIR/bin`).

### NGSadmix

Module (most clusters): `module avail ngsadmix`, then add `module load <name>` to `~/.bashrc`.

From source (needs git, make, a C/C++ compiler and zlib, bzip2, xz (lzma) and libcurl development files):

```bash
cd DIR/src
[ -d angsd ] || git clone --recurse-submodules https://github.com/ANGSD/angsd.git; cd angsd && make
```

Then put `angsd/misc/NGSadmix` on PATH (for example, link them into `DIR/bin`).

### bowtie2

Module (most clusters): `module avail bowtie2`, then add `module load <name>` to `~/.bashrc`.

Conda:

```bash
conda install -c conda-forge -c bioconda bowtie2
```

From source (needs git, make, a C/C++ compiler and perl (bowtie2's wrapper script), and unzip or python3):

```bash
cd DIR/src
V=2.5.4; D=bowtie2-$V-linux-$(uname -m); U=https://github.com/BenLangmead/bowtie2/releases/download/v$V/$D.zip; (curl -fsSL -o $D.zip $U || wget -q -O $D.zip $U) && rm -rf $D && (unzip -oq $D.zip || python3 -m zipfile -e $D.zip .) && chmod +x $D/bowtie2* && ln -sfn $D bowtie2
```

Then put `bowtie2/bowtie2`, `bowtie2/bowtie2-build` on PATH (for example, link them into `DIR/bin`).

### samtools

Module (most clusters): `module avail samtools`, then add `module load <name>` to `~/.bashrc`.

Conda:

```bash
conda install -c conda-forge -c bioconda samtools
```

From source (needs git, make, a C/C++ compiler and zlib, bzip2 and xz (lzma) development files):

```bash
cd DIR/src
V=1.21; U=https://github.com/samtools/samtools/releases/download/$V/samtools-$V.tar.bz2; (curl -fsSL -o samtools-$V.tar.bz2 $U || wget -q -O samtools-$V.tar.bz2 $U) && tar -xjf samtools-$V.tar.bz2 && cd samtools-$V && ./configure --without-curses --prefix=$PWD/inst && make -j4 && make install && cd .. && ln -sfn samtools-$V samtools
```

Then put `samtools/inst/bin/samtools` on PATH (for example, link them into `DIR/bin`).

### cutadapt

Module (most clusters): `module avail cutadapt`, then add `module load <name>` to `~/.bashrc`.

Conda:

```bash
conda install -c conda-forge -c bioconda cutadapt
```

pip:

```bash
python3 -m pip install --user cutadapt
```

### PCAngsd

Module (most clusters): `module avail pcangsd`, then add `module load <name>` to `~/.bashrc`.

Conda:

```bash
conda install -c conda-forge -c bioconda pcangsd
```

pip:

```bash
python3 -m pip install --user git+https://github.com/Rosemeis/pcangsd.git
```

### evalAdmix

Module (most clusters): `module avail evaladmix`, then add `module load <name>` to `~/.bashrc`.

From source (needs git, make, a C/C++ compiler and zlib development files):

```bash
cd DIR/src
git clone https://github.com/GenisGE/evalAdmix.git && cd evalAdmix && make
```

Then put `evalAdmix/evalAdmix` on PATH (for example, link them into `DIR/bin`).

### NGSremix

Module (most clusters): `module avail ngsremix`, then add `module load <name>` to `~/.bashrc`.

From source (needs git, make, a C/C++ compiler and zlib development files):

```bash
cd DIR/src
git clone https://github.com/KHanghoj/NGSremix.git && cd NGSremix/src && make -f CPP_Makefile
```

Then put `NGSremix/src/NGSremix` on PATH (for example, link them into `DIR/bin`).

### ngsParalog

Module (most clusters): `module avail ngsparalog`, then add `module load <name>` to `~/.bashrc`.

From source (needs git, make, a C/C++ compiler and zlib development files):

```bash
cd DIR/src
git clone https://github.com/tplinderoth/ngsParalog.git && cd ngsParalog && make
```

Then put `ngsParalog/ngsParalog` on PATH (for example, link them into `DIR/bin`).

### BaseSpace CLI

Module (most clusters): `module avail basespace`, then add `module load <name>` to `~/.bashrc`.

See https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview.

## Conda without conda

With no conda or mamba on PATH, setup installs Miniforge (conda-forge's installer; no admin
rights, conda-forge channels only); by hand that is:

```bash
curl -fsSL -o Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh
bash Miniforge3.sh -b -p DIR/miniforge
DIR/miniforge/bin/mamba create -y -p DIR/conda -c conda-forge -c bioconda angsd bowtie2 samtools
```

## Python packages

numpy, scipy, matplotlib, moments-popgen, demes, PCAngsd and cutadapt are installed with comRAD.
R is not needed.

