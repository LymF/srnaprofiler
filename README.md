# srnaprofiler

Small RNA-seq analysis and cluster classification pipeline (quality control,
mapping, coverage-based cluster classification into siRNA / piRNA, and miRNA
analysis).

> This README covers **installation only**, so the tool can be installed and
> tested on real data.

## Requirements

- **Python ≥ 3.9**
- **Python packages** (see `requirements.txt`): `pysam`, `pyyaml`, `biopython`, `pytest`
- **External tools** (only needed for the full pipeline; the `classify` module
  needs just the Python packages above):
  - [Bowtie 1](https://bowtie-bio.sourceforge.net/index.shtml) and `bowtie-build`
  - [SAMtools](http://www.htslib.org/)
  - [fastp](https://github.com/OpenGene/fastp) **or** [cutadapt](https://cutadapt.readthedocs.io/) (read trimming)
  - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (optional QC report)
  - [R](https://www.r-project.org/) with `Rscript` (plots)
  - [miRDeep2](https://github.com/rajewsky-lab/mirdeep2) (miRNA module)

## Installation

### 1. Get the code

```bash
git clone git@github.com:LymF/srnaprofiler.git
cd srnaprofiler
```

### 2. Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows (PowerShell): .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 3. External tools (recommended via conda / bioconda)

```bash
conda create -n srnaprofiler -c bioconda -c conda-forge \
    bowtie samtools fastp fastqc mirdeep2 r-base
conda activate srnaprofiler
pip install -r requirements.txt
```

### 4. Verify the installation

```bash
# run the test suite (all tests should pass)
pytest -q

# show the command-line help
python3 srnaprofiler.py --help
```

## Modules

| Command    | Purpose                                                         |
|------------|----------------------------------------------------------------|
| `qc`       | Read QC, trimming, and genome standardization                  |
| `mapping`  | Align small RNA reads to the genome with Bowtie 1              |
| `classify` | Coverage-based cluster calling and siRNA/piRNA classification  |
| `mirna`    | miRNA prediction and quantification with miRDeep2             |
| `plots`    | Distribution and pattern plots                                 |

Run `python3 srnaprofiler.py <command> --help` for per-module options.

### classify — calibrating `--mincov`

Read counts are weighted by `1/NH`, so honest per-locus abundances are much lower than
raw alignment pileups. With a permissive `--mincov`, shallow mapping tails seed tens of
thousands of 1–3 read islands that can never pass `decision.min_reads`, and most loci end
up flagged `review`. Raising `--mincov` removes that noise at the source — it does not
cost signal.

Sweep on DRR351382 (Drosophila ovary, 30.8M alignments; `min_reads: 10`):

| `--mincov` | clusters | review rate | kept loci | read mass kept | loci >=50 reads |
|-----------:|---------:|------------:|----------:|---------------:|----------------:|
| 2          | 83,464   | 79.3%       | 17,271    | 0.704          | 7,324           |
| 5          | 64,828   | 66.2%       | 21,926    | 0.697          | 8,511           |
| 10         | 58,909   | 56.4%       | 25,673    | 0.685          | 9,395           |
| 20         | 48,278   | 44.2%       | 26,932    | 0.680          | 10,120          |
| 30         | 41,015   | 35.2%       | 26,587    | 0.669          | 10,298          |
| 50         | 30,176   | 15.1%       | 25,619    | 0.664          | 9,787           |

The 1,000 most abundant loci called at `--mincov 2` are recovered in full (loci and read
mass) at every threshold, including 50. Kept loci and substantial loci peak at 20–30 while
the retained read mass stays flat, so the default is **10** (safe for shallow libraries);
deep libraries (>=30M alignments) do better with **20–30**.

Raising `decision.min_reads` instead is the wrong lever: it only relabels loci as `review`
without retaining more signal (read mass kept drops to 0.671 / 0.644 / 0.615 at
`min_reads` 20 / 50 / 100 for `--mincov 10`).

### classplots — per-class plots

Re-aggregates a classify run's alignments grouped by the assigned class and produces
per-class publication signal plots (size distribution, pattern, ping-pong, phasing,
sequence logos) plus a matplotlib QC summary. Each **chosen locus** — the `--top-n` most
abundant per class, and every `--region` — also gets an IGV-style coverage track (sense up,
antisense down; raw depth).

    srnaprofiler classplots --clusters run.clusters.tsv --bam sample.sorted.bam \
        -o classplots_out --prefix run --top-n 5

Plot a specific coordinate span (regardless of class) with one or more `--region`
(the `name` column of the `.clusters.tsv` is exactly `chrom:start-end`):

    srnaprofiler classplots --clusters run.clusters.tsv --bam sample.sorted.bam \
        -o classplots_out --region chr1:100-400 --region X:9000-9500

`--classes` is case-insensitive (`piRNA,siRNA` or `pirna`). The R signal plots (size
distribution, pattern, ping-pong, phasing) need `Rscript` (skipped with a warning if
absent). Sequence logos are drawn by **logomaker** (`pip install logomaker`), and the
coverage track and QC summary use matplotlib — all pure Python, no Ghostscript or external
CLI. `--logo-format` picks the logo format (pdf/svg/png/eps; default pdf).
