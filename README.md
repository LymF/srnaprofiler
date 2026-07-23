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

### classplots — per-class plots

Re-aggregates a classify run's alignments grouped by the assigned class and produces
publication signal plots (size distribution, pattern, ping-pong, phasing, coverage,
sequence logos) plus a matplotlib QC summary.

    srnaprofiler classplots --clusters run.clusters.tsv --bam sample.sorted.bam \
        -o classplots_out --prefix run --top-n 5

Plot a specific coordinate span (regardless of class) with one or more `--region`:

    srnaprofiler classplots --clusters run.clusters.tsv --bam sample.sorted.bam \
        -o classplots_out --region chr1:100-400 --region X:9000-9500

Signal plots need `Rscript` (skipped with a warning if absent) and the `weblogo` CLI for
logos (`pip install weblogo`); the QC summary is pure Python.
