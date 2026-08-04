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

Without `--norm`, the library size for RPM is the **weighted** one (sum of `1/NH`), on the
same scale as the per-locus counts. It costs one extra pass over the BAM (~70s for 30M
alignments); pass `--norm` to skip it. Counting alignment lines instead would deflate RPM
by the average multiplicity — 3.19x on DRR351382 (30.8M lines for 9.66M reads).

### Order of the steps

    qc -g genome.fa  ->  mapping  ->  mirna  ->  classify --mirna-loci <miRDeep2 BED>

`qc -g` is **not optional**: it standardizes the FASTA headers and maps IUPAC
ambiguity codes to N. miRDeep2 rejects a genome that has either ("not allowed
whitespaces in its first identifier", "allowed only to comprise characters
[ACGTNacgtn]"), and NCBI genomes have both.

`classify` has no miRNA class of its own — miRNA loci come from `--mirna-loci`.
Skipping it on the Arabidopsis run left 184 annotated miRNAs carrying a
confidently asserted wrong class; feeding the loci in fixed all 247. Filter the
miRDeep2 predictions by score first: predictions matching an annotated miRNA
score a median of 96.3, the rest 10.0.

### Plants vs arthropods (`--config`)

The engine carries no class names: `decider`/`features`/`reporter` only iterate over
whatever `classes` the config declares, so the biology lives entirely in the YAML.

| config | classes | validated on |
|---|---|---|
| `rules.yaml` (default) | siRNA, piRNA | Drosophila ovary (DRR351382): 94.4% recall over known piRNA clusters |
| `rules.plant.yaml` | phasiRNA, hcsiRNA, siRNA | Arabidopsis flower buds (SRR30831721): 6/6 TAS loci |

One caveat that bites when writing a config: `decision.review_confidence` is a
threshold on a **softmax**, so its strictness depends on how many classes compete.
At two classes, 0.6 means a margin of `ln(1.5) = 0.405` over the runner-up; at
three it means a margin of 1.10. The plant config uses 0.43 to get the same
0.41 margin — stricter-looking, identical in practice.

Plants have no PIWI proteins, hence no piRNA class and no ping-pong terms, and their
24-nt window belongs to heterochromatic siRNA (DCL3 -> AGO4, RdDM) rather than piRNA.
The 5' base sorts a plant small RNA into an Argonaute — 5'U to AGO1, 5'A to AGO4 —
which the `pct_5a` column and the `w_5a` weight cover.

Ping-pong and phasing are both **size-gated**: each only scores the class whose size
window already fits. That is what lets one phasing implementation serve both kingdoms —
a phased 21-nt locus reinforces a plant phasiRNA, a phased 27-nt one an arthropod piRNA.

### Validating the calls (Drosophila)

`extras/validate_piclusters.py` cross-checks a run against 90 conserved ovary piRNA
clusters (`extras/piC_ref_dm6.bed`, Release 6 coordinates — no liftover needed):

    python3 extras/validate_piclusters.py run.clusters.tsv

For a plant run, `extras/validate_plant_tas.py` does the equivalent against the
canonical Arabidopsis TAS loci (21-nt phased tasiRNA sources), reading their
coordinates straight out of a TAIR10 GFF3 by gene name:

    python3 extras/validate_plant_tas.py run.clusters.tsv TAIR10_GFF3_genes.gff

It reports recall over the reference clusters, the piRNA/siRNA read-mass enrichment
inside them, and the flamenco-vs-42AB signature contrast. On DRR351382: 94.4% recall,
3.0x piRNA enrichment, and flamenco comes out uni-strand and ping-pong-poor
(frac_top 0.900, pingpong_frac 0.059) against a dual-strand, ping-pong-rich 42AB
(0.519 / 0.132) — the published contrast, with nothing in the pipeline aware of either
locus.

### Reading the score

`score` is `log10(1/p)`: **2 means p = 0.01**, 3 means p = 0.001. It is the
weakest of the evidences that define the assigned class, so a score of 2 says
every defining evidence cleared p = 0.01.

Filter with it directly — there is no confirmed/unconfirmed flag, because the
threshold belongs to whoever is analysing.

Two caveats:

- **A low score means "not demonstrated", not "false".** The test has little
  power on shallow loci: on the Arabidopsis run the hcsiRNA loci that fail the
  5'A test still show 2.71x transposon enrichment at 19 median reads.
- **The 5' base score mixes biology with protocol.** T4 RNA ligase 1 prefers
  C > U >= A > G at the 5' end, and no tool in the field corrects for it.

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
