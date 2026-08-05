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

### classify — calibrating `--mincov` and `--pad`

`--mincov` is in **reads per million**, converted once per run into an absolute depth
using the library size. As a raw count it did not transfer between libraries: a 10
calibrated on a 30.8M-alignment run means something else at any other depth.

A locus is an island of positions above that depth, extended by `--pad` on **both** sides
(so the tolerated gap is `2*pad`). Neighbours merge only when their profiles are
compatible: dominant read lengths in the same class size window (from the config, so a
25-nt island joins a 26-nt one but not a 21-nt one) and strand fractions within 0.5.
Neighbours that do not merge split the gap between them — loci are always disjoint,
because a read is assigned to every locus it touches.

Compatibility is judged by the class window holding most of the read mass, not by
the modal length: a raw island has ~5 weighted reads, so its mode flips with one read
either way, while the dominant window is stable (76% of the islands under 10 reads have
one). Islands with no dominant window do not block a merge.

Measured at the defaults on DRR351382 (Drosophila ovary) and SRR30831721 (Arabidopsis
flower buds), against the 90 reference piRNA clusters and the six TAS loci:

| | Drosophila | | Arabidopsis | |
|---|---|---|---|---|
| | before (pad 75) | now | before | now |
| loci | ~186,000 | 45,666 | 23,248 | 46,492 |
| median span | 88 bp | **335 bp** | ~100 bp | **532 bp** |
| reads per locus | 34 | 36 | 34 | **74** |
| review rate | 79% | **29.2%** | 33.1% | **8.2%** |
| known loci | 85/90 piCs | 85/90 piCs | 6/6 TAS | 6/6 TAS |
| piRNA read mass in piCs | 3.0x | **3.4x** | — | — |
| siRNA read mass in piCs | 1.4x | **1.6x** | — | — |

**On piC recall as a criterion:** it asks whether a 10-40 kb reference region holds at
least one non-review piRNA locus, so it *rewards fragmentation* — the more loci a region
is cut into, the likelier one of them is piRNA-sized. A run that scored 87/90 did so with
111-bp loci; with larger loci the dominant population wins and genuinely 22-nt sub-regions
inside those references are correctly called siRNA (three of them hold 88-96% of their
read mass at 20-22 nt, with 0-2% at 24-31 nt). Prefer the enrichment figures, which sum
read mass instead of counting occurrences, when calibrating clustering parameters.

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
- **miRNA loci always score 0.** The class comes from `--mirna-loci`, defined by
  precursor structure rather than by read features, so it has no defining
  evidence to score. Filtering by score would drop every miRNA.
- **`s_pingpong` and `s_density` are diagnostic only — no shipped config uses
  them as defining evidence.** Both nulls are weak: the ping-pong test treats
  the 20 offset bins as if they were independent trials, but the denominator
  is actually a sum of 1/NH weight *products* across read pairs, not a count
  of independent draws, so the p-value understates the true uncertainty. The
  density test is conditioned on the locus existing in the first place — a
  locus is only in the table because reads piled up there, so testing that
  pile-up against a uniform-genome null is circular. Both will be revisited
  once clustering no longer selects loci by the same read pile-up being tested.

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
