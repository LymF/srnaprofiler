# srnaprofiler

Small RNA-seq analysis and cluster classification pipeline (quality control,
mapping, coverage-based cluster classification into siRNA / piRNA, and miRNA
analysis).

> This README covers **installation only**, so the tool can be installed and
> tested on real data.

## Installation

```bash
git clone git@github.com:LymF/srnaprofiler.git
cd srnaprofiler
conda env create -f environment.yml
conda activate srnaprofiler

python3 srnaprofiler.py setup --tools                       # toolchain
python3 srnaprofiler.py setup --data --rfam-dir ~/ncrna_db  # references
python3 srnaprofiler.py setup --check                       # verify
```

That is the whole installation. What each step does, and why it is a command
rather than a list of instructions, follows.

### `setup --check`

Reports, per module, whether a run would work here — and checks each requirement
through **the interface that actually breaks**, which is not the same as asking
whether a package is installed:

```
mirna
  [ok  ] bowtie       .../bin/bowtie
  [ok  ] RNAfold      .../bin/RNAfold
  [ok  ] viennarna    .../python3.13/site-packages/RNA/__init__.py
  [warn] miRScore     optional: validation is skipped without it

modules that can run here:
  classify   yes
  mirna      yes
```

`RNAfold` and `viennarna` are two rows for one package on purpose: the miRDeep
engines call the **binary**, miRScore imports the **Python binding**, and
`conda install viennarna` on Python 3.14 gives you the first without the second.
That is why `environment.yml` pins **python=3.13** — it is what lets one
environment serve every module.

It touches nothing, so it is safe to run at any time.

### `setup --tools`

Installs the toolchain into the **active** conda environment (it refuses when
none is active, rather than writing into `base`), fetches miRDeep-P2 into
`./tools` beside the code, patches it, and finishes by running `--check`.

miRDeep-P2 needs three fixes, and all three fail **silently** — the pipeline runs
to completion, writes `end...` in its log, and predicts zero:

- a path bug in its own script (writes `$out/$n/$n.ori.fa`, reads `$out/$n.ori.fa`);
- a missing `rfam_index`, which is a separate download;
- gzipped input, which it cannot read.

miRScore adds a sixth, and it is the one that bites hardest: it derives the
library key from the FASTQ name **cut at the first dot** while naming the BAM
with the whole stem, so `sample.ncrna_filtered.fastq` counts zero reads at every
locus with the reads sitting aligned in the BAM. That is the name `mapping`
emits, so chaining the two modules used to yield an empty result and no error.
The `mirna` module stages the reads under a dot-free name, and treats zero reads
at *every* locus as a failure rather than a finding.

`setup` handles the first two and the `mirna` module handles the third. The patch
keeps the original beside it, is idempotent, and **refuses a script whose known
bug it cannot find** — a later version may fix it differently, and rewriting
third-party code on a guess is worse than stopping.

### `setup --data`

Builds the Rfam structural RNA index (~1.8 GB downloaded, then `bowtie-build`;
skipped if already built), links it into miRDeep-P2 so the 2.9 GB is not held
twice, and records the paths in `~/.config/srnaprofiler/paths.yaml`. After this
`mapping` and `mirna` find their references with no flags.

### Reference data you provide

| data | used by | note |
|---|---|---|
| genome FASTA + bowtie index | `mapping`, `classify`, `mirna` | your assembly |
| miRBase GFF3 + hairpin | `classify --mirna-loci` | when the species has one |

**Check that a miRBase annotation belongs to your genome before using it.**
miRBase publishes each species against whatever assembly it used, and mixing
assemblies fails silently: loci are carved at the wrong coordinates, miRNA is
asserted over them, and every validation still returns a number.

```bash
python3 extras/check_mirbase_assembly.py genome.fa osa.gff3 hairpin.fa osa
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

Island depth is on the same `1/NH` scale as every other count, so `--mincov 1.0` means one
weighted read per million of library — genuinely ten reads deep on a 10M-read library.

Measured at the defaults on DRR351382 (Drosophila ovary) and SRR30831721 (Arabidopsis
flower buds), against the 90 reference piRNA clusters and the six TAS loci:

| | Drosophila | Arabidopsis |
|---|---|---|
| loci | 25,661 | 29,407 |
| median span | 546 bp | 605 bp |
| weighted reads per locus | 63 | 124 |
| per-base coverage in a locus | 4.2x | 6.4x |
| distinct sequences (the `n` of the tests) | 163 | 115 |
| review rate | 14.2% | 3.7% |
| known loci | 84/90 piCs | 6/6 TAS |
| piRNA read mass in the piCs | 3.3x | — |

`--mincov` sweep on Drosophila (piRNA enrichment is flat across it, so the threshold does
not discard piRNA mass — it concentrates the same mass into fewer, deeper loci):

| `--mincov` (RPM) | loci | reads/locus | per-base coverage | distinct seqs | review |
|---:|---:|---:|---:|---:|---:|
| 0.3 | 53,547 | 33 | 2.1x | 88 | 31.9% |
| 0.5 | 41,473 | 38 | 2.9x | 103 | 19.7% |
| **1.0** | 25,661 | 63 | **4.2x** | **163** | 14.2% |

**Multi-mapping reads** contribute `1/NH` at each place they align — the same scheme as
tiny-count's *normalize by genomic hits* and ShortStack's `f` mode. Tools that instead place
each read at a single position report higher per-locus counts, by roughly the average
multiplicity. `unique_frac` reports how much of a locus's mass comes from uniquely-mapping
reads: 57% of Drosophila piRNA loci have none, so their whole signature rests on repeats.

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

`classify` has no miRNA class of its own — miRNA loci come from `--mirna-loci`,
and each interval given there is **carved out of the cluster around it** before
the class is decided. A hairpin is ~90 bp and a locus is 600–1000 bp, so without
carving the miRNA is a seventh of a locus whose class comes from everything else
in it: every covered precursor was labelled siRNA or piRNA, and 21% of the plant
phasiRNA class sat on an annotated miRNA. With carving and a curated annotation,
all 244 covered precursors in Arabidopsis and all 134 in Drosophila come out
miRNA, none flagged for review, and the contamination of the other classes goes
to zero. The flanks survive as loci of their own rather than being discarded.

Two consequences worth knowing:

- **Give it curated annotation when one exists.** A carved interval becomes a
  locus with an asserted label, so a wrong prediction now produces a wrong
  locus, where before it merely failed to label one. Of the loci carved from raw
  miRDeep2 predictions only 57% matched an annotated miRNA. Filter by score
  first: predictions matching an annotation score a median of 96.3, the rest 10.0.
- **Read abundance figures accordingly.** miRNA loci held 33.7% of the read mass
  of the Drosophila library. Any per-class share computed without `--mirna-loci`
  counts that mass as siRNA or piRNA.

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

`phase_pvalue` is corrected for the window scan. The test slides a window of 11 phase
cycles across the locus and reports the most convincing one, so its p-value pays a
Bonferroni price for the windows actually tried — reporting the best of many as if it were
a single test made the false-positive rate grow with locus length. On the Arabidopsis
library the rate at `p <= 0.01` was 9.8% before the correction and 0.89% after, against
the 1% a calibrated null gives, and annotated PHAS regions went from 2.0x to 6.7x enriched
for significant phasing over the background.

### Filtering structural RNA (`--ncrna-db`)

tRNA, rRNA, snoRNA and snRNA fragments are abundant in any small RNA library and none of them
is siRNA or piRNA. Measured here, loci overlapping them are 2.26% of the Arabidopsis run and
1.56% of the Drosophila one, but carry **36.5%** and **25.7%** of the read mass, and in plants
they account for 16.1% of the phasiRNA class.

Rejecting those loci after the fact was measured and failed — the best threshold still lost a
quarter of the experimentally characterised TAS loci. Removing the *reads* instead means the
locus never exists and its abundance is never wrong, which is what the reference tools do.

#### Building the Rfam index (once per machine)

    THREADS=8 bash extras/fetch_ncrna_db.sh ~/ncrna_db

That is the whole tutorial; what follows is what it does and what to expect, so a run that
goes wrong is recognisable.

1. **Selects the families.** Downloads Rfam's `family.txt` and keeps every accession whose
   *type* field is rRNA, tRNA, snRNA or snoRNA — **829 families** at the release this was
   written against. Selection is by type, not by a fixed list, so families added later are
   picked up, and not by matching the words anywhere in the line, which would drag in RNase P,
   tmRNA, telomerase and group I/II introns because their descriptions mention the words.
   Fewer than 200 families aborts the run.
2. **Downloads one FASTA per family** from `fasta_files/`. Rfam publishes no single FASTA of
   its sequences — the whole-database files at the top of the FTP are covariance models.
   Roughly 1 h on a home connection; a family that 404s is reported and skipped, and fewer
   than 100,000 sequences in total aborts, so a truncated download cannot become a thin index.
3. **Runs `bowtie-build`.** About 1 h 45 with 8 threads.

What it costs, measured:

| | |
|---|---|
| `ncrna.fa` | 1.8 GB, 7.3 M sequences |
| index (6 files) | 2.9 GB |
| peak RAM to build | fits in 8 GB, no deduplication needed |
| filtering 24 M reads | ~6 min with 8 threads, ~30 s of it loading the index |

**Put the index on a local disk.** On a Windows drive mounted through WSL (`/mnt/...`) the
index is slow enough that `bowtie-inspect` appears to hang. Copying the six `.ebwt` files to
the Linux filesystem takes 20 s and removes the problem.

**The index is not deduplicated,** on purpose: 37.7% of the Rfam sequences are exact copies of
another, so it could be a third smaller, but it builds within the memory available as it is,
and a deduplication step is one more thing that can go wrong quietly.

Then filter while mapping:

    srnaprofiler mapping -p sample -o out --reads clean.fastq.gz \
        --index-prefix genome --ncrna-db ~/ncrna_db/ncrna

The log reports how many reads left, which is the number to check: 3% and 90% removed are
not the same run. Single-end only — `--un` splits pairs into separate files.

Measured on the two test libraries: **41.5%** of the Arabidopsis reads and **31.0%** of the
Drosophila ones are structural RNA. The read mass sitting on annotated structural loci drops
from 36.5% to 3.5% and from 25.7% to 4.4%, phasiRNA contamination from 16.1% to 9.0%, with
the reference loci untouched (TAS 6/6, piC 85/90, miRBase 242/242 and 134/134).

**Compare runs at the same denominator.** `--mincov` is RPM, so a library that lost a third
of its reads gets a lower absolute depth threshold and reports more, smaller, shallower loci.
That is the threshold, not the filter: with `--norm` set to the unfiltered library size the
loci come back to what they were (median span 594 vs 589 bp).

**Sequences come from Rfam rather than from your genome's annotation on purpose:** a non-model
organism has no structural RNA annotated to extract, while tRNA and rRNA are conserved enough
for a 20-30 nt read to match across species at one mismatch.

**miRNA families are excluded.** Rfam holds 1598 of them, 38% of the database, and miRNA is a
class this tool calls — filtering them out would destroy it. sRNAminer names its equivalent
parameter the same way: *ncRNA [rRNA, tRNA, snoRNA and snRNA] without miRNA database*.

Rfam ships no single FASTA of its sequences: the whole-database files at the top of its FTP
(`Rfam.cm.gz`, `Rfam.tar.gz`, both 44 MB) are *covariance models*, which serve `cmscan`
searching a genome for structure — hours on a plant genome. Sequences live one file per
family.

Two caveats worth knowing:

- The two snRNA families are less conserved across clades than tRNA and rRNA, so in a distant
  organism they will catch less. That is a silent false negative, not an error you will see in
  the log.
- The script does not pin an Rfam release: accessions are permanent (RF00005 is tRNA in every
  release) and the path follows `CURRENT`. What can still break is EBI moving the path, and
  the script is written to fail loudly if it does — it checks that each family arrives with
  at least 100 sequences, so a truncated download aborts instead of quietly building a thinner
  index.

### Leaving out the organelle (`--organelle-contigs`)

Chloroplast and mitochondrion give 336 loci and **7.9%** of the read mass in the Arabidopsis
run (27 loci and 0.6% in Drosophila), almost all called `hcsiRNA`. Plants have no organellar
siRNA pathway, so that is degradation being labelled. Those loci score low, but low score is
not absence: they are in the table, and in any abundance summed per class.

    srnaprofiler classify --sam sample.bam -g genome.fa -o out -p sample \
        --organelle-contigs NC_000932.1,NC_037304.1

Contigs are named, never guessed from the FASTA description. RefSeq headers do carry the word
(`>NC_000932.1 ... chloroplast, complete genome`), but miRDeep2 refuses spaces in the
identifier, so a genome prepared for it has already lost it, and an assembly from elsewhere
may carry no description at all — the guess would work here and fail silently on the next
genome. A name that is not in the BAM aborts the run instead of being ignored.

The library size still counts the excluded reads, so RPM stays comparable with runs that
kept them.

### Formulas

Every quantity the pipeline reports, with its formula, its null and its known limitations,
is in [`docs/FORMULAS.md`](docs/FORMULAS.md). How each one compares to what the reference
tools do is in `docs/superpowers/specs/2026-08-07-signature-literature-review.md`.

Two things worth knowing before reading any score:

- **Pass `-g/--genome`.** The 5' base scores (`s_1u`, `s_5a`) are tested against the locus's
  own DNA composition: a cut placed at random in a locus yields a 5'U as often as that locus
  carries T. Without the genome the null falls back to a uniform 0.25, which no genome has —
  measured, it is 0.272 in Drosophila and 0.332 in Arabidopsis — and a T-rich transposon then
  looks 1U-biased with no selection at all.
- **`score` ranks by evidence strength, in `[0,1]`.** It is the weakest of the class's
  defining evidences, each read as a normalised Wilson lower bound: 0.61 means that evidence
  is at least 61% of the way from its null to certainty. It grows with the strength of the
  signal, not with the size of the locus. The `s_*` columns still hold `-log10(p)` and answer
  the different question of whether the observation would be surprising under the null.

### Validating the calls (Drosophila)

`extras/validate_piclusters.py` cross-checks a run against 90 conserved ovary piRNA
clusters (`extras/piC_ref_dm6.bed`, Release 6 coordinates — no liftover needed):

    python3 extras/validate_piclusters.py run.clusters.tsv

`extras/validate_mirna.py` checks either kingdom against miRBase, whose GFF3 coordinates
match the assemblies already in use (`dme` on Release_6, `ath` on TAIR10 — no liftover):

    curl -sSLO https://www.mirbase.org/download/dme.gff3
    python3 extras/validate_mirna.py run.clusters.tsv dme.gff3

It reports recall and — the part the other validators cannot ask — what the covered
precursors are called when `--mirna-loci` is *not* given, plus how much of each class sits
on an annotated miRNA. Validating one class at a time hides a regression in another: the
miRNA route broke and both the piRNA and the PHAS validator passed clean, because neither
reference contains a miRNA.

`extras/validate_phas.py` checks a plant run against the PHAS loci annotated by sRNAanno
(`extras/phas21_ref_at_TAIR10.gff3`, 35 loci for Arabidopsis; the server also has PHAS24
sets for the grasses). It reports recall, how well our phasing test separates annotated
PHAS regions from the background, and the agreement with the reference's own p-value,
locus by locus. Six TAS genes were too narrow a base to calibrate phasing on.

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
