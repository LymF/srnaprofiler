# sRNAProfiler

Small RNA loci discovery and classification from sequencing data. It calls loci from
coverage, measures the signatures each biogenesis pathway leaves behind (size, strand
balance, 5′ nucleotide, ping-pong overlap, phasing), and assigns a class with a transparent,
editable rule set. Plants and animals share one engine, and the biology lives in a config file.

Classes come out of `rules.yaml` (animal: siRNA, piRNA) or `rules.plant.yaml` (plant:
phasiRNA, hcsiRNA, siRNA, and optional 24-nt phasiRNA). miRNA is not called by signature; it
is taken from annotation or prediction and carved out of the surrounding locus.

---

## Install

```bash
conda env create -f environment.yml
conda activate srnaprofiler
```

Two tools conda cannot install are handled by the `setup` module:

```bash
python3 srnaprofiler.py setup --check              # what is present, and which modules can run
python3 srnaprofiler.py setup --tools              # miRDeep-P2 (SourceForge tarball + a one-line patch)
python3 srnaprofiler.py setup --data DIR           # Rfam index for the ncRNA filter, ~2 h to build
```

`setup` records the paths it prepared in `~/.config/srnaprofiler/paths.yaml`, so `mapping` and
`mirna` find them without flags.

---

## Order of the steps

```
qc -g genome.fa  →  mapping  →  mirna  →  classify --mirna-loci <BED>  →  classplots
```

`qc -g` is **not optional**: it standardizes FASTA headers and maps IUPAC codes to N. The
miRNA engines reject a genome with either, and NCBI genomes have both.

`classify` has no miRNA class of its own: miRNA loci come from the BED passed to
`--mirna-loci`. A precursor that is not in that BED is classified by its signature like any
other locus, and since miRNAs share their size range and 5' bias with other pathways, it may
be assigned to a neighbouring class. Passing a curated annotation, or the validated
predictions from `mirna`, is what keeps them labelled correctly.

---

## Modules

| module | use it to |
|---|---|
| `setup` | check the toolchain and prepare reference data, once per machine |
| `qc` | trim reads and standardize the genome |
| `mapping` | align to the genome, optionally dropping structural RNA |
| `mirna` | find and validate miRNA loci, so `classify` can carve them out |
| `classify` | call loci, measure their signatures, assign a class |
| `classplots` | look at the result: per-class and per-locus figures |

### Main parameters

**`qc`** — `-p/--prefix`, `-o/--outdir`, `--reads1` (`--reads2` for paired), `-g/--genome`,
`--adapter` (auto-detected if omitted), `--min-len` (default 15), `--qtrim` (20), `--threads`.

**`mapping`** — `--reads`, `-g/--genome` or `--index-prefix`, `--mismatches`,
`--k` (max alignments per read, default 20), `--ncrna-db` (Bowtie index of rRNA/tRNA/snoRNA/
snRNA; drops those reads before clustering), `--threads`. A sorted, indexed BAM is the default
output.

**`mirna`** — `--engine {mirdeep2,mirdeep-p2}` (**no default**: mirdeep2 for animals,
mirdeep-p2 for plants; the wrong one halves precision), `-g/--genome`, `--index-prefix`
(required; without it the engine aligns nothing and reports success with zero predictions),
`--reads`, `--species`, `--mature`/`--hairpin` (miRBase FASTA, optional), `--no-validate`.

**`classify`** — `--sam` (SAM or BAM), `-g/--genome`, `--mincov` (seeding depth in **reads per
million**, default 1.0), `--pad` (island extension in bp), `--config` (`rules.yaml` or
`rules.plant.yaml`), `--mirna-loci` (BED), `--organelle-contigs` (comma-separated contig names
to drop), `--norm`, `--threads`.

**`classplots`** — `--clusters` (the `.clusters.tsv`), `--bam`, `--top-n`, `--classes`,
`--region`, `--logo-format`.

---

## The config

One file per kingdom, and it is meant to be edited. Everything below is per class.

**Locus calling**

| key | meaning |
|---|---|
| `clustering.pad` | how far a coverage island grows on each side before neighbours may merge |
| `extractor.dicermin` / `dicermax` | read lengths kept (18–35) |
| `extractor.strand_cutoff` | above this share on one strand, the locus is called single-stranded |

**Decision**

| key | meaning |
|---|---|
| `decision.review_confidence` | below this softmax confidence, the locus is flagged `review`. Its rigour depends on how many classes exist, so it is not comparable between configs |
| `decision.min_reads` | weighted read count below which a locus is flagged `review` (a flag, not a filter; use `--mincov` to remove noise) |
| `decision.min_margin` | minimum gap to the runner-up; a no-op with only two classes |

**Per class**

| key | meaning |
|---|---|
| `size: {lo, hi}` | the length window that defines the class |
| `w_size` | weight of falling inside that window |
| `w_1u` / `w_5a` | weight of a 5′ U (piRNA, miRNA) or 5′ A (AGO4-loaded hc-siRNA) |
| `w_fractop_mid` | rewards a balanced strand ratio (double-stranded precursor) |
| `w_pingpong` | rewards 10-nt 5′–5′ overlap between strands (secondary piRNA) |
| `w_pingpong_penalty` | penalizes it, for classes where it should be absent |
| `w_phase` | rewards periodic 5′ spacing (phasiRNA, Zucchini-processed piRNA) |
| `phase_x0`, `phase_scale` | midpoint and width of the sigmoid applied to `phase_score` |
| `phase_on: cycles` | score phasing by the run of cycles in register rather than by depth |
| `enabled: false` | switch a class off entirely. This also removes its window from clustering |
| `defining_evidence` | which signatures must be present for the reported `score`; a nested list is an OR (`[[1u, pingpong], size]` = 5′ U **or** ping-pong, **and** size) |

Every weight is multiplied by its own class's size window, so a signature cannot drag a locus
into a class whose length range it does not occupy.

---

## Output

`classify` writes `<prefix>.clusters.tsv`, `<prefix>.clusters.gff3` and `<prefix>.run.json`.

**Counts are 1/NH-weighted.** A read mapping to *N* places contributes 1/*N* to each, so
counts are *not* comparable with tools that assign a multi-read to a single position. Those
are larger by the multiplicity factor, 3.19× in one Drosophila library.

### Main columns

| column | meaning |
|---|---|
| `chrom`, `start`, `end`, `name` | locus position (0-based, half-open) |
| `reads`, `reads_rpm` | weighted abundance, raw and per million |
| `distinct_seqs` | distinct sequences. The sample size of the binomial tests, not an abundance |
| `major_rna`, `dicer_call` | dominant read length, and the length class called from it |
| `frac_in_window` | share of the mass inside the winning class's size window |
| `size_purity` | share of the mass within ±1 nt of the peak. Low means two populations share the locus |
| `dicer_call_2` | a second peak ≥2 nt away carrying ≥20% of the mass |
| `frac_top`, `strand_call` | share on the plus strand, and the resulting call |
| `pct_1u`, `pct_5a`, `pct_10a` | 5′ U, 5′ A, and adenine at position 10 |
| `gc` | GC content |
| `pingpong_z`, `pingpong_frac` | strength of the 10-nt 5′–5′ overlap signature |
| `pingpong_pairs10` | pairs available to measure it. **Filter on this before averaging `pingpong_frac`**, which is 0 below 10 pairs, meaning not measurable rather than absent |
| `phase_score`, `phase_number`, `phase_ratio` | phasing strength, cycles in register, and share of the mass in the dominant register |
| `phase_pvalue` | hypergeometric test, Bonferroni-corrected over the scanned windows |
| `multimap_frac`, `unique_frac` | share of alignments that are multi-mapping, and share of the mass from uniquely mapping reads. `unique_frac` near 0 marks a locus whose signature is borrowed from repeats |
| `cls`, `confidence` | assigned class and softmax confidence |
| `review` | true when the classes compete, the locus is too shallow, or evidence is thin. A flag to inspect, not a rejection |
| `score` | strength of the weakest defining evidence, in [0, 1]; a normalized Wilson lower bound against each signature's null, so it ranks reliability rather than sample size |
| `s_1u`, `s_5a`, `s_fractop_bid`, `s_pingpong`, `s_phase` | the individual p-values, as −log10 |
| `evidence` | which terms fired for the assigned class |

`score` and the `s_*` columns answer different questions: `s_*` say how surprising a signature
is, which grows with depth; `score` says how strong the effect is, which does not. Use `score`
to compare loci.

## License

MIT — see [LICENSE](LICENSE).
