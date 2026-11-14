# `te-library`: build a TE library for any genome — design

## Intent

`annotate-loci` needs a transposon library. For model clades the Dfam curated library
is enough, but measured on Arabidopsis it is not: Dfam 4.0 curated has 9 families for
*A. thaliana* and 0 lineage-specific, and RepeatMasker with it found 0 TEs in 22584 loci
(recall 0.2% against TAIR10). RepeatModeler2 on the genome (3h49, 32 threads) gave 631
families and, used in `annotate-loci`, 68.8% recall and 83.2% precision against TAIR10
(130 s for the loci).

The library has to be built per genome, and doing it by hand is error prone. `te-library`
runs the recipe and removes the known contaminants. It is optional and separate:
`classify` and `report` never call it. (This reverses the non-goal "not wrapping
RepeatModeler" of the annotation design; decided with the user, 2026-11-14.)

## Measured contamination (Arabidopsis, hcsiRNA loci)

Loci called TE by the RM2 library only: 1624, of which 20.8% lie >=50% on a gene exon,
against 0.9% of the loci where RM2 and TAIR10 agree. So multicopy gene families end up
in the library as "TE". A filter against gene sequence is needed.

## Interface

```
srnaprofiler te-library -g genome.fa -o OUTDIR -p PREFIX [--threads N] [--resume]
    [--organelle-contigs Pt,Mt] [--genes GENES.gff3|gtf[.gz]] [--exclude-fasta FASTA ...]
    [--dfam-species "Latin name"] [--families EXISTING.fa]
```

Output: `OUTDIR/PREFIX-families.fa` (the final library, ready for
`annotate-loci --te-library`), `OUTDIR/PREFIX.te-library.report.tsv` (one row per
input family: name, length, kept/dropped, reason, matched sequence) and the logs.

## Steps (in order)

1. **Nuclear genome.** Write `OUTDIR/PREFIX.nuclear.fa` without `--organelle-contigs`
   (same flag as `classify`). Skipped (and the genome used as is) when none is given.
2. **RepeatModeler2.** `BuildDatabase -name PREFIX nuclear.fa`, then
   `RepeatModeler -database PREFIX -threads N -LTRStruct`, run inside OUTDIR. With
   `--resume` and an existing `RM_*` directory from an interrupted run, pass
   `-recoverDir`. The RepeatModeler2 library is `RM_*/consensi.fa.classified`
   (RepeatModeler2 also copies it to `PREFIX-families.fa` in the working directory when
   it finishes; use that, else the classified consensi).
   `--families EXISTING.fa` skips this step and starts from that library (a library
   from elsewhere, or a previous run).
3. **Filter contaminants**, each family tested with `minimap2 -x asm20` against the
   exclusion sequences. A family is **dropped** when one alignment gives identity >= 0.80
   and covers >= 50% of the family length. Exclusion sequences:
   - `--genes`: the exons of the non-TE genes of that annotation, extracted from the
     genome (use `annotation.py`: the gene-layer rules already drop `transposable_element*`
     rows and their children, and keep exon-like rows only; contigs are those of the FASTA);
   - `--exclude-fasta` (repeatable): any FASTA, for example Rfam rRNA/tRNA/snoRNA
     sequences, or a protein-free CDS set of the user.
   Without either option no filtering happens and the report says so. Families classed
   `rRNA`, `tRNA`, `snRNA`, `scRNA` or `srpRNA` (the part of the name after `#` and before
   `/`) are dropped outright, filter options or not: they are not transposons. `Satellite`
   is kept.
4. **Dfam merge (optional).** `--dfam-species`: append the curated Dfam consensi for
   that taxon and its ancestors (`famdb.py families --curated -a -ad -f fasta_name TAXON`,
   found through `famdb.py` on PATH), with a `dfam|` prefix on the names and no
   duplicate names. Missing `famdb.py` or an empty result is a WARN, not an error.
5. **Write** `PREFIX-families.fa` and the report; log counts per class and per reason.

## Constraints

- Opt-in and isolated: no change to `classify`/`report`; no new hard dependency.
  External tools (`BuildDatabase`, `RepeatModeler`, `minimap2`, `famdb.py`) are looked up
  on PATH; a missing required one is an ERROR naming it and the install hint; `setup
  --check` lists them under `te-library` as optional (follow the existing `annotate-loci`
  entry in `setup_module.py`).
- Large runs: never hold the library in memory twice; stream FASTA.
- RepeatModeler2 is slow (hours): print what is running, where its log is, and how to
  resume. Do not capture its output in memory.
- Refuse to overwrite an existing `PREFIX-families.fa` unless it came from this tool's
  incomplete run (follow what `annotate_loci.py` does with an existing output).
- Tests use tiny FASTA fixtures and a fake `minimap2`/`RepeatModeler` on PATH; no test
  runs the real tools.

## Verification (on the server, not in unit tests)

On Arabidopsis, filter the existing RM2 library (`--families ~/annot-test8/ath-families.fa
--genes <Ensembl TAIR10.63 gff3>`), rerun `annotate-loci --te-library`, and report the
share of "RM2-only" hcsiRNA loci on gene exons (was 20.8%) and the recall against TAIR10
(was 68.8%): the filter must lower the first without losing much of the second.
