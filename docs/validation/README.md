# Validation tables

Tidy CSVs behind every validation number quoted in `docs/ROADMAP.md`. Regenerate
them instead of copying numbers by hand:

    TAIR_GFF=/path/TAIR10_genomic.gff TAIR_TE=/path/TAIR10_Transposable_Elements.txt \
      python3 extras/collect_validation.py docs/validation/manifest.tsv docs/validation/

`manifest.tsv` records what each run was (mincov, config, miRNA source). The
`.clusters.tsv` files themselves are too large to version and live outside the
repo; the manifest is what makes a run reproducible.

| file | one row per | holds |
|---|---|---|
| `summary.csv` | run | clusters, review rate, class counts |
| `size_distribution.csv` | run x DicerCall | size profile (the kingdom signature) |
| `pic_recall.csv` | Drosophila run | recall over 90 known piRNA clusters, read-mass enrichment |
| `tas_recall.csv` | plant run x TAS locus | class, size, phase score/p-value per landmark |
| `te_enrichment.csv` | plant run x class | transposon overlap vs genome baseline |

## External references

- piRNA clusters (`extras/piC_ref_dm6.bed`): Srivastav, Feschotte & Clark (2024),
  https://doi.org/10.1101/gr.278062.123
- TAS loci: extracted by gene name from the TAIR10 GFF3 (`extras/validate_plant_tas.py`)
- Transposons: TAIR10_Transposable_Elements.txt, https://www.arabidopsis.org
- miRNA annotation: `miRNA` features of the TAIR10 GFF3

## Not yet in the tables

The miRNA-route comparison (ShortStack `MIRNA=Y` at 41/41 precision against
miRDeep2 de novo at 50/94) is still computed ad hoc: it needs the ShortStack and
miRDeep2 outputs in the manifest before it can be regenerated here.
