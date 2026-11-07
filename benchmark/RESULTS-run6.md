# Run #6: locus floor, read-bounded edges and phased minimum span

Full benchmark run of the locus-floor changes (plan `docs/superpowers/plans/2026-11-06-locus-floor.md`), compared with run #4, the last full run before them. Run date: 2026-11-06/07, server `litrp3`, all 14 libraries of `benchmark/manifest.tsv`. Test libraries (Bombyx, maize) are reported only; no value below was tuned on them.

## 1. Provenance

| run | benchmark commit | pipeline hash | libraries | status |
|---|---|---|---|---|
| #4 | e3a91da (pipeline efb011a) | 86a4a154 | 14 | 14 × exit 0, all `ok` |
| #6 | a5c7615 | b89dedd6 | 14 | 14 × exit 0, all `ok` |

Run #5 (5edb1ee) was superseded before analysis: the final review found that the span floor dropped strongly expressed narrow stacks (fly: 25 fragments with ≥ 100 reads, one with ~13 000 reads over 29 bp), fixed in af75a58 by `locus_span_exempt_reads: 100`. Commit 815acc3, after run #6, changes YAML comments only; because the pipeline hash covers the bytes of `rules*.yaml`, a report generated on 815acc3 would mark run #6 stale although every configuration value is identical. Reports for run #6 are generated on a5c7615.

## 2. What changed between runs #4 and #6

1. **Locus floor.** A region whose 1/NH read mass is below `locus_min_reads` (20), or whose read span is below `locus_min_span` (50 bp) while its mass is below `locus_span_exempt_reads` (100), is not a locus: it is written to `<prefix>.fragments.tsv` and never classified. Regions contained in a supplied miRNA locus are exempt. The floor plays the role of ShortStack's minimum coverage for calling a cluster (Axtell, 2013; Shahid & Axtell, 2014), with values chosen on tuning data.
2. **Read-bounded edges** (`trim_to_reads: true`): locus bounds are the outermost read positions instead of the reads padded by `pad` (500 bp).
3. **Phased minimum span.** phasiRNA must span ≥ 84 bp and phasiRNA24 ≥ 96 bp of reads, i.e. four phase cycles, following the four-phased-read rule of PhaseTank (Guo et al., 2015); a locus that wins the class but is shorter is reported ambiguous (review).
4. **Benchmark accounting.** Tiny-locus counts now use the read span and exclude miRNA loci (from 66a6876 on), and the report adds `fragments` and `fragment_reads` (fragment mass / (fragment + locus mass)).

Values were selected on tuning libraries only (fly SRR30544370, Arabidopsis SRR24236124, rice DRR076318) by a rule fixed before calibration: the strictest setting with no class losing more than 0.01 recall or precision against trim alone, and `fragment_read_frac` ≤ 0.02.

## 3. Tuning libraries (run #4 → run #6)

| species | class | precision | recall | calls | false positives |
|---|---|---|---|---|---|
| Arabidopsis (3 libs) | phasiRNA | 0.223 / 0.237 / 0.197 → 0.247 / 0.288 / 0.233 | 0.744 / 0.674 / 0.744 → 0.744 / 0.698 / 0.767 | 352 / 372 / 389 → 313 / 333 / 344 | 146 / 151 / 167 → 131 / 126 / 145 |
| | hcsiRNA | 0.997 → 0.997–0.998 | 0.718 / 0.714 / 0.718 → 0.704 / 0.700 / 0.704 | ~20 800 → ~20 700 | 49 / 64 / 53 → 44 / 59 / 53 |
| Fly (3 libs) | piRNA | 0.997–0.998 (unchanged) | 0.892 / 0.895 / 0.892 → 0.895 / 0.895 / 0.892 | 4399 / 4585 / 4672 → 4468 / 4646 / 4726 | 2 / 2 / 4 (unchanged) |
| | siRNA | 0.161 / 0.141 / 0.122 → 0.165 / 0.139 / 0.125 | unchanged (0.474 / 0.447 / 0.158) | 753 / 809 / 258 → 731 / 794 / 250 | 182 / 214 / 79 → 177 / 210 / 77 |
| Rice (3 libs) | phasiRNA | 0.984 / 0.981 / 0.984 → 0.989 (all) | 0.882 / 0.917 / 0.921 → 0.881 / 0.916 / 0.919 | 2732 / 2955 / 3009 → 2618 / 2773 / 2824 | 34 / 42 / 36 → 25 / 26 / 26 |
| | phasiRNA24 | 1.000 (unchanged) | 0.171 / 0.130 / 0.131 → 0.171 / 0.130 / 0.132 | ~unchanged | 0 |

miRNA is unchanged in every library (its loci are carved and exempt). ECE moved by at most 0.006. Calls of classes with no truth in a species (hcsiRNA and siRNA in rice; siRNA in Arabidopsis) fell by 8–21%: rice hcsiRNA 570 / 894 / 876 → 509 / 743 / 697.

The only class that lost more than 0.01 is Arabidopsis hcsiRNA recall (−0.014). This is the cost of trimming recorded at calibration: the padding had credited annotated regions with no reads. It was accepted by ruling (2026-11-06) and reverts with `trim_to_reads: false`.

## 4. Test libraries (run #4 → run #6; reported, not tuned on)

| species | class | precision | recall | false positives |
|---|---|---|---|---|
| Bombyx DRR377498 | piRNA | 0.843 (unchanged) | 0.949 (unchanged) | 60 (unchanged) |
| Maize (4 libs) | phasiRNA | 0.959 / 0.941 / 0.967 / 0.967 → 0.975 / 0.956 / 0.976 / 0.977 | 0.731 / 0.814 / 0.736 / 0.697 → 0.734 / 0.809 / 0.739 / 0.697 | 21 / 31 / 17 / 16 → 14 / 23 / 13 / 12 |
| | phasiRNA24 | 1.000 (unchanged) | 0.121 / 0.092 / 0.079 / 0.085 → 0.125 / 0.096 / 0.079 / 0.090 | 0 |

- Maize phasiRNA false positives fell by 25–33% in every library with recall within ±0.005; phasiRNA calls fell 1248 / 1532 / 1210 / 1024 → 1034 / 1443 / 1160 / 905.
- Maize hcsiRNA calls overlapping other-class truth rose slightly in three libraries (268 → 280, 264 → 275, 280 → 291) and were flat in the fourth (186 → 185).
- ECE: Bombyx 0.117 → 0.119; maize within ±0.006.
- 24-PHAS recall in maize stays below 0.13 and remains the main open weakness.

## 5. Fragments and locus size

| species | fragments per library | fragment read share | loci ≤ 50 bp, run #4 → #6 | loci with < 20 reads, run #4 → #6 |
|---|---|---|---|---|
| Arabidopsis | 146–178 | < 0.1% | 376–416 → 44–50 | 64–73 → 0–1 |
| Fly | 83–276 | < 0.1% | 273–375 → 16–29 | 57–163 → 0 |
| Rice | 1556–2089 | 0.3–0.4% | 66–81 → 10–24 | 1505–2039 → 1–7 |
| Bombyx | 121 | < 0.1% | 177 → 17 | 96 → 0 |
| Maize | 181–2359 | 0.0–0.5% | 60–113 → 43–78 | 58–2150 → 0–2 |

The tiny-locus columns are not strictly comparable: run #4 measured the padded window and counted miRNA loci, run #6 measures the read span without miRNA loci. Under one definition the calibration showed the same drop (loci ≤ 50 bp: 212 / 180 / 130 → 4 / 8 / 0 in fly / Arabidopsis / rice). For the same reason, maize loci ≤ 100 bp rose (e.g. 282 → 404): a locus that was 1 kb of padded window can be under 100 bp of reads. The few remaining loci with fewer than 20 reads are expected to be carved miRNA loci, which are exempt from the floor (not checked locus by locus).

## 6. Delimitation

Mean recovered fraction of annotated clusters fell by 0.01–0.06 in every source (e.g. Arabidopsis sRNAanno 0.656 → 0.599; maize sRNAanno 0.639 → 0.587; Bombyx 0.721 → 0.699), and mean fragments per cluster fell slightly (e.g. rice Ono et al. 2018: 1.144 → 1.106). The drop in recovered fraction is the same loss of padding credit as in §3: with read-bounded edges a locus no longer covers annotated bases that have no reads.

## 7. Caveats

- **Per-locus precision depends on fragmentation.** Precision is counted per called locus, so merging fragments shrinks the denominator. Between BASELINE and run #4, Bombyx piRNA went from tp 1670 / fp 71 to tp 321 / fp 60 (12.6 → 3.6 fragments per cluster): precision fell 0.959 → 0.843 while the false positives themselves (piRNA loci on miRBase hairpins) barely changed. Run #6 leaves Bombyx fragmentation unchanged, so its run #4 → #6 comparison is not affected.
- **Density criterion.** The proTRAC piRNA-cluster density test (p ≤ 0.01; Rosenkranz & Zischler, 2012) is now computed over the read span rather than the padded window, which makes it more permissive; its end-to-end effect was checked only on fly piRNA and should be re-checked at the next calibration.
- **Arabidopsis phasiRNA precision remains low (0.23–0.29)** against sRNAanno PHAS21 and the experimental TAS seed, despite the gain.
- **Competitor tools** have not been run; no claim of superiority over other tools follows from this document.

## Referências

Axtell, M. J. (2013). ShortStack: Comprehensive annotation and quantification of small RNA genes. *RNA*, 19(6), 740–751. https://doi.org/10.1261/rna.035279.112

Guo, Q., Qu, X., & Jin, W. (2015). PhaseTank: genome-wide computational identification of phasiRNAs and their regulatory cascades. *Bioinformatics*, 31(2), 284–286. https://doi.org/10.1093/bioinformatics/btu628

Rosenkranz, D., & Zischler, H. (2012). proTRAC – a software for probabilistic piRNA cluster detection, visualization and analysis. *BMC Bioinformatics*, 13, 5. https://doi.org/10.1186/1471-2105-13-5

Shahid, S., & Axtell, M. J. (2014). Identification and annotation of small RNA genes using ShortStack. *Methods*, 67(1), 20–27. https://doi.org/10.1016/j.ymeth.2013.10.004
