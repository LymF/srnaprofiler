# Benchmark provenance

All metadata below was looked up on 2026-11-03 against the primary records (ENA/SRA portal and
BioProject XML, GEO series records, the source papers and the databases that serve the truth
files). No reads were downloaded. Run accessions in `manifest.tsv` are SRA/DDBJ run accessions.

Roles: tuning = *Drosophila melanogaster*, *Arabidopsis thaliana*, *Oryza sativa*;
test = *Bombyx mori* (arthropod) and *Zea mays* (plant). All libraries are wild-type,
single-end Illumina small-RNA-seq.

## Libraries

### Drosophila melanogaster (tuning), ovary

| Run | Replicate | Experiment | BioSample | Platform |
|---|---|---|---|---|
| SRR30544372 | 1 | SRX25967759 | SAMN43496128 | Illumina HiSeq 2500 |
| SRR30544371 | 2 | SRX25967760 | SAMN43496128 | Illumina HiSeq 2500 |
| SRR30544370 | 3 | SRX25967761 | SAMN43496128 | Illumina HiSeq 2500 |

- BioProject PRJNA1156618 (SRP530440), "Sequencing RNAs from fruit fly ovaries"; no GEO series.
- Genotype w1118, ovaries. ENA titles: "small-RNA-seq of Drosophila melanogaster: w1118, ovaries,
  biological replicate N". In the paper the same control is called sakura +/+.
- The three runs share one BioSample record; the submitter labels them biological replicates 1-3.
- Source paper: Azlan et al. (2025), eLife 13:RP103828, https://doi.org/10.7554/elife.103828.
  The Data availability statement cites PRJNA1156618. The BioProject also holds sakura mutant,
  heterozygote and rescue libraries; these are not used.

### Arabidopsis thaliana (tuning), inflorescence

| Run | Replicate | Experiment | BioSample | Platform |
|---|---|---|---|---|
| SRR24236126 | 1 | SRX20032294 | SAMN34276913 | NextSeq 500 |
| SRR24236125 | 2 | SRX20032295 | SAMN34276912 | NextSeq 500 |
| SRR24236124 | 3 | SRX20032296 | SAMN34276911 | NextSeq 500 |

- BioProject PRJNA958093; GEO series GSE230228 (GSM7207840, GSM7207841, GSM7207842).
- Genotype Col-0, inflorescence, total RNA, 17-40 nt gel fraction, NEBNext small RNA kit.
  GEO labels: "Wild-type inflorescence RNA, mock-treated, replicate 1/2/3", treatment "mock (CMC)".
  These are the mock arm of a CMC pseudouridine-probing experiment (the CMC reagent is omitted),
  i.e. the untreated control; no CMC-treated, IP or mutant library is used.
- Source paper: Herridge et al. (2025), Nat. Struct. Mol. Biol. 32(2):277-286,
  https://doi.org/10.1038/s41594-024-01392-6 (preprint: https://doi.org/10.1101/2023.05.27.542553).
  The BioProject-to-PubMed link (39242979, 37398006) was verified through NCBI elink.

### Oryza sativa (tuning), anther at early meiosis (ST.2)

| Run | Replicate | Experiment | BioSample | Platform |
|---|---|---|---|---|
| DRR076318 | 1 | DRX070162 | SAMD00066521 | Illumina HiSeq 2500 |
| DRR076319 | 2 | DRX070163 | SAMD00066522 | Illumina HiSeq 2500 |
| DRR076320 | 3 | DRX070164 | SAMD00066523 | Illumina HiSeq 2500 |

- BioProject PRJDB5291 (DRP003992); DDBJ submission DRA005258.
- Cultivar Nipponbare, wild-type ST.2 (early meiosis) anthers, total sRNA-seq; ENA titles
  "total sRNAseq of WT ST.2 (early meiosis) anthers, biological replicate1/2/3".
  The same BioProject holds ST.1/ST.4 wild-type and eat1-4 mutant libraries; not used.
- Source paper: Ono et al. (2018), PLoS Genet. 14(2):e1007238,
  https://doi.org/10.1371/journal.pgen.1007238. Its Data Availability statement cites DRA005258
  for the mRNA-seq, sRNA-seq and MEL1-RIP-seq data.

### Bombyx mori (test), ovary

| Run | Replicate | Experiment | BioSample | Platform |
|---|---|---|---|---|
| DRR377498 | 1 | DRX363399 | SAMD00492477 | Illumina HiSeq 4000 |

- BioProject PRJDB13571; DDBJ submission DRA014190.
- ENA title "small RNA library from ovaries of day 3 of 5th instar larvae (wild type)".
  The strain is not recorded in the ENA metadata.
- Source paper: Kiuchi et al. (2023), PLoS Genet. 19(9):e1010912,
  https://doi.org/10.1371/journal.pgen.1010912. Its Data Availability statement cites DRA014190.
- Only one wild-type ovary small-RNA library exists in this study (the others are BmAgo3 KO,
  or other tissues), so this species has a single replicate. See "Known limitations".

### Zea mays (test), anther (inbred W23)

| Run | Tissue | Replicate | Experiment | BioSample | Platform |
|---|---|---|---|---|---|
| SRR15530410 | anther 1.5 mm | 1 | SRX11828796 | SAMN20863897 | NextSeq 550 |
| SRR15530411 | anther 1.5 mm | 2 | SRX11828797 | SAMN20863898 | NextSeq 550 |
| SRR15530418 | anther 2.0 mm | 1 | SRX11828786 | SAMN20863905 | NextSeq 550 |
| SRR15530419 | anther 2.0 mm | 2 | SRX11828787 | SAMN20863906 | NextSeq 550 |

- BioProject PRJNA756318; GEO series GSE182476 (GSM5529924, GSM5529925, GSM5529932,
  GSM5529933); GEO records PubMed 35451503.
- Wild-type W23 inbred, fixed whole single anthers (GEO overall design), sample names
  "W23_1_5mm, replicate N" and "W23_2_0mm, replicate N". The same series holds laser-capture
  meiocyte, tapetum and outer-cell-layer libraries; not used.
- Source paper: Zhou et al. (2022), New Phytol. 235(2):488-501,
  https://doi.org/10.1111/nph.18167.
- The reads come from W23 but the reference genome is B73 (NAM-5.0); see "Known limitations".

## Genomes

The `genome_url` of each species points to a pinned release, so the file does not change.
Contig names in the Ensembl FASTA files are the bare chromosome names (`2L`, `1`, `Mt`, `scaf_445`);
the truth files were written with the same names (see the per-file notes).

| Species | Assembly (accession) | URL |
|---|---|---|
| D. melanogaster | BDGP6.54 = dm6 / Release 6 (GCA_000001215.4), Ensembl release 116 | https://ftp.ensembl.org/pub/release-116/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz |
| A. thaliana | TAIR10 (GCA_000001735.1), Ensembl Genomes release 63 | https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-63/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz |
| O. sativa | IRGSP-1.0 (GCA_001433935.1), Ensembl Genomes release 63 | https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-63/plants/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz |
| B. mori | BmF-NBRP_1.0 (GCA_030269925.1), NCBI GenBank; FASTA headers are the GenBank accessions BSVE01000001.1-BSVE01000030.1 (checked against the assembly report) | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/269/925/GCA_030269925.1_BmF-NBRP_1.0/GCA_030269925.1_BmF-NBRP_1.0_genomic.fna.gz |
| Z. mays | Zm-B73-REFERENCE-NAM-5.0 (GCA_902167145.1), Ensembl Genomes release 63 | https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-63/plants/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz |

Assembly names, accessions and contig names of the Ensembl assemblies were checked through the
Ensembl REST service (`/info/assembly`), including `Mt`, `AP008247` (rice) and `scaf_445` (maize);
the FASTA files themselves were not downloaded.

## Truth files

All truth files were retrieved on 2026-11-03.

### dmel

`benchmark/truth/dmel/dme_mirgenedb_v3.0.gff3`
- Content: miRNA genes (pre_miRNA and mature miRNA features), 420 features, 2L/2R/3L/3R/X.
- Source: https://mirgenedb.org/gff/dme (MirGeneDB v3.0, curated).
- Assembly: header `genome-build-id: Drosophila_melanogaster_BDGP6` (dm6); contig names match
  the Ensembl BDGP6.54 FASTA. File stored unmodified.
- Paper: Clarke et al. (2025), Nucleic Acids Res. 53(D1):D116-D128.

`benchmark/truth/dmel/dme_piRNAclusterDB2_proTRAC.gff3`
- Content: 196 proTRAC piRNA clusters (class `piRNA`), contigs 2L, 2R, 3L, 3R, 4, X, Y and 6 unplaced
  scaffolds (all present in the Ensembl BDGP6.54 FASTA).
- Source: https://www.pirnaclusterdb.org/data/FASTA/Drosophila_melanogaster.piRNAclusters.gtf
  (piRNAclusterDB 2.0, table row "Drosophila melanogaster, Release_6, 196 clusters"), retrieved 2026-11-03.
- Assembly: Release_6 (GCA_000001215.4 = dm6); contig names already equal the Ensembl names. No liftover.
- Changes: GTF converted to GFF3 (coordinates untouched, 1-based) and `ID=dme_piC_N` added by row order, because the
  attributes column of the GTF is empty.
- Papers: Rosenkranz et al. (2022), piRNAclusterDB 2.0; proTRAC (Rosenkranz & Zischler 2012).
- The tuning seed `truth/seed/piC_ref_dm6.bed` (Srivastav et al. 2024) is a separate source and is listed in the
  manifest for the tuning libraries only.

### atha

`benchmark/truth/atha/ath_mirbase_v23.1.gff3`
- Content: miRNA primary transcripts and mature miRNAs, 894 features, chromosomes 1-5.
- Source: https://mirbase.org/download/CURRENT/genomes/ath.gff3 (header: miRBase v23.1,
  `##date 2026-9-16`).
- Assembly: TAIR10 (GCA_000001735.1), same as `genome_url`.
- Changes: the server returned the GFF3 wrapped in HTML (`<p>`, `<br>`, escaped entities); the
  wrapper was removed and contig names `chrN` were rewritten to `N` to match the Ensembl FASTA.
  A comment line states this in the file. Coordinates are unchanged.
- Paper: Kozomara et al. (2019), Nucleic Acids Res. 47(D1):D155-D162.

`benchmark/truth/atha/ath_sRNAanno_PHAS24.gff3` (79 loci, class `phasiRNA24`) and
`benchmark/truth/atha/ath_sRNAanno_hcsiRNA.gff3` (17482 loci, class `hcsiRNA`)
- Source: sRNAanno (http://www.plantsrnas.org/PHAS24/Arabidopsis_thaliana.PHAS24.gff3 and
  http://www.plantsrnas.org/hc-siRNA/Arabidopsis_thaliana.hc-siRNA.gff3), retrieved 2026-11-03; loci called with
  sRNAminer from public sRNA-seq (Chen et al. 2021). The database gives no version number.
- Assembly: the sRNAanno genome FASTA `Arabidopsis_thaliana.genome.fasta` was checked against the Ensembl TAIR10
  FASTA (28 random 100-bp windows over all 7 sequences, read through HTTP range requests; 0 mismatches), so the
  assembly is TAIR10. No liftover. Contigs `ChrN` rewritten to `N` (`ChrM` to `Mt`); 79/79 and 17482/17482 features kept.
- The sRNAanno 21-PHAS file of Arabidopsis (35 loci) is byte-identical (columns 1-9) to the seed
  `phas21_ref_at_TAIR10.gff3`, so it was not added a second time; the seed stays tuning-only.

### osat

`benchmark/truth/osat/osa_mirbase_v23.1.gff3`
- Content: 1498 miRNA features on chromosomes 1-12, Mt and AP008247.
- Source: https://mirbase.org/download/CURRENT/genomes/osa.gff3 (miRBase v23.1).
- Assembly: IRGSP-1.0 (GCA_001433935.1), same as `genome_url`.
- Changes: HTML wrapper removed; contig `dbj|AP008247|` rewritten to `AP008247`. Coordinates unchanged.
- Paper: Kozomara et al. (2019).

`benchmark/truth/osat/osa_24PHAS_ono2018.bed`
- Content: 104 rice 24-nt PHAS loci (BED6: contig, start 0-based, end, locus ID, 0, strand), from
  S4 Table of Ono et al. (2018) (file `pgen.1007238.s018.xlsx`): section A-1 (93 loci detected in
  that study, EAT1-dependent), A-2 (8 previously reported, EAT1-dependent) and B (3 previously
  reported, EAT1-independent). Sections C (24-PHAS silent in ST.1-ST.4 anthers) and D (non-PHAS
  lincRNA) were left out.
- Source: https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=10.1371/journal.pgen.1007238.s018
- Assembly: the Methods map all reads to "rice IRGSP1.0"; the table itself does not name the
  assembly, so IRGSP-1.0 is inferred from the Methods. Changes: `ChrN` rewritten to `N`; 1-based
  inclusive `from`/`to` converted to BED (start - 1).
- The loci were called from the same anther libraries as the rice libraries above (PRJDB5291), so
  they are valid for tuning only, never as an independent test.
- Paper: Ono et al. (2018).

`benchmark/truth/osat/osa_sRNAanno_PHAS21.gff3` (2595 loci, `phasiRNA`),
`benchmark/truth/osat/osa_sRNAanno_PHAS24.gff3` (804 loci, `phasiRNA24`),
`benchmark/truth/osat/osa_sRNAanno_hcsiRNA.gff3` (194873 loci, `hcsiRNA`, 14.7 MB)
- Source: sRNAanno, http://www.plantsrnas.org/PHAS21/Oryza_sativa.PHAS21.gff3, `.../PHAS24/Oryza_sativa.PHAS24.gff3`,
  `.../hc-siRNA/Oryza_sativa.hc-siRNA.gff3`; retrieved 2026-11-03 (Chen et al. 2021). These sRNAanno PHAS loci are
  a separate source from the Ono et al. 24-PHAS table above and are not merged with it.
- Assembly: 48 random 100-bp windows over Chr1-Chr12 of the sRNAanno `Oryza_sativa.genome.fasta` matched the Ensembl
  IRGSP-1.0 FASTA (0 mismatches; the offsets were computed from the Ensembl lengths, so equal sequence at Chr12 also
  implies equal lengths of Chr1-Chr11). No liftover.
- Changes: `ChrN` rewritten to `N`. The 733 hc-siRNA loci on `ChrUn` (405) and `ChrSy` (328) were dropped because
  these concatenated unplaced/Syngenta pseudo-chromosomes are not in the Ensembl assembly (which carries the
  scaffolds AP008246/7, Syng_TIGR_*, AC*). PHAS21 and PHAS24: 2595/2595 and 804/804 kept; hc-siRNA 194873 of 195606 kept.
- The sRNA-seq libraries that sRNAanno used for rice were not checked, so overlap with the tuning libraries
  (PRJDB5291) cannot be excluded; rice is a tuning species and these files are tuning truth only.

### zmay

`benchmark/truth/zmay/zma_mirbase_v23.1.gff3`
- Content: 501 miRNA features on chromosomes 1-10 and scaffold scaf_445.
- Source: https://mirbase.org/download/CURRENT/genomes/zma.gff3 (miRBase v23.1).
- Assembly: Zm-B73-REFERENCE-NAM-5.0 (GCA_902167145.1), same as `genome_url`.
- Changes: HTML wrapper removed; coordinates and contig names unchanged.
- Paper: Kozomara et al. (2019).

`benchmark/truth/zmay/zma_sRNAanno_PHAS21_NAM5.gff3` (587 loci, `phasiRNA`) and
`benchmark/truth/zmay/zma_sRNAanno_PHAS24_NAM5.gff3` (480 loci, `phasiRNA24`)
- Source: sRNAanno (http://www.plantsrnas.org/PHAS21/Zea_mays.PHAS21.gff3, `.../PHAS24/Zea_mays.PHAS24.gff3`),
  retrieved 2026-11-03 (Chen et al. 2021). Independent of the test libraries (sRNAanno is a separate, pooled
  annotation; no seed material is used).
- Source assembly: B73 RefGen_v4 (AGPv4 = Zm-B73-REFERENCE-GRAMENE-4.0; contigs `1`-`10` and `B73V4_ctgN`). Verified: 40
  random 100-bp windows over chromosomes 1-10 of the sRNAanno `Zea_mays.genome.fasta` equal the MaizeGDB
  Zm-B73-REFERENCE-GRAMENE-4.0 FASTA (0 mismatches).
- Liftover: CrossMap 0.7.3 (`CrossMap gff`), chain
  https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/chain_files/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain
  (md5 492eb7f383b66f2cca2c9349be980eaa, matches `chain_files.md5` of the same directory; assemblies: Jiao et al.
  2017; Hufford et al. 2021). The chain already uses the contig names of the Ensembl FASTA (`1`-`10`, `scaf_N`).
  CrossMap writes strand `+` for features without strand; `.` was restored.
- Counts (in / lifted / dropped): PHAS21 590 / 587 / 3; PHAS24 488 / 480 / 8. All lifted features lie inside their
  contig of the Ensembl NAM-5.0 FASTA index.

`zma_sRNAanno_hcsiRNA_NAM5.gff3` (hcsiRNA, 769615 loci, 59.6 MB) is NOT committed (above the 20 MB limit) and is
therefore not in `truth_index.tsv` or the manifest. Recipe: download
http://www.plantsrnas.org/hc-siRNA/Zea_mays.hc-siRNA.gff3 (790569 loci, v4), run the CrossMap command above
(790569 / 769615 / 20954), restore strand `.`. A built copy (md5 9a5fb9672dac136f6116cf4628ca90fd) is on litrp3 in
`~/srnaprofiler/data/truth_large/`.

### bmor

`benchmark/truth/bmor/bmo_mirbase_v23.1.gff3`
- Content: 1025 miRNA features on the 30 contigs BSVE01000001.1-BSVE01000030.1.
- Source: https://mirbase.org/download/CURRENT/genomes/bmo.gff3 (miRBase v23.1).
- Assembly: BmF-NBRP_1.0 (GCA_030269925.1), same as `genome_url`.
- Changes: HTML wrapper removed; `dbj|BSVE...|` rewritten to the bare GenBank accession; two
  features on `NW_004581752.1` (precursor MI0014305, bmo-mir-3250, a contig that is not in this
  assembly) were dropped.
- Paper: Kozomara et al. (2019).

`benchmark/truth/bmor/bmo_piC_p50ma_ovary_lee2025.gff3`
- Content: 136 ovary piRNA clusters (class `piRNA`, proTRAC, GFF3, contigs BSVE01000001.1...).
- Source: Lee et al. (2025), Sci. Data 12:359, https://doi.org/10.1038/s41597-025-04679-5; file
  `p50ma_Ovary_piCs.gtf` (243 clusters) of the figshare record https://doi.org/10.6084/m9.figshare.26303761 (v2,
  2024-12-21), https://ndownloader.figshare.com/files/51340121.
- Assembly: the clusters are on chromosomes `Bomo_Chr01`...`Bomo_ChrW` of the p50ma assembly, which the paper does not
  tie to an accession. The manifest genome (BmF-NBRP_1.0, strain p50T, GCA_030269925.1) uses the same names
  (`Bomo_ChrNN` = BSVE010000NN.1 in its assembly report) but splits Bomo_Chr24 in `Bomo_Chr24_A/B`. No chain exists, so
  identity was tested on sequence: for the p50ma gene models of the same figshare record (BRAKER, `...longest.renamed.gff3`
  and `...filtered.renamed.cds.fasta`) the CDS was rebuilt from the BmF-NBRP_1.0 FASTA at the same coordinates and compared
  with the published CDS. Chromosomes 01-04, 06, 08-10, 12-18, 21-23, 25-28 (22) gave 100% identical CDS (every gene model
  tested, 290-793 per chromosome); Chr05 (87/727 identical), Chr07 (420/488), Chr11 (250/891), Chr19 (94/571), Chr20 (375/508)
  and ChrW (60/310) are partly shifted; Chr24 does not exist as one sequence in BmF-NBRP_1.0.
- Filter (not a liftover): only clusters on the 22 verified chromosomes were kept (136 of 243); the 107 on Chr05, 07, 11, 19,
  20, 24 and W (14, 11, 9, 6, 10, 32, 25) were dropped. Contigs were renamed to GenBank accessions. Coordinates are unchanged
  (1-based GTF kept as GFF3 start/end).
- Limitation: identity is shown for the gene-bearing sequence only; the piRNA-cluster intervals themselves (largely
  repeats) were not independently compared. The small-RNA data behind these clusters (p50ma ovary) are not the benchmark
  library (DRR377498, Kiuchi et al. 2023).

## Truth expansion (Task 8b): species x class coverage

Every truth file is listed in `benchmark/truth/truth_index.tsv` (path, cls, source, feature_types, aliases) and checked
by `benchmark/tests/test_truth_index.py`. Counts are features loaded by `load_regions` (miRNA: precursor features only).
`seed` = tuning material in `truth/seed/`, listed only for tuning libraries.

| Species (role) | miRNA | piRNA | siRNA | phasiRNA (21) | phasiRNA24 | hcsiRNA |
|---|---|---|---|---|---|---|
| D. melanogaster (tuning) | 140 MirGeneDB | 196 piRNAclusterDB 2.0; 90 seed | none obtained | n/a (animal) | n/a | n/a |
| A. thaliana (tuning) | 296 miRBase | n/a (plant) | none obtained | 35 seed (sRNAanno) + 8 seed TAS | 79 sRNAanno | 17482 sRNAanno |
| O. sativa (tuning) | 508 miRBase | n/a | none obtained | 2595 sRNAanno | 804 sRNAanno; 104 Ono 2018 | 194873 sRNAanno |
| B. mori (test) | 403 miRBase | 136 Lee 2025 (ovary, 22 verified chromosomes) | none obtained | n/a | n/a | n/a |
| Z. mays (test) | 172 miRBase | n/a | none obtained | 587 sRNAanno (lifted) | 480 sRNAanno (lifted) | 769615 sRNAanno (lifted), not committed (60 MB) |

Why nothing usable exists (siRNA): the endo-siRNA / hp-RNA loci of Drosophila (Czech et al. 2008; Ghildiyal et al. 2008;
Okamura et al. 2008) are given in the supplements of those papers, which Europe PMC does not serve as open-access files
and the publisher/PMC pages block automated download; no machine-readable coordinate table on dm6 was found. No B. mori
endo-siRNA locus set was found. Plant non-heterochromatic siRNA loci have no curated public set; the 24-nt hc-siRNA loci
are given as `hcsiRNA`.

### Liftovers and filters

| File | Source assembly | Method | In | Out | Dropped |
|---|---|---|---|---|---|
| zma_sRNAanno_PHAS21_NAM5.gff3 | B73 RefGen_v4 | CrossMap 0.7.3 + MaizeGDB v4-to-NAM-5.0 chain | 590 | 587 | 3 |
| zma_sRNAanno_PHAS24_NAM5.gff3 | B73 RefGen_v4 | same | 488 | 480 | 8 |
| zma_sRNAanno_hcsiRNA_NAM5.gff3 (not committed) | B73 RefGen_v4 | same | 790569 | 769615 | 20954 |
| bmo_piC_p50ma_ovary_lee2025.gff3 | p50ma (Bomo_Chr names) | filter to sequence-verified chromosomes, rename | 243 | 136 | 107 |
| osa_sRNAanno_hcsiRNA.gff3 | IRGSP-1.0 | contig filter (ChrUn, ChrSy) | 195606 | 194873 | 733 |
| all other new files | already on the manifest assembly | contig rename only | all kept | | 0 |

Tools: CrossMap 0.7.3 (conda env `srnabench` on litrp3, bioconda; Zhao et al. 2014), pysam for the sequence checks.
Downloaded genomes and intermediates stay in `~/srnaprofiler/data/` on litrp3 (git-ignored).

## Known limitations

- Bombyx mori has one library only; no replicate exists in the verified study.
- Classes with no usable truth (see the coverage table in "Truth expansion"): siRNA (endo-siRNA/hp-RNA loci) for D. melanogaster
  and B. mori; piRNA for the plants is not applicable; plant `siRNA` (non-hc) has no curated locus set; maize hc-siRNA is
  built but not committed (size).
- piRBase (bigdata.ibp.ac.cn) and the Mainz piRNAclusterDB host could not be reached from the Mac or from litrp3 (DNS/connection
  failures, the piRBase page earlier also returned 403); piRNAclusterDB 2.0 itself was reached at www.pirnaclusterdb.org and has
  no Lepidoptera entry for B. mori (only Plutella xylostella, Trichoplusia ni, Pararge aegeria).
- The Axtell-lab portal (Lunardon et al. 2020) is down for maintenance (notice on the page) and its alternative host returns a
  JavaScript bot check, which was not circumvented.
- Maize reads are W23 but the genome is B73; the B73 NAM-5.0 truth is the miRBase assembly.
- The miRBase file header reads "miRBase v23.1" with date 2026-9-16, as served; it is recorded as is.

## Referências

Azlan, A., Zhu, L., & Fukunaga, R. (2025). Female-germline specific protein Sakura interacts with Otu and is crucial for germline stem cell renewal and differentiation and oogenesis. *eLife*, 13, RP103828. https://doi.org/10.7554/elife.103828

Chen, C., Li, J., Feng, J., Liu, B., Feng, L., Yu, X., Li, G., Zhai, J., Meyers, B. C., & Xia, R. (2021). sRNAanno—a database repository of uniformly annotated small RNAs in plants. *Horticulture Research*, 8(1), 45. https://doi.org/10.1038/s41438-021-00480-8

Clarke, A. W., Høye, E., Hembrom, A. A., Paynter, V. M., Vinther, J., Wyrożemski, Ł., Biryukova, I., Formaggioni, A., Ovchinnikov, V., Herlyn, H., Pierce, A., Wu, C., Aslanzadeh, M., Cheneby, J., Martinez, P., Friedländer, M. R., Hovig, E., Hackenberg, M., Umu, S. U., Johansen, M., Peterson, K. J., & Fromm, B. (2025). MirGeneDB 3.0: improved taxonomic sampling, uniform nomenclature of novel conserved microRNA families and updated covariance models. *Nucleic Acids Research*, 53(D1), D116–D128. https://doi.org/10.1093/nar/gkae1094

Czech, B., Malone, C. D., Zhou, R., Stark, A., Schlingeheyde, C., Dus, M., Perrimon, N., Kellis, M., Wohlschlegel, J. A., Sachidanandam, R., Hannon, G. J., & Brennecke, J. (2008). An endogenous small interfering RNA pathway in Drosophila. *Nature*, 453(7196), 798–802. https://doi.org/10.1038/nature07007

Ghildiyal, M., Seitz, H., Horwich, M. D., Li, C., Du, T., Lee, S., Xu, J., Kittler, E. L. W., Zapp, M. L., Weng, Z., & Zamore, P. D. (2008). Endogenous siRNAs derived from transposons and mRNAs in Drosophila somatic cells. *Science*, 320(5879), 1077–1081. https://doi.org/10.1126/science.1157396

Herridge, R. P., Dolata, J., Migliori, V., de Santis Alves, C., Borges, F., Schorn, A. J., van Ex, F., Lin, A., Bajczyk, M., Parent, J.-S., Leonardi, T., Hendrick, A., Kouzarides, T., & Martienssen, R. A. (2025). Pseudouridine guides germline small RNA transport and epigenetic inheritance. *Nature Structural & Molecular Biology*, 32(2), 277–286. https://doi.org/10.1038/s41594-024-01392-6

Hufford, M. B., Seetharam, A. S., Woodhouse, M. R., Chougule, K. M., Ou, S., Liu, J., Ricci, W. A., Guo, T., Olson, A., Qiu, Y., Della Coletta, R., Tittes, S., Hudson, A. I., Marand, A. P., Wei, S., Lu, Z., Wang, B., Tello-Ruiz, M. K., Piri, R. D., Wang, N., Kim, D. W., Zeng, Y., O'Connor, C. H., Li, X., Gilbert, A. M., Baggs, E., Krasileva, K. V., Portwood, J. L., Cannon, E. K. S., Andorf, C. M., Manchanda, N., Snodgrass, S. J., Hufnagel, D. E., Jiang, Q., Pedersen, S., Syring, M. L., Kudrna, D. A., Llaca, V., Fengler, K., Schmitz, R. J., Ross-Ibarra, J., Yu, J., Gent, J. I., Hirsch, C. N., Ware, D., & Dawe, R. K. (2021). De novo assembly, annotation, and comparative analysis of 26 diverse maize genomes. *Science*, 373(6555), 655–662. https://doi.org/10.1126/science.abg5289

Jiao, Y., Peluso, P., Shi, J., Liang, T., Stitzer, M. C., Wang, B., Campbell, M. S., Stein, J. C., Wei, X., Chin, C.-S., Guill, K., Regulski, M., Kumari, S., Olson, A., Gent, J., Schneider, K. L., Wolfgruber, T. K., May, M. R., Springer, N. M., Antoniou, E., McCombie, W. R., Presting, G. G., McMullen, M., Ross-Ibarra, J., Dawe, R. K., Hastie, A., Rank, D. R., & Ware, D. (2017). Improved maize reference genome with single-molecule technologies. *Nature*, 546(7659), 524–527. https://doi.org/10.1038/nature22971

Kiuchi, T., Shoji, K., Izumi, N., Tomari, Y., & Katsuma, S. (2023). Non-gonadal somatic piRNA pathways ensure sexual differentiation, larval growth, and wing development in silkworms. *PLOS Genetics*, 19(9), e1010912. https://doi.org/10.1371/journal.pgen.1010912

Kozomara, A., Birgaoanu, M., & Griffiths-Jones, S. (2019). miRBase: from microRNA sequences to function. *Nucleic Acids Research*, 47(D1), D155–D162. https://doi.org/10.1093/nar/gky1141

Lee, J., Fujimoto, T., Sahara, K., Toyoda, A., & Shimada, T. (2025). Comprehensive genome annotation of Bombyx mori p50ma strain, a newly developed standard strain. *Scientific Data*, 12(1), 359. https://doi.org/10.1038/s41597-025-04679-5

Okamura, K., Chung, W.-J., Ruby, J. G., Guo, H., Bartel, D. P., & Lai, E. C. (2008). The Drosophila hairpin RNA pathway generates endogenous short interfering RNAs. *Nature*, 453(7196), 803–806. https://doi.org/10.1038/nature07015

Ono, S., Liu, H., Tsuda, K., Fukai, E., Tanaka, K., Sasaki, T., & Nonomura, K.-I. (2018). EAT1 transcription factor, a non-cell-autonomous regulator of pollen production, activates meiotic small RNA biogenesis in rice anther tapetum. *PLOS Genetics*, 14(2), e1007238. https://doi.org/10.1371/journal.pgen.1007238

Rosenkranz, D., & Zischler, H. (2012). proTRAC - a software for probabilistic piRNA cluster detection, visualization and analysis. *BMC Bioinformatics*, 13, 5. https://doi.org/10.1186/1471-2105-13-5

Rosenkranz, D., Zischler, H., & Gebert, D. (2022). piRNAclusterDB 2.0: update and expansion of the piRNA cluster database. *Nucleic Acids Research*, 50(D1), D259–D264. https://doi.org/10.1093/nar/gkab622

Zhao, H., Sun, Z., Wang, J., Huang, H., Kocher, J.-P., & Wang, L. (2014). CrossMap: a versatile tool for coordinate conversion between genome assemblies. *Bioinformatics*, 30(7), 1006–1007. https://doi.org/10.1093/bioinformatics/btt730

Zhou, X., Huang, K., Teng, C., Abdelgawad, A., Batish, M., Meyers, B. C., & Walbot, V. (2022). 24-nt phasiRNAs move from tapetal to meiotic cells in maize anthers. *New Phytologist*, 235(2), 488–501. https://doi.org/10.1111/nph.18167
