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
- No piRNA-cluster table was added (see "Known limitations").

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
- No PHAS/TAS table was added (see "Known limitations").

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

### zmay

`benchmark/truth/zmay/zma_mirbase_v23.1.gff3`
- Content: 501 miRNA features on chromosomes 1-10 and scaffold scaf_445.
- Source: https://mirbase.org/download/CURRENT/genomes/zma.gff3 (miRBase v23.1).
- Assembly: Zm-B73-REFERENCE-NAM-5.0 (GCA_902167145.1), same as `genome_url`.
- Changes: HTML wrapper removed; coordinates and contig names unchanged.
- Paper: Kozomara et al. (2019).
- No PHAS table was added (see "Known limitations").

### bmor

`benchmark/truth/bmor/bmo_mirbase_v23.1.gff3`
- Content: 1025 miRNA features on the 30 contigs BSVE01000001.1-BSVE01000030.1.
- Source: https://mirbase.org/download/CURRENT/genomes/bmo.gff3 (miRBase v23.1).
- Assembly: BmF-NBRP_1.0 (GCA_030269925.1), same as `genome_url`.
- Changes: HTML wrapper removed; `dbj|BSVE...|` rewritten to the bare GenBank accession; two
  features on `NW_004581752.1` (precursor MI0014305, bmo-mir-3250, a contig that is not in this
  assembly) were dropped.
- Paper: Kozomara et al. (2019).

## Known limitations

- Bombyx mori has one library only; no replicate exists in the verified study.
- No curated piRNA-cluster table was added for D. melanogaster or B. mori, and no PHAS/TAS table
  for A. thaliana or Z. mays. The piRBase download page returned 403 and no table on the right
  assembly could be verified. Truth for these classes is therefore missing for those species.
- Maize reads are W23 but the genome is B73; the B73 NAM-5.0 truth is the miRBase assembly.
- The miRBase file header reads "miRBase v23.1" with date 2026-9-16, as served; it is recorded as is.

## Referências

Azlan, A., Zhu, L., & Fukunaga, R. (2025). Female-germline specific protein Sakura interacts with Otu and is crucial for germline stem cell renewal and differentiation and oogenesis. *eLife*, 13, RP103828. https://doi.org/10.7554/elife.103828

Clarke, A. W., Høye, E., Hembrom, A. A., Paynter, V. M., Vinther, J., Wyrożemski, Ł., Biryukova, I., Formaggioni, A., Ovchinnikov, V., Herlyn, H., Pierce, A., Wu, C., Aslanzadeh, M., Cheneby, J., Martinez, P., Friedländer, M. R., Hovig, E., Hackenberg, M., Umu, S. U., Johansen, M., Peterson, K. J., & Fromm, B. (2025). MirGeneDB 3.0: improved taxonomic sampling, uniform nomenclature of novel conserved microRNA families and updated covariance models. *Nucleic Acids Research*, 53(D1), D116–D128. https://doi.org/10.1093/nar/gkae1094

Herridge, R. P., Dolata, J., Migliori, V., de Santis Alves, C., Borges, F., Schorn, A. J., van Ex, F., Lin, A., Bajczyk, M., Parent, J.-S., Leonardi, T., Hendrick, A., Kouzarides, T., & Martienssen, R. A. (2025). Pseudouridine guides germline small RNA transport and epigenetic inheritance. *Nature Structural & Molecular Biology*, 32(2), 277–286. https://doi.org/10.1038/s41594-024-01392-6

Kiuchi, T., Shoji, K., Izumi, N., Tomari, Y., & Katsuma, S. (2023). Non-gonadal somatic piRNA pathways ensure sexual differentiation, larval growth, and wing development in silkworms. *PLOS Genetics*, 19(9), e1010912. https://doi.org/10.1371/journal.pgen.1010912

Kozomara, A., Birgaoanu, M., & Griffiths-Jones, S. (2019). miRBase: from microRNA sequences to function. *Nucleic Acids Research*, 47(D1), D155–D162. https://doi.org/10.1093/nar/gky1141

Ono, S., Liu, H., Tsuda, K., Fukai, E., Tanaka, K., Sasaki, T., & Nonomura, K.-I. (2018). EAT1 transcription factor, a non-cell-autonomous regulator of pollen production, activates meiotic small RNA biogenesis in rice anther tapetum. *PLOS Genetics*, 14(2), e1007238. https://doi.org/10.1371/journal.pgen.1007238

Zhou, X., Huang, K., Teng, C., Abdelgawad, A., Batish, M., Meyers, B. C., & Walbot, V. (2022). 24-nt phasiRNAs move from tapetal to meiotic cells in maize anthers. *New Phytologist*, 235(2), 488–501. https://doi.org/10.1111/nph.18167
