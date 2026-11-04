# Benchmark do sRNAProfiler — plano de implementação

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Um benchmark reproduzível (dados públicos, gabarito curado e simulado, métricas de classificação e de delimitação, simulação e relatório) que dá à ferramenta atual seu primeiro número de falso positivo por classe.

**Architecture:** Pacote Python `benchmark/` ao lado do pipeline, que chama `srnaprofiler.py` por subprocesso (a ferramenta não é alterada). Cada métrica é um módulo puro, testável com dados sintéticos pequenos. Um manifesto versionado descreve as bibliotecas; os dados baixados ficam em `data/`, fora do git.

**Tech Stack:** Python 3.13, pytest, pysam, pandas, pyyaml (já no env `srnaprofiler`); `iseq` num env separado `benchmark/environment.yml`.

**Spec:** `docs/superpowers/specs/2026-11-03-benchmark-design.md`

## Global Constraints

- Somente dados públicos; as bibliotecas de *Apis mellifera* em `~/Libs_smallRNA_LSV_eric` não entram.
- Entrada da ferramenta no benchmark: somente genoma e reads. Sem anotação de genes, TE ou outras camadas opcionais.
- Espécies de **ajuste**: *Drosophila melanogaster*, *Arabidopsis thaliana*, arroz (*Oryza sativa*). Espécies de **teste**: pelo menos um artrópode além de *Drosophila* e pelo menos uma planta além de *Arabidopsis* e arroz.
- `piC_ref_dm6.bed`, `phas21_ref_at_TAIR10.gff3` e `phas_experimental_at_TAIR10.bed` são material de ajuste e nunca gabarito de teste.
- Sem mutantes: o gabarito é só anotação curada e simulação.
- O passo `mirna` roda em todas as bibliotecas, sem `--mature`: `mirdeep2` em artrópodes, `mirdeep-p2` em plantas.
- Download de SRA com iSeq.
- Limiares reportados: loci minúsculos ≤50 pb e ≤100 pb; poucas reads <20, <50 e <100. Sem meta de aprovação.
- Dados baixados em `data/` (já no `.gitignore`); a lista de acessos é versionada.
- Python 3.13 é fixado por `environment.yml`; sem sudo no litrp3.
- `tests/` e `extras/` estão no `.gitignore`: os testes do benchmark ficam em `benchmark/tests/`.
- Qualquer reajuste de regras após ver o resultado de teste invalida esse resultado.

## Review Focus

1. Biblioteca sem nenhuma read mapeada, ou `clusters.tsv` vazio: as métricas devolvem zeros, não falham (Task 6).
2. Nomes de contig diferentes entre gabarito e genoma (`chr2L` versus `2L`, `NC_...`): hoje daria sobreposição zero em silêncio; deve dar erro listando os contigs sem correspondência (Task 4).
3. Espécie sem nenhum miRNA predito: o passo `mirna` gera um BED vazio e o `classify` segue sem `--mirna-loci`, sem falhar (Task 5).
4. Execução SRA pareada, ou várias execuções por experimento: sRNA é single-end, então rejeitar com mensagem clara (Task 2).
5. Arquivo de regras alterado depois da avaliação de teste: o relatório marca o resultado como inválido (Task 9).

## Estrutura de arquivos

- `benchmark/config.yaml` — limiares de tamanho, de reads e de sobreposição.
- `benchmark/manifest.py`, `benchmark/manifest.tsv` — bibliotecas, papéis e validação de cobertura.
- `benchmark/fetch.py`, `benchmark/environment.yml` — download SRA (iSeq) e genoma.
- `benchmark/truth.py`, `benchmark/truth/seed/` — regiões de gabarito e arquivos recuperados de `ca6ddc2^`.
- `benchmark/run.py` — executa o pipeline e grava metadados.
- `benchmark/metrics_class.py`, `metrics_delim.py` — as duas famílias de métricas.
- `benchmark/simulate.py`, `benchmark/sim_params.yaml` — simulação.
- `benchmark/report.py`, `benchmark/__main__.py` — relatório e CLI.
- `benchmark/tests/` — um arquivo de teste por módulo.

---

### Task 1: Esqueleto, configuração e manifesto

**Files:**
- Create: `benchmark/__init__.py`, `benchmark/__main__.py`, `benchmark/config.yaml`, `benchmark/manifest.py`, `benchmark/tests/test_manifest.py`

**Interfaces:**
- Produces:
  - `@dataclass Library(accession: str, organism: str, role: Literal["tuning","test"], taxon_group: Literal["arthropod","plant"], tissue: str, replicate: int, genome_url: str, truth_sources: list[str])`
  - `load_manifest(path: Path) -> list[Library]`
  - `validate_coverage(libs: list[Library]) -> None` (levanta `ValueError`)
  - `@dataclass Config(tiny_bp: list[int], low_reads: list[int], min_overlap_frac: float)` e `load_config(path: Path) -> Config`
  - CLI `python -m benchmark check-manifest`

- [ ] **Step 1: Escrever os testes que falham** em `benchmark/tests/test_manifest.py`: `test_config_defaults_match_spec` (`tiny_bp == [50, 100]`, `low_reads == [20, 50, 100]`), `test_load_manifest_roundtrip`, `test_validate_coverage_rejects_missing_test_species` (sem artrópode ou planta de teste → `ValueError`), `test_rice_cannot_be_test_role` (arroz com `role="test"` → `ValueError`), `test_tuning_species_required` (faltar *Drosophila*, *Arabidopsis* ou arroz em ajuste → `ValueError`).
- [ ] **Step 2: Rodar e ver falhar:** `python -m pytest benchmark/tests/test_manifest.py -v` → FAIL (módulo inexistente).
- [ ] **Step 3: Implementar** `manifest.py` e `config.yaml`. Valor inicial a conferir (decisão do assistente): `min_overlap_frac: 0.5`.
- [ ] **Step 4: Rodar e ver passar:** mesmo comando → PASS.
- [ ] **Step 5: Commit:** `git add benchmark && git commit -m "feat(benchmark): manifest, config and coverage validation"`

### Task 2: Download de bibliotecas e genomas

**Files:**
- Create: `benchmark/fetch.py`, `benchmark/environment.yml`, `benchmark/tests/test_fetch.py`

**Interfaces:**
- Consumes: `Library` (Task 1).
- Produces: `fetch_library(lib: Library, dest: Path) -> Path` (FASTQ `.fastq.gz` single-end), `fetch_genome(lib: Library, dest: Path) -> Path` (FASTA).

- [ ] **Step 1: Testes que falham** com `subprocess.run` substituído por um falso: `test_fetch_library_skips_existing` (arquivo já presente → nenhuma chamada), `test_fetch_library_calls_iseq_with_accession`, `test_paired_run_rejected` (a execução gera `_1` e `_2` → `ValueError` mencionando "paired").
- [ ] **Step 2: Rodar e ver falhar.**
- [ ] **Step 3: Implementar** chamando `iseq` (conferir as flags com `iseq --help` antes; pedir saída gzip). Idempotente. `benchmark/environment.yml` instala `iseq` (bioconda) num env separado, sem tocar o `environment.yml` principal.
- [ ] **Step 4: Rodar e ver passar.**
- [ ] **Step 5: Commit:** `git commit -m "feat(benchmark): idempotent SRA and genome fetch"`

### Task 3: Curar o manifesto

**Files:**
- Create: `benchmark/manifest.tsv`, `benchmark/PROVENANCE.md`

**Interfaces:** Produces `manifest.tsv` no formato de `Library`.

- [ ] **Step 1:** Para cada espécie exigida pelo spec, escolher bibliotecas selvagens com 2 a 3 réplicas, começando por *Drosophila*. Espécies de teste: um artrópode (*Aedes aegypti* ou *Bombyx mori*) e milho. Usar busca na SRA e conferir cada acesso contra o artigo de origem.
- [ ] **Step 2:** Registrar em `PROVENANCE.md`, por biblioteca, o artigo de origem (autores, título, periódico, ano, DOI) e a fonte do gabarito. Nenhum acesso entra sem artigo conferido.
- [ ] **Step 3: Verificar:** `python -m benchmark check-manifest` → `OK: N libraries, coverage satisfied`.
- [ ] **Step 4: Commit:** `git commit -m "data(benchmark): curated manifest with provenance"`

### Task 4: Gabarito curado

**Files:**
- Create: `benchmark/truth.py`, `benchmark/truth/seed/` (arquivos recuperados), `benchmark/tests/test_truth.py`

**Interfaces:**
- Consumes: `Library` (Task 1).
- Produces: `@dataclass Region(chrom: str, start: int, end: int, cls: str, source: str, name: str)`, `load_regions(path: Path, cls: str, source: str) -> list[Region]` (BED ou GFF3), `check_contigs(regions: list[Region], fai: Path) -> None`, `guard_seed_use(lib: Library, path: Path) -> None`.

- [ ] **Step 1: Testes que falham:** `test_load_bed_regions`, `test_load_gff3_regions`, `test_contig_mismatch_raises` (região em `chr2L` com genoma em `2L` → `ValueError` listando o contig), `test_seed_files_rejected_for_test_species` (`role="test"` e caminho em `truth/seed/` → `ValueError`).
- [ ] **Step 2: Rodar e ver falhar.**
- [ ] **Step 3: Recuperar o material do histórico:** `git show ca6ddc2^:extras/piC_ref_dm6.bed`, `phas21_ref_at_TAIR10.gff3` e `phas_experimental_at_TAIR10.bed` para `benchmark/truth/seed/`, mais um `README` de uma linha dizendo "somente ajuste". Implementar `truth.py`.
- [ ] **Step 4: Rodar e ver passar.**
- [ ] **Step 5: Commit:** `git commit -m "feat(benchmark): curated truth regions and seed guard"`

### Task 5: Execução do pipeline por biblioteca

**Files:**
- Create: `benchmark/run.py`, `benchmark/tests/test_run.py`

**Interfaces:**
- Consumes: `Library` (Task 1).
- Produces: `@dataclass RunResult(clusters_tsv: Path, bam: Path, meta: Path)` e `run_library(lib: Library, reads: Path, genome: Path, workdir: Path, rules: Path) -> RunResult`. `meta` é `benchmark_run.json` com `tool_git_sha`, `rules_sha256`, `accession`, `mirna_engine: str`.

- [ ] **Step 1: Testes que falham** com subprocesso falso: `test_command_sequence` (ordem `qc`, `mapping`, `mirna`, `classify`; flags `-p`, `-o`, `--reads1`, `-g`, `--index-prefix`, `--sam`, `--config`), `test_mirna_runs_for_every_library_without_mature` (sem `--mature`; `--species` = `lib.organism`), `test_engine_by_taxon_group` (artrópode→`mirdeep2`, planta→`mirdeep-p2`), `test_empty_mirna_bed_does_not_fail`, `test_meta_records_rules_sha256`.
- [ ] **Step 2: Rodar e ver falhar.**
- [ ] **Step 3: Implementar.** Nomes de saída seguem o README: `{prefix}.trimmed.R1.fastq.gz`, `{prefix}.genome.std.fa`, `{prefix}.bowtie.bam`, `{prefix}.clusters.tsv`. O passo `mirna` sempre roda, com `--index-prefix` do `mapping`; o BED de loci de miRNA resultante vai para `classify --mirna-loci`.
- [ ] **Step 4: Rodar e ver passar.**
- [ ] **Step 5: Commit:** `git commit -m "feat(benchmark): per-library pipeline runner with provenance"`

### Task 6: Métricas de classificação

**Files:**
- Create: `benchmark/metrics_class.py`, `benchmark/tests/test_metrics_class.py`

**Interfaces:**
- Consumes: `Region` (Task 4), `Config` (Task 1); `loci` são os dicts lidos de `clusters.tsv` (colunas `chrom,start,end,cls,confidence,reads,...`).
- Produces: `match_loci(loci: list[dict], regions: list[Region], min_overlap_frac: float) -> list[tuple[dict, Region | None]]` e `classification_metrics(loci: list[dict], regions: list[Region], cfg: Config) -> dict` com chaves `per_class` (`precision`, `recall`, `fdr`), `ambiguous_coverage`, `unknown` (`count`, tamanho e reads em quantis), `ece`.

- [ ] **Step 1: Testes que falham:** `test_precision_recall_per_class`, `test_ambiguous_is_abstention_not_error`, `test_unannotated_loci_reported_as_unknown_not_fp`, `test_ece_perfect_calibration_is_zero`, `test_empty_loci_returns_zero_metrics`.
- [ ] **Step 2: Rodar e ver falhar.**
- [ ] **Step 3: Implementar.** `ambiguous` sai do denominador da precisão e entra em `ambiguous_coverage`; locus sem região correspondente vai para `unknown`.
- [ ] **Step 4: Rodar e ver passar.**
- [ ] **Step 5: Commit:** `git commit -m "feat(benchmark): classification metrics"`

### Task 7: Métricas de delimitação

**Files:**
- Create: `benchmark/metrics_delim.py`, `benchmark/tests/test_metrics_delim.py`

**Interfaces:**
- Consumes: `Region` (Task 4), `Config` (Task 1).
- Produces: `delimitation_metrics(loci: list[dict], clusters: list[Region], cfg: Config) -> dict` com `per_cluster` (`recovered_fraction`, `n_fragments`), `merged_distinct`, `tiny_counts` (por limiar de pb) e `low_read_counts` (por limiar de reads).

- [ ] **Step 1: Testes que falham:** `test_cluster_split_into_three_reports_three_fragments`, `test_recovered_fraction_is_union_of_overlaps`, `test_locus_spanning_two_clusters_counts_as_merged`, `test_tiny_counts_use_config_thresholds`.
- [ ] **Step 2: Rodar e ver falhar.**
- [ ] **Step 3: Implementar.** `tiny_counts` e `low_read_counts` consideram todos os loci, não só os que casam com o gabarito, para expor o problema dos loci minúsculos.
- [ ] **Step 4: Rodar e ver passar.**
- [ ] **Step 5: Commit:** `git commit -m "feat(benchmark): delimitation metrics"`

### Task 8: Simulação

**Files:**
- Create: `benchmark/simulate.py`, `benchmark/sim_params.yaml`, `benchmark/tests/test_simulate.py`

**Interfaces:**
- Consumes: `Region` (Task 4).
- Produces: `@dataclass SimSpec(classes: dict[str, int], seed: int)`, `@dataclass SimResult(fastq: Path, truth_bed: Path)` e `simulate(spec: SimSpec, background_fastq: Path, genome: Path, out_dir: Path) -> SimResult`.

- [ ] **Step 1: Testes que falham:** `test_simulate_is_deterministic_for_same_seed`, `test_truth_coordinates_inside_chromosome`, `test_pirna_reads_have_10nt_overlap_on_opposite_strands`, `test_simulate_does_not_import_rules` (o módulo não importa `classify` nem lê `rules*.yaml`).
- [ ] **Step 2: Rodar e ver falhar.**
- [ ] **Step 3: Implementar.** Injeta loci sintéticos em uma biblioteca real usada como fundo. Cada parâmetro de `sim_params.yaml` traz o DOI da fonte na literatura; nenhum valor vem das regras da ferramenta.
- [ ] **Step 4: Rodar e ver passar.**
- [ ] **Step 5: Commit:** `git commit -m "feat(benchmark): spike-in simulation independent of rules"`

### Task 9: Relatório e CLI

**Files:**
- Create: `benchmark/report.py`, `benchmark/tests/test_report.py`
- Modify: `benchmark/__main__.py`

**Interfaces:**
- Consumes: `RunResult` (Task 5), as duas funções de métricas (Tasks 6 e 7), `Library`.
- Produces: `@dataclass LibraryResult(lib: Library, metrics: dict, meta: dict)`, `build_report(results: list[LibraryResult], rules: Path) -> str` (markdown) e CLI `python -m benchmark run --library ACC` e `python -m benchmark report`.

- [ ] **Step 1: Testes que falham:** `test_report_keeps_tuning_and_test_separate` (seções distintas "Ajuste" e "Teste", nunca somadas), `test_report_flags_rules_changed_after_test_run` (`rules_sha256` gravado ≠ hash atual → resultado de teste marcado INVÁLIDO).
- [ ] **Step 2: Rodar e ver falhar.**
- [ ] **Step 3: Implementar** o relatório e os subcomandos `run` e `report`. Só resultados de espécies de teste aparecem sob a seção que sustenta a afirmação de que a ferramenta classifica melhor.
- [ ] **Step 4: Rodar a suíte inteira:** `python -m pytest benchmark/tests -v` → todos PASS.
- [ ] **Step 5: Commit:** `git commit -m "feat(benchmark): role-separated report and CLI"`

### Task 10: Linha de base em *Drosophila*

**Files:**
- Create: `benchmark/BASELINE.md`
- Modify: `.gitignore` (adicionar `benchmark/results/`)

- [ ] **Step 1:** No litrp3, baixar as bibliotecas de *Drosophila* do manifesto com `python -m benchmark run --library ACC` para cada uma.
- [ ] **Step 2: Verificar:** existe `benchmark/results/<run_id>/metrics.json` por biblioteca, com as chaves `per_class`, `unknown`, `per_cluster`, `tiny_counts` e `low_read_counts`.
- [ ] **Step 3:** Gerar `python -m benchmark report` e copiar os números de falso positivo por classe e de delimitação para `BASELINE.md`, rotulados "Ajuste".
- [ ] **Step 4: Commit:** `git add benchmark/BASELINE.md .gitignore && git commit -m "docs(benchmark): Drosophila baseline"`

---

## Self-review

- **Cobertura do spec:** táxons e papéis (Tasks 1, 3); anotações curadas e simulação (4, 8); falso positivo de classe e de delimitação (6, 7); `ambiguous` e desconhecidos (6); loci minúsculos reportados (7); miRNA em todas as bibliotecas (5); circularidade e invalidação por regras (4, 9); dados e reprodutibilidade (2, 3, 5); primeiro número de falso positivo (10). Mutantes retirados por decisão do usuário. **Fora do plano por decisão do spec:** melhorar o classificador, ML, outras ferramentas e camadas opcionais.
- **Consistência de tipos:** `Library`, `Region`, `Config`, `RunResult` e `LibraryResult` têm os mesmos nomes e campos em todas as tasks que os consomem.
- **Proporção:** o plano tem uma task por módulo, com assinaturas e nomes de teste, sem corpos de função.
- **Pontos que dependem de fonte externa e não foram inventados:** acessos SRA, DOIs dos artigos e valores de `sim_params.yaml` são preenchidos nas Tasks 3 e 8, a partir da literatura.
