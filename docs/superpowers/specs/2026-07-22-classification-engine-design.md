# Design — Motor de classificação de clusters de pequenos RNAs (srnaprofiler v1)

**Data:** 2026-07-22
**Status:** aprovado no brainstorming, aguardando revisão do spec
**Autor:** Lucas

## 1. Motivação

O srnaprofiler atual separa reads por tamanho e mostra padrões para o **usuário** rotular a classe. A lógica de classificação (`SAMStatistics._classify_loci` em `clusters.py`) é uma heurística de tamanho + balanço de fita com **cortes absolutos** (>50/100/300 reads), sobre **tiles fixos de k-mer** (2000 nt) como locus. Uma auditoria (ver memória `srnaprofiler-fidelity-audit`) mostrou que essa lógica: (a) não é reprodutível entre bibliotecas (cortes absolutos, sem RPM); (b) não mede as assinaturas que definem as classes (viés 1U/10A, ping-pong, phasing) — elas existem em `phasing.py` mas estão **desacopladas** da decisão; (c) usa tiles arbitrários que partem/misturam clusters biológicos.

Este design substitui essa heurística por um **motor de classificação** que **mede** sinais robustos por cluster e **decide** a classe com confiança, seguindo os métodos estabelecidos das ferramentas de referência (ShortStack, piRAT, PhasiHunter, tinyRNA), reimplementados — não copiados.

## 2. Escopo

**Neste spec (v1):**
- Clustering **por cobertura** (ilhas), removendo o tiling.
- Camada **MEDIR** (`FeatureExtractor`) produzindo um `ClusterRecord` por locus.
- Camada **DECIDIR** (`RuleScorer`, modelo de score transparente) → classe + confiança + flag "revisar".
- Classes: **siRNA, piRNA, ambíguo/revisar**. miRNA via **entrada opcional** de loci (interface).
- Novo subcomando `classify`. Ajuste em `qc` (remover tiling).

**Fora deste spec (specs seguintes):**
- Pós-miRDeep (curar `result_*.csv` do miRDeep2 e alimentar `--mirna-loci`).
- `MLDecider` (Random Forest/XGBoost/DL) sobre a matriz de features.
- Clustering por DBSCAN (a fonte de cluster é plugável; entra se a biologia do piRNA exigir).

## 3. Princípio central: separar MEDIR de DECIDIR

Duas camadas com interface bem definida, para que o Decisor seja **plugável** (regras agora, ML depois) sem tocar na medição:

```
SAM/BAM (genoma inteiro) + genoma  (+ regiões BED opcional)  (+ miRNA loci opcional)
        │
        ▼
[clustering por cobertura] ── regiões (default; BED opcional sobrescreve)
        │
        ▼
[FeatureExtractor] ── por locus → ClusterRecord (vetor de features + coords)
        │
        ▼
[Decider (interface)] ── ClusterRecord → classe + confiança + revisar + evidência
     v1: RuleScorer          futuro: MLDecider (mesma interface)
        │
        ▼
[Reporter] ── TSV (= matriz de treino do ML) + GFF + plots + proveniência
```

Duas propriedades garantem o futuro:
- **Fonte de cluster plugável:** default = ilhas de cobertura; um BED sobrescreve; DBSCAN pode entrar depois como outra fonte, sem mudar o resto.
- **Decisor plugável:** trocar `RuleScorer` → `MLDecider` não toca em medição nem no reporter. A matriz de features do v1 **é** o dataset de treino do ML depois.

## 4. Componentes (módulos novos, ex.: pacote `srnaprofiler/classify/`)

| Módulo | Responsabilidade | Depende de |
|---|---|---|
| `clustering.py` | pileup do BAM → ilhas de cobertura → regiões (BED) | pysam |
| `features.py` | regiões + SAM + genoma → `list[ClusterRecord]` (refatora `SAMStatistics`) | pysam, `signals.py` |
| `signals.py` | ping-pong z-score + phase score/p (em cima de `phasing.py`) | `phasing.py` |
| `decider.py` | interface `Decider` + `RuleScorer` | `rules.yaml` |
| `rules.yaml` | config por organismo (janelas de tamanho, pesos, limiares) | — |
| `classify.py` | orquestrador + subcomando `classify` | todos acima |

**Reaproveitados:** `SAMStatistics` (contagem, refatorada para dentro de `features.py`), `phasing.py` (`calculate_offsets` + cobertura), `mirna_filter.py` (exclusão/rotulagem opcional de miRNA por overlap), `utils.py`. O antigo `clusters`/`_classify_loci` fica **intacto** durante a transição (comparação A/B), removido após validação.

## 5. Contrato de dados: `ClusterRecord`

Um registro por locus. É a saída da medição, a entrada do Decisor e — emitido como TSV — a matriz de treino do ML.

| Grupo | Campos |
|---|---|
| Coords | `chrom, start, end, length, name` |
| Abundância | `reads, reads_rpm, distinct_seqs, major_rna, major_rna_reads` |
| Tamanho | `dicer_call, frac_in_window, short, long, counts_por_tamanho[]` |
| Fita | `frac_top, strand_call` |
| Composição | `pct_1U, pct_10A, gc` |
| Biogênese | `pingpong_z, pingpong_pairs10, phase_score, phase_pvalue, phase_number, phase_ratio, overlap19_frac` |
| Repetição | `multimap_frac` |
| Classificação | `class, confidence, review, class_scores{}, evidence` |

## 6. Clustering por cobertura (`clustering.py`)

Método **ilhas de cobertura** (estilo ShortStack) — escolhido por ser a definição biológica consagrada de "locus produtor de sRNA" para múltiplas classes (si/pi/mi), não piRNA-only:
- Posições com cobertura ≥ `mincov`; fundir ilhas cujo gap < `pad`.
- **Agnóstico de fita** — o cluster pode ser bifita (e aí medimos FracTop/ping-pong dentro dele).
- Parâmetros: `--mincov`, `--pad`.
- Pileup por cromossomo via pysam (streaming; não estoura memória em genoma grande).
- A **fonte de cluster é uma interface** — `--regions BED` sobrescreve; DBSCAN pode ser adicionado depois.

## 7. FeatureExtractor (`features.py`) — seguindo os métodos das ferramentas

Refatora o `SAMStatistics`, corrigindo os bugs da auditoria e computando cada sinal pelo método estabelecido do tool de referência:
1. **`-F 256`** — exclui alinhamentos secundários (hoje contamos multimapeador como fita +).
2. **Strand por `flag & 16`; revcomp da fita −** antes de medir composição → **1U/10A corretos**.
3. **RPM** (de `--norm` ou do total mapeado) — elimina os cortes absolutos.
4. **DicerCall + FracTop** por **frações**, janela configurável — como ShortStack.
5. **MajorRNA / DistinctSequences** — foco do locus.
6. Delega a `signals.py`:
   - **ping-pong z-score** — distribuição de overlaps 5′–5′ entre pares de fitas opostas; `z = (C[10] − média(fundo)) / desvio(fundo)`, **fundo = offsets {0–9, 11–20}** (canônico, exclui o 10) — como piRAT/Antoniewski.
   - **phase score** = `phaseRatio × phaseNumber × log(phaseAbun)` e **p hipergeométrico** (Chen 2007), com bins por registro de fase varrendo todos os registros — como PhasiHunter. `phase_length` amarrado ao DicerCall (21 = tasi/phasi; ~26–27 = piRNA head-to-tail).

Interface: `FeatureExtractor(sam, genome, regions=None, lib_size=None, config).extract() -> list[ClusterRecord]`.

## 8. Decider (`RuleScorer`) — modelo de pontuação transparente (Abordagem A)

Para cada `ClusterRecord`, computa **um score por classe** = soma ponderada de **termos de evidência** (cada termo normalizado em [0,1] por uma função de forma). É um modelo linear — o ancestral direto do ML.

**Termos por classe:**
- **siRNA** ↑: DicerCall na janela si (~20–22); FracTop perto de 0.5 (bifita); *bônus* phase score em 21. ↓: ping-pong z alto.
- **piRNA** ↑: DicerCall na janela pi (~24–31); ping-pong z alto; 1U alto; *bônus* multimap alto (TE) e phasing ~26–27.

**Funções de forma (na config, por organismo):** ping-pong `sigmoid((z−z0)/escala)` (z0≈1.65); tamanho = pertinência à janela com bordas suaves; FracTop (siRNA) = gaussiana em 0.5; 1U = a própria fração.

**Confiança e "revisar":**
- `softmax(scores)` → probabilidades; `class = argmax`; `confidence = prob do topo`; `margem = topo − 2º`.
- Flag **revisar** se `confidence < limiar` OU `margem < mínimo` OU score bruto máximo baixo (locus fraco).

**miRNA (opcional):** locus com overlap a uma anotação de miRNA de alta confiança recebe `class=miRNA` (curto-circuito) e não entra no scoring si/pi (via `mirna_filter.py`).

**Evidência:** guardamos `class_scores{}` e os termos que mais pesaram → report mostra *"piRNA (conf 0.86) — ping-pong z=8.2, 1U=0.91, tamanho=27"*.

**Config (`rules.yaml`):** janelas, pesos e limiares por organismo. Começa mínima (2–3 classes) e cresce sem tocar em código. `MLDecider` futuro = mesma interface, pesos aprendidos.

## 9. CLI e fluxo de dados

Impacto em módulos existentes:
- **`qc`**: remove o passo de tiling (`--kmer`/`--no-tiling`); segue padronizando o genoma.
- **`mapping`**: sem mudança de lógica — recebe o genoma inteiro em vez do tileado.
- **`mirna`**: sem mudança — o `classify` só define o contrato `--mirna-loci`.

Novo subcomando:
```
srnaprofiler classify --sam mapped.sam -g genome.fa -p proj -o out/ \
  [--mincov N --pad N]     # clustering por cobertura
  [--regions clusters.bed] # sobrescreve o clustering
  [--mirna-loci mir.gff]   # opcional: miRNA de alta conf p/ rotular+excluir
  [--config rules.yaml]    # default embutido
  [--norm LIBSIZE]         # RPM (senão, total mapeado)
```

**Saídas:**
- `<prefix>.clusters.tsv` — um ClusterRecord por locus (features + classe + conf + revisar + evidência).
- `<prefix>.clusters.gff3` — loci com a classe (genome browser).
- `<prefix>.report/` — plots por classe (tamanho, ping-pong, phasing), reusando `plots.py`/R adaptado.
- `<prefix>.run.json` — proveniência (config, params, versões).

## 10. Tratamento de erros / casos-limite

- Locus vazio / sem reads na janela → `class=NA`, não pontua.
- Uma fita só / ping-pong/phase indefinidos (poucos reads) → termo = 0 (não é erro).
- Sem `--norm` → RPM do total mapeado, com aviso.
- Clusters minúsculos → filtrados no `mincov`.
- Genoma grande → pileup por cromossomo via pysam (streaming).
- Sem config → `rules.yaml` default embarcado.

## 11. Testes

- **Unit por sinal** com SAMs sintéticos de resposta calculável à mão: DicerCall (frac_ok), FracTop, ping-pong z (overlap conhecido → z conhecido), phase score (locus perfeitamente fásico → score alto), 1U/10A.
- **Clustering:** cobertura sintética → ilhas/fronteiras esperadas (mincov/pad).
- **Decider:** ClusterRecords montados à mão → classe+conf esperadas (ping-pong+1U+27nt → piRNA alta conf; 21nt bifita → siRNA).
- **Regressão:** rodar nos dados de teste (`files/finalPIRNA.sam` etc.) e comparar com `_classify_loci` antigo — divergências esperadas e documentadas (bugs corrigidos).
- Ferramentas: `pysam` + `pytest`.

## 12. Referências (métodos reimplementados)

- Axtell, M. J. (2013). ShortStack: Comprehensive annotation and quantification of small RNA genes. *RNA*, 19(6), 740–751. https://doi.org/10.1261/rna.035279.112
- Johnson, N. R., et al. (2016). Improved Placement of Multi-mapping Small RNAs. *G3*, 6(7), 2103–2111. https://doi.org/10.1534/g3.116.030452
- Robak, D., & Ylla, G. (2026). piRAT: piRNA Annotation Tool. *BMC Bioinformatics*. https://doi.org/10.1186/s12859-026-06380-9
- Zhang, C., et al. (2023). PhasiHunter: a robust phased siRNA regulatory cascade mining tool. *Bioinformatics*, 39(11), btad676. https://doi.org/10.1093/bioinformatics/btad676
- Tate, A. J., et al. (2023). tiny-count: hierarchical classification and quantification of small RNA-seq reads. *Bioinformatics Advances*, 3(1), vbad065. https://doi.org/10.1093/bioadv/vbad065
- Antoniewski, C. (2014). Computing siRNA and piRNA overlap signatures. *Methods Mol Biol*, 1173, 135–146. https://doi.org/10.1007/978-1-4939-0931-5_12
- Chen, H.-M., et al. (2007). Bioinformatic prediction and validation of a miRNA-directed tasiRNA cascade in Arabidopsis. *PNAS*, 104(9), 3318–3323. https://doi.org/10.1073/pnas.0611119104
