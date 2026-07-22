# srnaprofiler — Roadmap / Próximos passos

**Autor:** Lucas
**Atualizado:** 2026-07-22

Documento vivo para alinhar sessões. Estado atual + backlog priorizado. Detalhes de COMO a
pipeline funciona (todas as fórmulas) estão na memória `srnaprofiler-pipeline-reference`.
Design/planos: `docs/superpowers/specs/2026-07-22-classification-engine-design.md` e
`docs/superpowers/plans/2026-07-22-classification-engine.md`.

---

## Estado atual (PRONTO)

Motor **classify** v1 completo, na branch `master` (local — falta `git push`). **37 testes** passando.
- Clustering por **ilhas de cobertura** (removeu o tiling do `qc`).
- `FeatureExtractor` por locus → `ClusterRecord`: DicerCall, FracTop, 1U/10A/GC, ping-pong (z guardado+capado e `pingpong_frac`), phasing (phaseScore + hipergeométrico), multimap.
- **Ponderação 1/NH** (via tag `XM` do bowtie1, detectado pelo `@PG`) — funciona: nos dados reais 96% são multimapeadores.
- `RuleScorer` (score transparente): classes **siRNA/piRNA**; "ambíguo" = flag `review`. miRNA só por `--mirna-loci` (overlap).
- Threaded **por cromossomo** (`--threads`), lê de BAM ordenado+indexado.
- Reporter: `.clusters.tsv` (= matriz de treino do ML), `.gff3`, `.run.json`.
- Rodado ponta a ponta em dados reais: **DRR351382** (ovário Drosophila WT). Resultado biológico coerente (piRNA dominante, 24–27nt, 1U alto, ping-pong nos loci bifita).

**Aprendizados operacionais** (na memória): `/mnt/d` é drvfs (Windows) — `samtools sort` quebra ("truncated file") escrevendo temporários lá; rodar sort/index em disco LOCAL (`~`), só TSV/GFF finais vão pro D:.

---

## Próximos passos

### Tier 1 — imediato (pequeno, alto valor)

1. **`--tmpdir` no sort (mapping + classify).**
   - Porquê: hoje `samtools sort` grava temporários na pasta de saída; se for `/mnt/d` (drvfs) → "truncated file". Força pré-ordenar na mão.
   - Onde: `samio.ensure_sorted_indexed_bam` (passar `-T <tmpdir>`), `mapping.convert_sam_to_bam`, e um arg `--tmpdir` (default = disco local, ex. `tempfile.gettempdir()`) nos subcomandos `mapping`/`classify`.
   - Pronto quando: `mapping --create-bam` e `classify --sam x.sam` rodam com saída em `/mnt/d` sem erro.

2. **Calibrar `mincov`/`min_reads` para a escala PONDERADA.**
   - Porquê: com 1/NH, contagens honestas caem → **79% dos clusters viram `review`**. Muitos são ilhas de baixa abundância real.
   - Onde: defaults em `srnaprofiler.py` (subparser classify) e `rules.yaml` (`decision.min_reads`); talvez subir `mincov` default (2→5). Decisão empírica: varrer valores e ver a taxa de review.
   - Pronto quando: taxa de review razoável (<~40%) sem perder loci reais de piRNA/siRNA.

3. **`git push` da master.** Local está à frente do `origin` (remoto foi limpo com filter-repo). Publicar.

### Tier 2 — validação biológica + refinamento

4. **Validar as chamadas contra verdade conhecida.**
   - Porquê: confirmar que nossos loci de piRNA/siRNA batem com a biologia.
   - Como: cruzar os loci `piRNA` (GFF) com clusters de piRNA conhecidos de Drosophila (piRNAclusterDB; flamenco em `NT_037436`/X), e os `siRNA` 21nt bifita com loci de endo-siRNA/TE. Usar `bedtools intersect`.
   - Pronto quando: sobreposição alta dos clusters conhecidos e sanidade das assinaturas (1U em piRNA, ping-pong em loci germinativos).

5. **Endurecer a guarda de `pingpong_frac` em loci minúsculos.**
   - Porquê: `pingpong_frac=1.0` espúrio em loci de 11–42 reads (poucos pares que por acaso caem no offset 10).
   - Onde: `features.build_record` (guarda de mínimo de pares — hoje `total_pairs≥10`, subir/tornar depth-aware).

6. **RPM consistente com a ponderação.**
   - Porquê: `lib_size` = `total_mapped` (linhas de alinhamento ≈ 4× reads únicos com `-k5`); RPM fica sub-escalado (consistente entre loci, mas escala absoluta off).
   - Onde: `classify.run_classify` — `lib_size` = nº de reads ÚNICOS (Σ 1/NH). Baixa prioridade (RPM é só reportado, não usado na decisão).

7. **Refinar/validar o phasing hipergeométrico.** Parametrização atual é aproximação v1 (janela 10 ciclos). Validar com locus fásico conhecido (tasiRNA) ou marcar coluna como experimental.

### Tier 3 — features maiores (specs próprios)

8. **Módulo pós-miRDeep (rota miRNA).**
   - Escopo: curar `result_*.csv` do miRDeep2 (score, star>0, randfold p, descartar rfam alert; MFEI/AMFE via RNAfold), isomiRs, e **alimentar `--mirna-loci`** do classify (reusa `mirna_filter.py` por overlap). miRDeep2 continua o motor; nós curamos/integramos.
   - Ver decisões em `srnaprofiler-tools-research` (memória).

9. **Plots por classe no reporter.** Distribuição de tamanho, ping-pong, phasing por classe (reusar `plots.py`/R). Estava no spec, adiado no v1.

10. **`MLDecider` (fase 2 de ML).** RandomForest/XGBoost sobre o `.clusters.tsv` (que já é a matriz de features); interpretar com SHAP. Precisa de loci rotulados (da validação #4). Mesma interface do `RuleScorer` (plugável) — troca sem tocar em medição.

11. **Clustering por DBSCAN (fonte plugável alternativa).** Se as ilhas de cobertura fragmentarem clusters grandes/esparsos de piRNA. `clustering.py` já é a fonte; adicionar DBSCAN por parâmetro.

---

## Ponteiros internos
- **Como funciona + fórmulas:** memória `srnaprofiler-pipeline-reference`.
- **Ferramentas de referência estudadas** (em `ref_tools/`, gitignored) e decisões: memória `srnaprofiler-tools-research`.
- **Auditoria da heurística Perl antiga:** memória `srnaprofiler-fidelity-audit`.
- **Dados de teste:** DRR351382 (ENA), BAM local ordenado em `~/dros.sorted.bam`; resultados em `/mnt/d/srnaprofiler-teste/results` e `~/classify_out`.
- **Convenções:** código/comentários em inglês; autor sempre "Lucas" (nunca IA, nem em commits).
