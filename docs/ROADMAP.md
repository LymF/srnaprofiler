# srnaprofiler — Roadmap / Próximos passos

**Autor:** Lucas
**Atualizado:** 2026-08-01

Documento vivo para alinhar sessões. Estado atual + backlog priorizado. Detalhes de COMO a
pipeline funciona (todas as fórmulas) estão na memória `srnaprofiler-pipeline-reference`.
Design/planos: `docs/superpowers/specs/2026-07-22-classification-engine-design.md` e
`docs/superpowers/plans/2026-07-22-classification-engine.md`.

---

## Estado atual (PRONTO)

Motor **classify** v1 completo + subcomando **`classplots`**, na branch `master` (local — falta `git push`). **73 testes** passando.
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

> Nota: o problema do drvfs (`/mnt/d`) NÃO será corrigido na tool (decidido). Workaround: ordenar/indexar
> em disco local (`~`), só os TSV/GFF finais vão pro D:. Não é item de backlog.

### Tier 1 — imediato (pequeno, alto valor)

1. ~~**Calibrar `mincov`/`min_reads` para a escala PONDERADA.**~~ **FEITO (2026-08-01).**
   - Diagnóstico: dos 79% em `review`, **75,9% caíam por `reads < 10`** e só 3,4% por confiança sozinha —
     o gargalo era o clustering, não o scorer (mediana de 3 reads ponderadas por cluster com `mincov=2`).
   - Varredura em DRR351382 (`mincov` 2→50, tabela no README): loci mantidos e loci ≥50 reads têm pico em
     20–30 com massa de reads retida praticamente plana; os 1000 loci mais abundantes de `mincov=2` são
     recuperados 100% em todos os cortes. Subir `min_reads` é a alavanca errada (só rerrotula).
   - Mudado: `--mincov` default 2→**10** (review 79%→56%; para libs ≥30M alinhamentos usar 20–30, que
     bate a meta <40%). `min_margin` documentado como no-op enquanto houver só duas classes
     (com softmax de 2 vias, `conf<0,6` ⟺ `margin<ln(1,5)≈0,405`).

2. **`git push` da master.** Local está à frente do `origin` (remoto foi limpo com filter-repo). Publicar.
   - A chave `~/.ssh/id_ed25519` ("git wsl", já cadastrada no GitHub) tem **passphrase** e o ambiente do
     agente não tem `ssh-agent`. Push tem que sair de um terminal do usuário:
     `eval "$(ssh-agent -s)" && ssh-add ~/.ssh/id_ed25519 && git push origin master`.

### Tier 2 — validação biológica + refinamento

3. **Validar as chamadas contra verdade conhecida.**
   - Porquê: confirmar que nossos loci de piRNA/siRNA batem com a biologia.
   - Como: cruzar os loci `piRNA` (GFF) com clusters de piRNA conhecidos de Drosophila (piRNAclusterDB; flamenco em `NT_037436`/X), e os `siRNA` 21nt bifita com loci de endo-siRNA/TE. Usar `bedtools intersect`.
   - Pronto quando: sobreposição alta dos clusters conhecidos e sanidade das assinaturas (1U em piRNA, ping-pong em loci germinativos).

4. ~~**Endurecer a guarda de `pingpong_frac` em loci minúsculos.**~~ **FEITO (2026-08-01)** — mas o problema
   real era outro.
   - A guarda em si virou não-problema: com `--mincov 10` (item 1), `pingpong_frac=1.0` espúrio caiu para
     4 loci na corrida inteira. `features.build_record` ficou como está (`total_pairs≥10`).
   - O que a investigação achou: **1.756 loci eram piRNA APENAS pelo termo de ping-pong**, e nenhum deles
     tinha tamanho de piRNA (813 de 21nt, 275 de 19nt, 251 de 22nt; 43% bifita). O 1U deles (0,61–0,66)
     é o dos siRNA (0,606), não o dos piRNA sólidos (0,741) — ou seja, endo-siRNA clássico de Drosophila
     sendo chamado de piRNA. Causa: reads de 21nt em fitas opostas pegam sobreposição 5'-5' de 10nt por
     geometria, e `w_pingpong: 1.5` (o maior peso) + `w_pingpong_penalty` no siRNA viravam a decisão
     contra uma chamada de tamanho inequívoca.
   - Correção: os termos de ping-pong agora são **multiplicados pela janela de tamanho da própria classe**
     (`decider._score_class`). Impacto medido: 644 loci (1,09%) mudam de classe, piRNA de 18–23nt cai de
     9.100 → 8.456, e os **37.978 piRNA de 24–31nt ficam intactos**.
   - Sobra: ~7.300 loci de 18–23nt ainda são piRNA por `w_1u` (que também não é condicionado ao tamanho).
     Avaliar o mesmo gate no 1U — decidir junto com a validação (#3).
   - Ferramenta nova: simulador do `RuleScorer` sobre o `.clusters.tsv` (reproduz classes e `review` com
     fidelidade 100%) permite avaliar variantes de regra em segundos, sem re-rodar o BAM.

5. **RPM consistente com a ponderação.**
   - Porquê: `lib_size` = `total_mapped` (linhas de alinhamento ≈ 4× reads únicos com `-k5`); RPM fica sub-escalado (consistente entre loci, mas escala absoluta off).
   - Onde: `classify.run_classify` — `lib_size` = nº de reads ÚNICOS (Σ 1/NH). Baixa prioridade (RPM é só reportado, não usado na decisão).

6. **Refinar/validar o phasing hipergeométrico.** Parametrização atual é aproximação v1 (janela 10 ciclos). Validar com locus fásico conhecido (tasiRNA) ou marcar coluna como experimental.

### Tier 3 — features maiores (specs próprios)

7. **Módulo pós-miRDeep (rota miRNA).**
   - Escopo: curar `result_*.csv` do miRDeep2 (score, star>0, randfold p, descartar rfam alert; MFEI/AMFE via RNAfold), isomiRs, e **alimentar `--mirna-loci`** do classify (reusa `mirna_filter.py` por overlap). miRDeep2 continua o motor; nós curamos/integramos.
   - Ver decisões em `srnaprofiler-tools-research` (memória).

8. **Plots por classe no reporter.** Implementado via subcomando `classplots` (distribuição de tamanho, ping-pong, phasing, logos de sequência e cobertura por classe, reutilizando `plots.py`/R); também inclui resumo QC em matplotlib.

9. **`MLDecider` (fase 2 de ML).** RandomForest/XGBoost sobre o `.clusters.tsv` (que já é a matriz de features); interpretar com SHAP. Precisa de loci rotulados (da validação #3). Mesma interface do `RuleScorer` (plugável) — troca sem tocar em medição.
   - Features extras disponíveis pra adicionar aqui (estavam no catálogo da pesquisa, ainda não extraídas): **k-mers**, **entropia de Shannon** da sequência, e **MFE/AMFE/MFEI** (estrutura via ViennaRNA — sobretudo p/ candidatos a miRNA). Baratas de somar ao `FeatureExtractor` quando o ML pedir.

10. **Clustering por DBSCAN (fonte plugável alternativa).** Se as ilhas de cobertura fragmentarem clusters grandes/esparsos de piRNA. `clustering.py` já é a fonte; adicionar DBSCAN por parâmetro.

11. **Alternativa de multimapeador estilo piRAT (read representante por locus).** Além do 1/NH (feito). Relevante porque o NH via XM é capado no `-k` (peso aproximado quando a multiplicidade real > k). Baixa prioridade.

---

## Ponteiros internos
- **Como funciona + fórmulas:** memória `srnaprofiler-pipeline-reference`.
- **Ferramentas de referência estudadas** (em `ref_tools/`, gitignored) e decisões: memória `srnaprofiler-tools-research`.
- **Auditoria da heurística Perl antiga:** memória `srnaprofiler-fidelity-audit`.
- **Dados de teste:** DRR351382 (ENA), BAM local ordenado em `~/dros.sorted.bam`; resultados em `/mnt/d/srnaprofiler-teste/results` e `~/classify_out`.
- **Convenções:** código/comentários em inglês; autor sempre "Lucas" (nunca IA, nem em commits).
