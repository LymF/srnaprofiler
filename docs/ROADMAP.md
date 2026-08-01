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

3. ~~**Validar as chamadas contra verdade conhecida.**~~ **FEITO (2026-08-01).**
   - Referência: `extras/piC_ref_dm6.bed` — 90 clusters de piRNA de ovário conservados (detectados em
     ≥6 das 8 cepas) já remapeados para ISO1/dm6, que é o mesmo Release 6 do nosso genoma, então **não
     precisa de liftover** (Srivastav et al. 2024, ver Referências). Script: `extras/validate_piclusters.py`
     (interseção em numpy — sem dependência de `bedtools`).
   - **Recall: 85/90 (94,4%)** dos piCs de referência contêm ≥1 locus nosso de piRNA não-review.
     Os 5 não recuperados são pequenos ou teloméricos.
   - **Enriquecimento:** os piCs ocupam 5,66% dos braços principais e concentram **16,9% da massa de
     reads de piRNA (3,0×)**, contra 1,4× para siRNA.
   - **Teste discriminativo (o mais forte):** flamenco (uni-fita, somático) vs 42AB (bifita, germinativo)
     — assinaturas opostas conhecidas, e nada no código sabe onde essas regiões ficam:

     | locus | frac_top | \|0,5−frac_top\| | pingpong_frac | 1U |
     |---|---|---|---|---|
     | flamenco (X) | 0,900 | 0,411 | 0,059 | 0,754 |
     | 42AB (2R) | 0,519 | 0,152 | **0,132** | 0,719 |

     O pipeline reproduz a assimetria de fita do flamenco e o ping-pong 2,2× maior do 42AB.
   - **Ressalva sobre o gate do item 4:** dos 478 loci reclassificados piRNA→siRNA nos braços principais,
     33,3% caem dentro de piCs de referência contra 26,7% de fundo (enriquecimento fraco, 1,25×), e
     custaram 1 piC de recall (`2L:2017-22057`, telomérico). Eles têm `frac_top` ~0,55 e ping-pong ~0,07 —
     compatível com endo-siRNA produzido por cluster bifita (Dcr-2 processa dsRNA dessas regiões), mas
     esse é o ponto a reexaminar se o gate for estendido ao `w_1u`.

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

5. ~~**RPM consistente com a ponderação.**~~ **FEITO (2026-08-01).**
   - `samio.weighted_library_size()` (Σ 1/NH num passe pelo BAM) vira o `lib_size` default do
     `classify` quando não se passa `--norm`. Custo: ~70s para 30M alinhamentos.
   - Em DRR351382: **9.658.002 reads ponderadas** contra 30.782.645 linhas de alinhamento — o RPM
     estava sub-escalado por **3,19×** (fator confirmado na rodada real, mediana 3,187).

6. ~~**Refinar/validar o phasing hipergeométrico.**~~ **FEITO (2026-08-01).** Eram 3 bugs + 1 decisão.
   - *(Nota: a versão antiga deste item sugeria validar com tasiRNA — **tasiRNA é de planta**, não existe
     em Drosophila. Em artrópode o phasing relevante é o de piRNA processado pela Zucchini.)*
   - **Bug 1:** `observed` (`phase_number`) era contado no locus inteiro, mas a população era uma janela
     de 10 ciclos. Quando o locus tinha mais ciclos que a janela, a soma do hipergeométrico ficava vazia
     e devolvia **p = 0,0 exato** — matematicamente impossível. Eram 22 loci; agora são **0**.
   - **Bug 2:** `phased_positions` recebia `pl` (18–35). O PhasiHunter passa `phase_length` porque, com
     janela de **11 ciclos**, as posições fásicas nas duas fitas são `2×11−1 = 21` — igual ao `phase_length`
     de 21nt *por coincidência numérica*. Copiamos o número sem a coincidência (nossa janela era de 10
     ciclos, e aplicávamos a qualquer tamanho). Agora é `2×ciclos−1`.
   - **Bug 3:** o registro era montado com TODAS as reads do locus. O PhasiHunter (`phaseFunction.py:460`)
     só usa reads **do tamanho da fase** — uma read de 26nt num locus de fase 21 é ruído no registro.
   - **Novo `signals.phase_test()`:** varre janelas de 11 ciclos em passos de 1 ciclo (slice por `bisect`)
     e fica com a mais significativa, então phasing local não se dilui em ilha grande. População, sucessos,
     draws e observed saem todos da MESMA janela.
   - **Decisão (gate de tamanho no `w_phase`):** antes o peso era idêntico nas duas classes, então somava a
     mesma constante aos dois scores e **se cancelava** — confirmado removendo `w_phase`: 0 mudanças em
     58.909 loci. O phasing era decorativo. Com o gate, ele aponta para a classe cujo tamanho bate, o que
     serve **os dois reinos com um único config**: em planta um locus fásico de 21/24nt é phasiRNA/tasiRNA
     (siRNA), em artrópode um de 24–31nt é piRNA.
   - **Impacto medido:** 546 loci mudam de classe (543 piRNA→siRNA, majoritariamente 18–22nt — a direção
     certa). Phasing significativo (p<0,001) em 427 loci (0,7%), concentrados em 24–27nt (piRNA por
     Zucchini) e 21nt. Recall dos piCs e landmarks flamenco/42AB inalterados.
   - **Limitação conhecida:** `phase_score` = `ratio × number × log(abun)` (fórmula PhasiHunter) ainda
     correlaciona 0,53 com `log(reads)` — mede abundância junto com fase, então não é comparável entre
     loci. Use `phase_pvalue` para comparar. Normalizar o score fica para o `MLDecider`.

### Tier 2b — suporte a plantas (o segundo reino da tool)

O motor já é agnóstico: `decider`/`features`/`reporter` não têm nome de classe embutido (só
`classify.py:27,70` trata "miRNA", que vem de `--mirna-loci`). Toda a biologia está no `rules.yaml`
— que hoje é **de artrópode**. O que falta:

12. **`rules.plant.yaml` com a taxonomia certa.** Planta **não tem piRNA**, e a janela `24–31nt` está
    ocupada por **hc-siRNA** (heterocromático, DCL3/AGO4, via RdDM) — abundantíssimo. Com o config atual
    esses loci sairiam rotulados "piRNA". Também não há ping-pong em planta (sem PIWI), então
    `w_pingpong`/`w_pingpong_penalty` não deveriam pontuar.
13. **Feature `pct_5a` (primeira base A).** Em planta o carregamento em Argonauta é guiado pela base 5':
    5'U → AGO1 (miRNA/tasiRNA 21nt), 5'A → AGO4 (hc-siRNA 24nt). Medimos `pct_1u` (1ª base U) e `pct_10a`
    (**10ª** base A, assinatura de ping-pong), mas não a 1ª base A — a discriminante em planta. Barata de
    somar ao `FeatureExtractor`.
14. **Validar em dados de planta.** Sem um dataset vegetal, os itens 12–13 são especulação informada.

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

---

## Referências

Srivastav, S. P., Feschotte, C., & Clark, A. G. (2024). Rapid evolution of piRNA clusters in the
*Drosophila melanogaster* ovary. *Genome Research*, 34(5), 711–724.
https://doi.org/10.1101/gr.278062.123
