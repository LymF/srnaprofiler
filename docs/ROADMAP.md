# srnaprofiler — Roadmap / Próximos passos

**Autor:** Lucas
**Atualizado:** 2026-08-05

Documento vivo para alinhar sessões. Estado atual + backlog priorizado. Detalhes de COMO a
pipeline funciona (todas as fórmulas) estão na memória `srnaprofiler-pipeline-reference`.
Design/planos: `docs/superpowers/specs/2026-07-22-classification-engine-design.md` e
`docs/superpowers/plans/2026-07-22-classification-engine.md`.

---

## Clustering reformulado — FEITO (2026-08-04)

Spec `2026-08-04-clustering-design.md`, plano `2026-08-04-clustering.md`, 8 tarefas, na
branch `evidence-scores`. 200 testes passando.

- Ilha cresce `pad` **dos dois lados** e só funde com vizinho de perfil compatível (mesmo
  tamanho dominante, `frac_top` dentro de 0,5). Vão tolerado passa a ser 2×pad.
- `--mincov` agora em **RPM** (default 1,0), convertido uma vez por corrida pelo
  `lib_size`. A calibração de 2026-08-01 (contagem bruta, `--mincov 10`) fica **superada**.
- `clustering.pad` vive no config; **500 nos dois reinos** (era 75).
- miRNA só toma o rótulo do locus se cobrir ≥50% dele.
- Colunas novas: `size_purity` e `dicer_call_2`.

**Antes → depois** (nos defaults, `mincov 1.0` RPM e `pad 500`):

| | Drosophila antes | depois | Arabidopsis antes | depois |
|---|---|---|---|---|
| span mediano | 88 bp | **546 bp** | ~100 bp | **605 bp** |
| reads ponderadas/locus | 34 | **63** | 34 | **124** |
| cobertura por base no locus | ~2× | **4,2×** | — | **6,4×** |
| sequências distintas (`n` dos testes) | ~40 | **163** | — | **115** |
| review | 79% | **14,2%** | 33,1% | **3,7%** |
| loci conhecidos | 85/90 piC | 84/90 | 6/6 TAS | 6/6 |
| massa piRNA nos piCs | 3,0× | **3,3×** | — | — |

**Unidade:** contagens NÃO são comparáveis com ferramentas que alocam a multi-read a uma
única posição. Aqui a read vai 1/NH para cada lugar; o ShortStack põe tudo num só, então
as contagens dele são maiores pelo fator de multiplicidade (3,19× nesta biblioteca). Os 63
reads ponderados equivalem a ~200 naquela unidade. A comparação "34 contra 444" que estava
neste documento misturava as duas escalas.

**Correção do que foi reportado em 2026-08-04:** os ganhos anunciados (span 1.025 bp,
122 reads/locus) eram **artefato de bug** — vizinhos incompatíveis se sobrepunham em até
2×pad (95% dos loci) e `_assign_reads` contava a mesma read em vários loci. Corrigido em
`bf63d22`. Igualmente, o `|frac_top−0,5|` do siRNA **não** caiu para 0,113; medido sem o
bug é ~0,20 contra 0,262 do baseline. A tese de que o clustering resolveria o termo de
fita não se sustentou — a branch foi mesclada pelos ganhos de review e enriquecimento,
que são reais, não pelo motivo alegado.

**A cadeia de causa da fragmentação** (três correções, nesta ordem):
1. `bf63d22` — loci sobrepostos: vizinhos que não fundem dividem o vão no meio.
2. `cba0def` — gate de tamanho por igualdade exata da moda barrava 87,7% dos pares
   alcançáveis no X, com diferença mais comum de **1 nt**; trocado por janela de classe.
3. `d82f3f9` + o estimador por massa — a **moda é o estimador errado** com ~5 reads
   (p25 = 1,4) por ilha crua, e 16,3% das ilhas tinham moda fora de toda janela. Passou a
   ser a janela de classe que concentra >50% da massa de reads, estável em 76% das ilhas
   com <10 reads. `clustering.profile_min_reads` = 2.

**Lição metodológica — o recall de piC premia fragmentação.** Ele pergunta se uma região
de 10–40 kb contém ≥1 locus piRNA não-review, então quanto mais picados os loci, maior a
chance de sobrar um fragmento com tamanho de piRNA. A corrida que deu 87/90 tinha loci de
111 bp. Com loci maiores, a população dominante ganha e sub-regiões genuinamente de 22 nt
dentro dessas referências são corretamente chamadas siRNA — três delas têm **88–96% da
massa em 20–22 nt e 0–2% em 24–31 nt**, com 14k–40k reads. Não era regressão, era acerto.
**Para calibrar clustering, usar o enriquecimento** (soma massa) e não o recall (conta
ocorrências). Um limiar de perfil alto (10) chegou a borrar as classes: piRNA 2,8× contra
siRNA 2,7×, separação ~1,0.

### O termo de fita — RESOLVIDO (2026-08-06)

`s_fractop_bal` era um TOST de equivalência afirmando `frac_top` em 0,5 ± 0,15. Tinha
**dois** defeitos, ambos medidos nos loci limpos:

- **Construto errado:** bifita não é 50/50 — um cluster pode ser genuinamente bidirecional
  a 70/30 (viés de transcrição e de carregamento em AGO). 68,5% dos loci de siRNA têm ≥20%
  na fita minoritária, mas só 36,8% caíam dentro da margem.
- **Sem poder:** um locus perfeitamente equilibrado com 50 sequências distintas tirava
  **1,68**, abaixo do 2 que ele deveria passar. Disparava em 4,4% dos siRNA contra 5,5% dos
  piRNA — não discriminava nada e era ~0 para todos, então o `score` da classe siRNA era
  zero **por construção**. Era este o sintoma que travou a branch `evidence-scores`; a
  causa era o termo, não os loci.

Substituído por `s_fractop_bid`: binomial unilateral da fita minoritária contra nulo de
0,10. Dispara em 68,2% dos siRNA e 57,5% dos piRNA. **Score mediano do siRNA: 0 → 4,79**,
com 61,5% dos loci passando `score >= 2`. As duas colunas seguem no TSV; só a
`defining_evidence` mudou, nos dois configs.

Fica em aberto (do mesmo diagnóstico): o piRNA não tem `w_fractop` porque `frac_top` é
bimodal na classe (flamenco 0,90 uni-fita somático; 42AB 0,50 bifita germinativo). A saída
seria `defining_evidence` com duas rotas — bifita COM ping-pong, ou uni-fita SEM — que a
sintaxe de slots OR já suporta.

### O score do phasiRNA — ENCERRADO, não é defeito (2026-08-06)

Com o clustering corrigido, **4 dos 6 TAS ficam significativos**: TAS1A 5,7e-4, TAS1C
7,7e-4, TAS2 2,0e-3, TAS1B 8,0e-3 (com os bugs davam 0,08–0,13).

Os dois fracos — TAS3 (0,034) e TAS3B (0,455) — foram diagnosticados lendo o registro de
fase direto do BAM, e **a hipótese da janela de 11 ciclos foi refutada**:

| locus | massa 21 nt | massa no registro dominante | ciclos ocupados | composição |
|---|---|---|---|---|
| TAS1C | 40.681 | **84,7%** | 12 (em 22) | 69% de 21 nt |
| TAS2 | 6.620 | 34,8% | 9 (em 13) | — |
| TAS3B | 1.106 | 35,0% | 4 (em 5) | 43% 21 nt + 42% 22 nt |
| TAS3 | 793 | **15,1%** | 4 (em 12) | 32% 21 nt, resto 22/24/27 |

No TAS3 os três registros mais fortes têm massa quase igual (120/116/109) — **não há
registro dominante**, e a massa total do locus é 24× menor que a do TAS1C. Não é a janela
diluindo: nesta biblioteca o TAS3 tem pouca expressão e seus 21-mers não estão em fase.
O `phase_test` está certo e o `score` baixo é a leitura correta da evidência — que é
exatamente o trabalho da coluna: separar locus chamado por peso de locus com evidência.

*Em aberto (menor):* se essa baixa expressão do TAS3 é esperada em botão floral (o tecido
de SRR30831721) ou indica problema de amostra — precisa de referência de expressão por
tecido, não verificada.

*Lembrete conceitual:* TAS **não** são miRNA. O miRNA é o gatilho, codificado em outro
locus (MIR173 para TAS1/TAS2, MIR390 para TAS3); o locus TAS produz tasiRNA pela via
RDR6/DCL4. A classe `phasiRNA` está certa.

**Ampliar o conjunto de referência** (sugestão do usuário): usar o **sRNAanno**
(Chen et al., 2021, https://doi.org/10.1038/s41438-021-00480-8) para ter dezenas de loci
PHAS anotados em várias espécies em vez de seis TAS — seis loci, dos quais dois pouco
expressos, é base estreita demais para calibrar peso.

### Multimapeamento — investigado e encerrado (2026-08-06)

O `-k 5` do BAM antigo foi trocado por `-k 20` (default da tool). Multiplicidade média
4,36 → **14,36**, e ainda assim **48,6% dos alinhamentos ficam no teto de 20** — a
multiplicidade real de metade das reads é maior que isso, o que é a cara de cluster de
piRNA (transpóson com centenas de cópias).

**Três coisas que se pensava e que a medição desmentiu:**

1. **Não há viés de contagem no 1/NH.** Com `-k` 5 ou 20 a biblioteca ponderada é
   *idêntica* (9.658.002) e `reads` por locus também (63). O teto reduz a chance de o locus
   ser sorteado mas aumenta a fatia na mesma proporção — a esperança se conserva. O teto
   adiciona **ruído**, não viés. (Chegou a ser afirmado aqui o contrário, em ambas as
   direções; ambas erradas.)
2. **O ganho do `-k` maior foi no `n` dos testes:** `distinct_seqs` 163 → **365**. Os testes
   binomiais (fita, base 5') usam `distinct_seqs`, não `reads`, então o poder dobrou. Era
   justamente o que se supunha que NÃO mudaria.
3. **Alocação única guiada por reads únicas não se sustenta neste organismo.** Ela precisa
   de sinal único como guia, e **57,3% dos loci de piRNA têm ZERO reads únicas** (fração
   única mediana 0,0000; só 17,1% chegam a 5%). No conjunto que importa o guia é nulo e o
   algoritmo cai no fallback 1/NH. Pior: nos 17% com sinal, o peso seria puxado para eles e
   drenado dos 57% sem — e a ausência de read única ali é a **norma biológica** da região
   densa em transpóson, não sinal de locus falso. Seria introduzir viés, não corrigir.
   Custo adicional: duas passadas em 72M alinhamentos e um dicionário de milhões de
   sequências por worker, numa máquina que já teve processo morto por OOM (item 16).

O que sobrou de valor é uma coluna, não uma arquitetura: **`unique_frac`** (`59346cd`),
a fração da massa do locus vinda de reads únicas, na escala 1/NH — ao contrário de
`multimap_frac`, que conta linhas de alinhamento e por isso pesa mais justamente as reads
menos confiáveis. Marca os loci cuja assinatura é toda emprestada de repetição.

*Contexto de método:* 1/NH não é invenção nossa — é o `normalize by genomic hits` do
tiny-count e o modo `f` do ShortStack. O modo `u` do ShortStack (colocação guiada,
Johnson et al. 2016, https://doi.org/10.1534/g3.116.030452) é o que foi avaliado e
descartado acima.

**Pendência achada na corrida:** `evidence.py:18` emite
`RuntimeWarning: overflow encountered in scalar divide` quando p ≈ 0 (o `1/p` estoura antes
do `min` com o teto). Não altera resultado; limpar.

---

## RETOMAR AQUI (sessão de 2026-08-05)

Tudo em **`master`** (a branch `evidence-scores` foi mesclada em `811390b` e apagada).
214 testes passando.

### O que está pronto

1. **Score de evidência**: `evidence.py` + colunas `score` e `s_*` no TSV. Reportado,
   não decide. Ressalva do `s_fractop_bal` na seção acima.
2. **Clustering**: reformulado + 3 bugs corrigidos depois (offset do dúplex, loci
   sobrepostos, gate de tamanho). Ver a seção acima, com os números honestos.

### Ordem da próxima sessão

1. ~~Remapear com `-k` maior~~ **FEITO (2026-08-06)** e ~~alocação única de multi-read~~
   **DESCARTADA com evidência** — ver a seção "Multimapeamento" abaixo. O BAM de referência
   agora é `~/dros_k20.sorted.bam` (`-k 20`).

2. **Puxar os PHAS do sRNAanno** para sair das seis referências de TAS (das quais duas são
   pouco expressas neste tecido). Só com dezenas de loci faz sentido mexer em `w_phase`.
3. Rota bimodal do piRNA no `defining_evidence` (seção acima).
3. Testar o `pad` em mais organismos: `pirnaclusterDB` mudou para **pirnaclusterdb.org**
   (o domínio uni-mainz saiu do ar), 51 espécies, clusters em GTF preditos com proTRAC
   2.4.4. Ressalva: são preditos, não validados — servem para testar generalização entre
   clados, não para provar acerto.

### Pendências menores

- `git push` (chave `~/.ssh/id_ed25519` tem passphrase; precisa de `ssh-add` num terminal seu).
- Minors diferidos no ledger: guarda de `evidence.density_score` (ordem dos checks com
  `lib_size=None`), e o teste de threads não empurra `genome_bp` não-zero pelo pool.
- **18% dos loci têm duas populações de tamanho** (pureza < 0,5), já em `pad=75`. A tarefa 7
  do plano passa a reportar (`size_purity`, `dicer_call_2`); decidir depois se vale
  classificar as duas.

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
   - **Impacto medido** (re-medido em 2026-08-01 após o fix `7057776`, que também é do phasing e
     contaminava estes números): **591 loci** mudam de classe (588 piRNA→siRNA, majoritariamente
     18–22nt — a direção certa). Phasing significativo (p<0,001) em **352 loci** (0,6%), concentrados
     em 24–26nt (piRNA por Zucchini) e 21nt. Recall dos piCs e landmarks flamenco/42AB inalterados.
     *(Os valores antes publicados aqui — 546 loci e 427 fásicos — foram medidos com o bug da janela
     vazia ainda ativo. O fix isolado move 45 loci de classe no Drosophila, contra a recalibração
     inteira em planta: lá `w_phase` é 0,5 e gated, aqui é 3,0 e define a classe.)*
   - **Limitação conhecida:** `phase_score` = `ratio × number × log(abun)` (fórmula PhasiHunter) ainda
     correlaciona 0,53 com `log(reads)` — mede abundância junto com fase, então não é comparável entre
     loci. Use `phase_pvalue` para comparar. Normalizar o score fica para o `MLDecider`.

### Tier 2b — suporte a plantas (o segundo reino da tool)

O motor já é agnóstico: `decider`/`features`/`reporter` não têm nome de classe embutido (só
`classify.py:27,70` trata "miRNA", que vem de `--mirna-loci`). Toda a biologia está no `rules.yaml`
— que hoje é **de artrópode**. O que falta:

12. ~~**`rules.plant.yaml` com a taxonomia certa.**~~ **FEITO (2026-08-01), NÃO VALIDADO.**
    Classes: `phasiRNA` (21–22nt, `w_phase: 1.5` como evidência principal, 5'U), `hcsiRNA` (23–24nt, 5'A,
    bifita) e `siRNA` (21–22nt bifita sem registro). Sem classe piRNA e sem termos de ping-pong — planta
    não tem PIWI. `--config rules.plant.yaml` no `classify`. 5 testes fixam a discriminação pretendida
    (inclusive que ping-pong não move nenhum score em planta).
13. ~~**Feature `pct_5a` (primeira base A).**~~ **FEITO (2026-08-01).** Nova coluna no `.clusters.tsv` e
    peso `w_5a` no decider. Ponderada por 1/NH e lida na orientação da read, como `pct_1u`.
14. ~~**Validar em dados de planta.**~~ **FEITO (2026-08-01).**
    - **Dados:** SRR30831721 (Col-0, botões florais, WT do BioProject PRJNA1011929 — estudo de Pol IV/CLSY,
      ou seja, a via RdDM dos 24nt). 25,7M reads; adaptador `AGATCGGAAGAGCACACGTCTGAACTCCAGT` **lido das
      próprias reads** (o default da tool, `TGGAATTCTCGGGTGCCAAGG`, é outro e teria destruído a
      distribuição de tamanhos). 24,1M reads após trimming, 82,7% mapeadas no TAIR10 (`GCF_000001735.4`).
      Ferramentas no env conda `srnaprofiler` (bowtie/fastp/samtools não estão no PATH padrão).
    - **Resultado:** **6/6 loci TAS** chamados `phasiRNA` de 21nt; genoma 82,2% de 24nt, 94,5% `hcsiRNA`;
      `pct_5a` 0,512 em 23–24nt contra 0,330 em 21–22nt (AGO4 vs AGO1). Taxa de review 33,1%.
    - **Calibração** (os pesos de livro-texto do item 12 estavam errados em dois pontos):
      `w_size` do phasiRNA 1.0→**0.5** e `w_phase` 1.5→**3.0** — com `w_size` alto a janela de tamanho
      carregava a classe sozinha e todo locus de 21–22nt virava phasiRNA (6.557 loci, 2,8% com phasing;
      a classe `siRNA` ficava com 311). `phase_scale` 3.0→**1.5**: a sigmoide larga devolvia 0,159 para
      `phase_score` ZERO, dando 0,24 de score de graça. Depois: 1.165 phasiRNA, 2.701 siRNA, 6/6 TAS.
    - **`review_confidence` 0.6→0.43** — e isso é MAIS estrito, não menos. O limiar é sobre softmax, logo
      seu rigor depende do nº de classes: com 2 classes `conf≥0.6` = margem `ln(1,5)=0,405`; com 3 classes
      = margem **1,10**, quase 3× mais severo sem ninguém ter decidido. Era a causa de 39,6% dos loci irem
      a review só por confiança. Review 66,9% → 33,1%.
    - **`mincov`: manter 10 em planta** (o default). Ao contrário do Drosophila, subir PERDE sinal:
      mincov 30 mantém 17.396 loci contra 23.248, com a mesma massa de reads (0,164 vs 0,162). hc-siRNA
      é difuso por natureza (espalhado por TEs), então cortar por profundidade descarta locus real.
    - **Bug encontrado graças a esses dados** (não aparecia no Drosophila): `phase_test` ranqueava as
      janelas por **menor p-valor**, e o hipergeométrico premia janela vazia (`P(X>=1)` com 1 draw = 0,046
      < o p de uma janela com 7 ciclos reais). O TAS1C, locus mais abundante, tinha phasing 91,2 reportado
      como 7,8, e o mesmo locus oscilava 4,15 / 0,69 / 11,84 conforme o `mincov`. Corrigido em `7057776`
      (ranqueia por ciclos em registro, depois score); os valores ficaram estáveis entre mincov.
    - **Bug no `simulate_rules.py`:** `pct_5a` faltava em `FEATURES`, então toda simulação de planta
      subestimava o `hcsiRNA` e inventava 3.080 mudanças de classe. Agora o script **aborta** se falta
      coluna, em vez de assumir 0 em silêncio.

### Rota miRNA validada ponta a ponta (2026-08-02, Arabidopsis)

Fluxo correto: **`qc -g` → `mapping` → `mirna` → `classify --mirna-loci`**. O `qc -g` NÃO é opcional
(padroniza headers e agora mapeia IUPAC→N; o miRDeep2 recusa genoma sem isso).

- **Especificidade do `classify` por cenário** (247 loci que SÃO miRNA anotado):

  | `--mirna-loci` | loci miRNA | anotados corretos |
  |---|---|---|
  | nenhum | 0 | 0/247 (0%) — 184 recebem classe errada AFIRMADA |
  | anotação TAIR10 (428) | 243 | **247/247 (100%)** |
  | de novo miRDeep2 (87) | 85 | 63/247 (25,5%) |

- **miRDeep2 de novo puro** (sem `--mature`/`--hairpin`): 94 predições, precisão 53,2% contra a anotação,
  recall 16,1%. O score separa bem (mediana 96,3 nos que batem, 10,0 nos que não) → **filtrar por score é
  obrigatório**; é a evidência empírica que faltava para o item 7.
- **A recuperação não melhora com abundância** (25% nos loci >10k reads, 34% no top-50), então o gargalo
  não é sensibilidade a baixa expressão: é o modo de novo puro. **Para organismo não-modelo, passar
  miRBase de espécie relacionada** (`--mature`/`--hairpin`, que o miRDeep2 aceita como "other species").
  Falta medir esse ganho.

### Achados operacionais desta rodada

15. **`--create-bam` não dimensiona a memória do `samtools sort`.** `mapping.convert_sam_to_bam` não passa
    `-m`, então usa o default de 768 MB **por thread**: `--threads 12` pede 9,2 GB e morre com
    "couldn't allocate memory for bam_mem" em máquina de 8 GB, em qualquer disco. O que ele produz está
    correto (BAM `SO:coordinate` + `.bai`, verificado) — o defeito é só o dimensionamento.
16. **O `classify` também não tem teto de memória.** Cada thread carrega os alinhamentos de um cromossomo
    inteiro; com 4 threads em Arabidopsis o processo é morto pelo OOM killer numa máquina de 8 GB
    (`--threads 1` usa ~4 GB de pico). Vale um aviso ou um teto automático.
17. **O `/mnt/d` (drvfs) NÃO quebra o `samtools sort`** — testado forçando spill com `-m 40M`: criou 7
    temporários no D:, mesclou e terminou íntegro. A nota antiga do ROADMAP não se sustentou no teste.
    O workaround de manter o BAM em disco Linux talvez seja desnecessário; falta testar o `classify`
    lendo BAM direto do D: (acesso aleatório via `fetch`, que é o padrão de I/O ainda não testado lá).

### Score de evidência estatística — FEITO (2026-08-03)

Spec: `docs/superpowers/specs/2026-08-02-evidence-scores-design.md`.
Plano: `docs/superpowers/plans/2026-08-02-evidence-scores.md`. Branch `evidence-scores`.

- Cada assinatura passa a ser testada contra um nulo explícito e vira
  `score = log10(1/p)` (binomial para base 5' e fita, Poisson para densidade,
  hipergeométrico já existente para fase). Nulo da base 5' = 0,25, como o proTRAC.
- O score da classe é o **mínimo** entre as evidências que a definem
  (`defining_evidence` no `rules.yaml`): score 2 significa que TODAS as
  definidoras têm p ≤ 0,01. Mínimo e não soma, para que um sinal forte não
  encubra um ausente.
- Novas colunas no `.clusters.tsv`: `score` e `s_1u`, `s_5a`, `s_fractop`,
  `s_pingpong`, `s_phase`, `s_density`.
- **Validação:** piC 85/90 e TAS 6/6 mantidos (o score muda a confiança, não
  quais loci existem). Em planta, organela tem score mediano **0,95** contra
  **3,155** nuclear — o falso positivo que a via RdDM proíbe cai sozinho na
  faixa baixa, sem a ferramenta saber o que é organela.
- **`min_evidence` removido** (era de 2026-08-02): fazia o mesmo trabalho com
  limiar sobre soma de pesos e custava o piC `2R:8443186-8476895` (recall
  84/90). Os 8.794 loci que ele descartava têm score mediano 1,17 contra 2,77
  dos demais, então `score >= 2` reproduz o filtro sem perder o cluster.
  Lição: a validação dele usou `extras/simulate_rules.py`, que **não reproduz
  `review`** — o simulador serve para comparar pesos, não mudanças que mexem no
  review; use os validadores reais.

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
