# srnaprofiler — Roadmap / Próximos passos

**Autor:** Lucas
**Atualizado:** 2026-08-08

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

### Fase: viés de seleção de janela — CORRIGIDO (2026-08-07)

`phase_test` varre uma janela de 11 ciclos por ciclo do locus e reporta a melhor, mas
devolvia o **p cru** daquela janela. Reportar o mínimo de W testes como se fosse um teste
único infla significância — isso é propriedade do algoritmo, não achado empírico, e vale
para qualquer organismo.

Medido (taxa de p ≤ 0,01 contra o 1% de um nulo calibrado), **nos dois reinos**:

| ciclos | Arabidopsis | Drosophila |
|---|---|---|
| <15 | 2,1× | 2,8× |
| 15–30 | 8,6× | 11,9× |
| 30–60 | 13,2× | 15,1× |
| >120 | 32,7× | 44,6× |

**Teste de permutação (2026-08-07)** — o nulo, onde a verdade é conhecida por construção:
posições 5' sorteadas uniformemente, nenhum locus fásico, tudo o mais preservado (nº de
reads, tamanho, nº de janelas). Taxa de p ≤ 0,01, 400 réplicas por ponto:

| ciclos | sem correção | com correção | nominal |
|---|---|---|---|
| 12 | 18,8% | 2,75% | 1% |
| 30 | 40,3% | 0,25% | 1% |
| 60 | 65,5% | 1,50% | 1% |
| 120 | **76,3%** | 2,50% | 1% |

Ou seja: **num locus de 120 ciclos sem fase alguma, o teste antigo dizia "fásico, p ≤ 0,01"
em 76% das vezes.** O dado real subestimava o estrago (~10×) porque a maioria dos loci é
pequena. A correção elimina a dependência do tamanho — que era o defeito estrutural.

**Limitação residual conhecida:** sobra anticonservação de ~2–3×, que **não** vem da
varredura e sim do próprio nulo hipergeométrico (population/successes são idealização de
onde uma read pode começar). É limitada e independente do tamanho do locus. A afirmação
anterior de que "a calibração caiu em cima do nominal" era otimista: os 0,89% do dado real
são em parte porque loci reais diferem deste nulo. Fixado em `tests/test_signals_phase_null.py`.

**Ressalvas de escopo** (importantes): as magnitudes são destas duas bibliotecas e não devem
ser citadas como número geral — dependem da distribuição de tamanho dos loci. A relação
dose-resposta é *consistente com* o viés mas não é prova isolada (locus maior poderia ser
genuinamente mais fásico); a prova está no código. E o "1% esperado" pressupõe que a maioria
dos loci não é fásica, hipótese razoável e não verificada.

Corrigido com **Bonferroni sobre as janelas efetivamente testadas** (conservador, já que
janelas vizinhas se sobrepõem — direção certa para evidência). Efeito em Arabidopsis:

| | antes | depois | esperado |
|---|---|---|---|
| p ≤ 0,05 | 31,7% | 6,35% | 5% |
| p ≤ 0,01 | 9,8% | 0,89% | 1% |
| p ≤ 0,01 dentro de PHAS × fora | 19,4% × 9,8% (**2,0×**) | 6,0% × 0,9% (**6,7×**) | — |

Recall de PHAS inalterado (20/35) e TAS 6/6 chamados: tirou falso positivo, não locus real.

**MAS o custo em sensibilidade é grande, e foi medido só depois (2026-08-07).** "TAS 6/6"
significa *chamados phasiRNA*, não *fase significativa*. Contra o conjunto experimental de 8
TAS (`extras/phas_experimental_at_TAIR10.bed` — clonados, verificados por mutante, com
gatilho conhecido), a detecção a p ≤ 0,01 é de **1 em 8** (só TAS3c, 2,3e-4); TAS4 não é
sequer chamado. Números por locus citados antes nesta seção (TAS1A 5,7e-4, TAS1C 7,7e-4)
vieram da corrida `at_wd`, **anterior** ao Bonferroni — registrei valores pré-correção numa
seção que descreve a correção.

**Como isso corrige uma conclusão anterior desta sessão:** o item "score do phasiRNA fraco"
tinha sido encerrado como não-defeito, com base nos 6 TAS. Havia defeito — só que do lado
oposto ao suposto: o teste não era conservador demais nos TAS, era **permissivo demais em
todo o resto**, o que achatava a discriminação. O diagnóstico do TAS3 em si continua válido
(pouco expresso, 15,1% da massa no registro dominante contra 84,7% do TAS1C), mas a
comparação com o fundo estava distorcida. Foram os 35 loci do sRNAanno que tornaram isso
visível — 6 loci não bastavam.

*Lembrete conceitual:* TAS **não** são miRNA. O miRNA é o gatilho, codificado em outro locus
(MIR173 para TAS1/TAS2, MIR390 para TAS3); o locus TAS produz tasiRNA via RDR6/DCL4.

*Referência nova:* `extras/validate_phas.py` + `extras/phas21_ref_at_TAIR10.gff3` (35 loci,
sRNAanno — Chen et al., 2021, https://doi.org/10.1038/s41438-021-00480-8). O servidor serve
`http://www.plantsrnas.org/PHAS21/<Especie>.PHAS21.gff3` e `PHAS24/` para gramíneas
(Oryza_sativa tem 804 loci PHAS24), o que permite testar generalização entre espécies.

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

### O que está aberto (revisado 2026-08-07)

**A0. Rota miRNA — quebrou e foi CORRIGIDA (2026-08-07)**

Validando contra o **miRBase** (curado, e cujas coordenadas batem com as montagens já
usadas: `dme` em Release_6, `ath` em TAIR10 — `extras/validate_mirna.py`), apareceu o que
as validações de piRNA e PHAS não podiam ver, porque cada uma olha uma classe só:

- **Só 1 locus** do genoma recebeu rótulo miRNA, com `--mirna-loci` passado. Antes da
  reforma do clustering eram **243 loci, 247/247 corretos** (anotação TAIR10).
- **21,05% da classe phasiRNA está sobre miRNA anotado** (131 de 627). Em Drosophila a
  contaminação é pequena (0,23% do piRNA, 0,99% do siRNA), mas **134 precursores cobertos
  recebem 100% de rótulo errado** (70 siRNA, 64 piRNA; 78 deles em review).

**Causa — interação de duas mudanças desta sessão, ambas minhas:** `mark_mirna` passou a
exigir que o miRNA cubra **≥50% do locus** (tarefa 6 do plano de clustering), regra
calibrada quando os loci tinham ~100 bp. Com loci de 600–1.000 bp e precursores de ~91 bp,
o precursor ocupa **14% do locus (mediana; 6,9% em Drosophila)** e o limiar nunca dispara:

| `min_fraction` | loci que seriam rotulados |
|---|---|
| 0,50 (atual) | 28 (9,2%) |
| 0,25 | 61 (20,1%) |
| 0,10 | 225 (74,0%) |

**Baixar o limiar é remendo.** O conserto certo é o que o ShortStack faz: um locus de miRNA
conhecido deve ser **recortado** do cluster e virar locus próprio, em vez de ser engolido —
"if one or more MIRNA loci overlap the initial cluster, the initial cluster is removed and
only the refined MIRNA region is reported" (README do ShortStack). Isso resolve os dois
lados: o miRNA deixa de ser diluído e a classe phasiRNA deixa de ser contaminada.

**Corrigido com recorte** (`clustering.carve_out`): o `--mirna-loci` passou a entrar
**antes** do merge, cada intervalo vira locus próprio e os flancos seguem como loci deles.
Baixar o limiar seria um número que quebra na próxima vez que o tamanho dos loci mudar — e
ele mudou três vezes só nesta sessão.

**Resultado com anotação curada (miRBase; `dme` em Release_6 e `ath` em TAIR10 batem com as
montagens em uso):**

| | Arabidopsis antes → depois | Drosophila antes → depois |
|---|---|---|
| precursores cobertos rotulados miRNA | 0% → **100%** (244/244) | 0% → **100%** (134/134) |
| loci de miRNA em review | 42 → **0** | 78 → **0** |
| classe phasiRNA sobre miRNA anotado | 21,05% → **0,00%** | — |
| classes piRNA/siRNA sobre miRNA | — | 0,23%/0,99% → **0,00%** |
| loci conhecidos | TAS 6/6, PHAS 20→19/35 | piC 84 → **85**/90 |

**O insumo importa mais do que antes.** Com `mirdeep_pred.bed` (87 predições de novo) só
19,7% dos precursores foram rotulados e a contaminação do phasiRNA mal se moveu (20,77%);
com os 326 do miRBase, foi a zero. E como o recorte cria um locus **por intervalo
informado**, uma predição falsa agora vira locus falso com rótulo afirmado — antes ela só
falhava em rotular. Dos 91 loci criados a partir do miRDeep2, apenas 57,1% caem sobre miRNA
anotado. **Passar anotação curada quando existir; predição de novo só filtrada por score**
(item 7 do Tier 3 ganhou consequência direta na saída).

**CUIDADO ao ler enriquecimento depois desta mudança.** Os 136 loci de miRNA em Drosophila
carregam **33,7% da massa de reads da biblioteca**. Tirá-los das classes encolhe o
denominador, e o enriquecimento "sobe" sozinho: piRNA 3,3×→5,1× e siRNA 2,4×→4,2×. Os dois
subirem juntos é assinatura de mudança de denominador, não de discriminação melhor — não é
ganho. O que é ganho limpo: rotulagem 100% correta, review a zero, contaminação a zero e
piC 84→85. *(Corolário: toda medida de massa feita antes nesta espécie estava distorcida por
um terço de miRNA contabilizado como siRNA/piRNA.)*

**Lição de método:** validar uma classe por vez esconde regressão em outra. A tarefa 6 foi
validada contra piC e TAS, que não têm miRNA, e passou. Cada referência nova desta sessão
revelou um defeito que as anteriores não podiam ver — 35 PHAS acharam o viés de janela,
o miRBase achou esta regressão, e o nulo por permutação dimensionou o estrago (76% de falso
positivo em locus de 120 ciclos). Nenhuma delas era "mais dados da mesma coisa".

> **Fórmulas:** `docs/FORMULAS.md` tem todas, com nulo e limitação de cada uma.
> **Comparação com o campo:** `docs/superpowers/specs/2026-08-07-signature-literature-review.md`.

**A-score. Score por tamanho de efeito — FEITO (2026-08-08)**

O `score` era `−log10(p)` e **não ordenava confiabilidade**: p mistura efeito com tamanho
amostral, então correlacionava 0,629 com `distinct_seqs` e só 0,323 com o efeito. Agora é o
mínimo das **forças** das evidências definidoras, cada uma o limite inferior de Wilson
normalizado pelo nulo, em `[0,1]`. As `s_*` seguem como p-valores.

| | antes | depois |
|---|---|---|
| corr. com o efeito (`pct_1u`) | 0,323 | **0,883** |
| corr. com `distinct_seqs` | 0,629 | **0,137** |
| 1U por faixa de score (Dros.) | não-monotônico | 0,38 → 0,56 → 0,70 → **0,82** |
| 5'A por faixa de score (Arab.) | — | 0,43 → 0,61 → 0,74 → **0,86** |
| piso | — | piC 85/90, TAS 6/6, miRBase 134/134 e 244/244, **0,00% de mudança de classe** |

**O critério de aceite original estava errado, e isso quase passou batido.** Ele media a
fração de loci de score alto dentro dos 90 piC de referência, e o desenho **reprovou**
(23,5/32,0/30,5/27,0% contra 28,3% de fundo). Investigando em vez de aceitar: os **4.542**
loci de score alto **fora** dos piC são indistinguíveis dos 1.689 de dentro (1U 0,82 nos
dois; ping-pong 0,068 × 0,079), e o controle de score baixo tem 1U 0,38 e ping-pong **zero**.
Os piC cobrem 5,66% do genoma e só os clusters conservados entre cepas — exigir concentração
ali **penaliza a tool por achar piRNA legítimo em outro lugar**, que é o trabalho dela.
Critério substituído por "a evidência definidora sobe com a faixa de score", que passa nos
dois reinos.

**Escalas diferem por reino, e é honesto:** mediana 0,47 em Drosophila contra **0,127** em
Arabidopsis. O `hcsiRNA` é definido só por `5a`, e 0,43 contra um nulo genômico de 0,332
deixa pouco espaço; o piRNA tem a rota `1u OU pingpong`. O score está dizendo que a evidência
de hc-siRNA é mais fraca — verdade que a escala de p-valor escondia.

**Também endurecido:** `defining_scores` devolvia 0,0 em silêncio para nome de evidência
desconhecido; agora `RuleScorer` recusa o config na construção. *Pendência:* `KNOWN_EVIDENCE`
ainda aceita `fractop` e `fractop_bal`, que não têm `e_*` — um config com esses nomes passa
e pontua 0,0, exatamente o furo que a guarda queria fechar.

**A-5p. Base 5' com nulo genômico — FEITO (2026-08-07)**

O nulo era 0,25 fixo, atribuído aqui ao proTRAC. **A atribuição era falsa** — o proTRAC usa
limiares de fração (`-1Tor10A 0.75`) e não faz teste nenhum de base. O 0,25 era nosso, e
pressupõe uma composição que nenhum genoma tem: o real é **0,272** (Drosophila) e **0,332**
(Arabidopsis).

Agora o nulo é a composição do DNA do próprio locus, ponderada por `frac_top`
(`samio.base_null`) — *"se o corte caísse em qualquer posição daqui, com que frequência o
produto começaria com esta base?"*. É o que faz o `-g/--genome` deixar de ser ignorado: ele
era declarado "optional in v1" e nunca lido.

| | Drosophila | Arabidopsis |
|---|---|---|
| loci conhecidos | piC 85/90, miRBase 134/134 | TAS 6/6, miRBase 244/244 |
| **classe mudou** | **0,00%** | **0,00%** |
| evidência | `s_1u` 61,4 → 52,0 | `s_5a` 7,5 → 4,0 |
| `score` mediano | 57,5 → 50,8 | 7,38 → 3,93 |

Nenhum locus troca de classe porque o classificador usa `pct_1u` cru (`decider.py:63`), nunca
`s_1u` — a mudança é de **honestidade da medida**, não de acurácia da chamada. Em Drosophila a
queda tem dose-resposta com a composição do DNA (0,0 em loci com T/A < 0,28; 24,1 em
0,32–0,36), ou seja, remove falso positivo específico em vez de apertar tudo.

*Atalho rejeitado por medição:* usar as posições 3–5 das próprias reads como controle interno
concorda com o DNA em Drosophila (0,269 × 0,272) mas **falha em Arabidopsis** (0,255 × 0,332,
r = 0,13). Sem o genoma para conferir, eu teria implementado o atalho.

*Erro de processo registrado:* a primeira tentativa passou o genoma por todo o caminho **menos**
a chamada em `run_classify`; o `replace` mirou texto que outro commit já tinha mudado e falhou
em silêncio. As duas validações voltaram byte a byte idênticas e só por isso o bug apareceu —
se eu tivesse olhado apenas o piso (piC/TAS/miRBase, que não mudam mesmo), teria dado por boa
uma feature desligada. Agora há teste de fiação ponta a ponta.

**A-pp. Ping-pong com nulo empírico — FEITO (2026-08-07)**

`s_pingpong` deixou de usar `binomtest(k, n, 1/20)` (nulo uniforme, invenção nossa) e passa a
estimar a taxa dos **próprios offsets do locus**, 1–20 exceto o 10 — o fundo que já
calculávamos para o z-score e descartávamos. É o que o PingPongPro faz
(`pingpongpro.cpp:695`), e pelo mesmo motivo: offsets não são uniformes.

**Efeito medido (Drosophila):** correlação 0,995 com o anterior; 88,9% dos loci mudam menos
de 1. Landmarks preservados (flamenco 7,6→8,2; 42AB 77,7→79,3), piC 85/90 e enriquecimento
inalterados. **Não é ganho agregado** — quando o espectro é plano, Poisson(média) e
binomial(1/20) coincidem.

Onde importa: nos **682 loci que divergem >10, o empírico é maior em 100%**, e são os de
`pingpong_frac` mediano **0,560** — os de ping-pong mais forte. O binomial diluía a
evidência ali porque o `n` incluía o próprio pico. Ou seja, a mudança é neutra na massa e
melhora na cauda que interessa.

**Por que NÃO seguimos o campo até o fim (medido).** O PingPongPro usa z-score (média **e
desvio**) e estratifica por profundidade. Testando essa forma exata sobre o nosso
`pingpong_z` já calculado, ela **discrimina pior** nos landmarks: flamenco 4,53 × 42AB 13,80
(**3,0×**) contra 8,18 × 79,29 do Poisson (**9,7×**).

O motivo é de contexto, não de qualidade: o z e a estratificação existem porque **ele compara
sítios entre loci no genoma inteiro**, onde a profundidade varia e precisa ser controlada.
**Nós comparamos dentro do próprio locus**, onde ela já é constante por construção. Importar
o método deles traz robustez desnecessária e custa poder.

**Ressalva conhecida:** o Poisson pressupõe variância = média. Medido no fundo de offsets:
mediana **1,3** (aceitável), mas **22,5% dos loci acima de 5** e 2,9% acima de 50 — nesses
ele é anticonservador. A síntese correta seria **binomial negativa**, que absorve
superdispersão sem perder poder. Refinamento, não conserto.

**A0-srna. Ideias do sRNAminer que valem implementar (lidas do jar descompilado)**

Ver `docs/superpowers/specs/2026-08-07-signature-literature-review.md` seção 6 para a
comparação completa. O que vale considerar, em ordem de valor:

1. **Massa total mínima por locus na clusterização** (`sRNAlociMinerCLI`: `totalAbundance >=
   minReadDepth`, default 1000). Nós só temos profundidade por posição e `min_reads` como
   *flag de review* — o locus raso existe do mesmo jeito e depois falta `n` para os testes.
   Um mínimo de massa **na formação do locus** garante poder por construção.
   *Adaptação:* em **RPM** (o 1000 absoluto não transfere entre bibliotecas) e **por reino**
   no `rules.yaml`, como o `clustering.pad` — planta tem hc-siRNA difuso onde subir corte
   PERDE sinal (medido em 2026-08-01), Drosophila tem cluster denso que tolera corte alto.
2. ~~Tolerância de ±1 nt no registro de fase~~ **TESTADA E REPROVADA (2026-08-08).**
   Aplicada ao p-valor, a detecção nos 8 TAS experimentais cai de **1/8 para 0/8**: somar os
   registros vizinhos triplica os sucessos (`3M`), então o nulo fica três vezes mais fácil de
   acertar por acaso e o ganho em `k` não compensa. No sRNAminer ela entra só na **razão de
   abundância reportada** (`inPhasRatio`), não no p-valor — eu a apliquei onde ela não
   pertence. Se for reaproveitada, é como métrica descritiva, não como teste.
3. **Deduplicação por (fita, 5') antes de contar** — eles fazem via `seenStringHits`; nós
   contamos posições distintas, que é equivalente. **Já alinhados**, sem ação.
4. **Encadeamento por vão em vez de ilha de cobertura** — filosofia diferente da nossa e da
   do ShortStack. **Não** adotar: nosso merge condicionado por perfil já resolve o mesmo
   problema com mais informação, e foi medido (fundiu 26% sem perder cluster conhecido).

*Nota de método:* a chamada de classe por **fração de tamanho** é padrão universal (proTRAC
`-pisize`, ShortStack `DicerCall`, sRNAminer `minRepeatRatio`, nosso `frac_in_window`).
Nisso não há divergência a resolver.

**A. Pronto para fazer**

1. ~~Rota bimodal do piRNA no `defining_evidence`~~ **FECHADO sem mudança (2026-08-07).**
   A premissa estava errada: `[1u OR pingpong]` **já é** a regra bimodal — 1U é a assinatura
   do piRNA primário (somático, uni-fita) e ping-pong a do secundário (germinativo, bifita).
   O que se propôs acrescentar foi exigir coerência de fita em cada rota
   (`min(s_fractop, s_1u)` OU `min(s_fractop_bid, s_pingpong)`), e isso **mede pior**:
   `score >= 2` cai de 97,0% para 87,1% e 2.040 loci perdem evidência. O perfil deles é de
   piRNA legítimo — 25 nt medianos, 1U 0,67, 43 reads — e **36,2% caem dentro de piC
   conhecido contra 29,4% de fundo**, ou seja, são enriquecidos em cluster real. São
   clusters bifita pouco abundantes falhando a rota germinativa por não exibir ping-pong,
   que com 43 reads não tem pares para ser detectado. Mesma armadilha do `s_fractop_bal`:
   teste sem poder virando evidência de ausência.
   *O que continua válido:* não dar ao piRNA um **peso** único de fita (`w_fractop`), porque
   aí sim a bimodalidade quebraria um dos dois modos (flamenco 0,90 × 42AB 0,50).
2. **Decidir as duas populações de tamanho.** 18% dos loci têm `size_purity` < 0,5. Hoje
   são reportadas (`size_purity`, `dicer_call_2`) mas o locus recebe um rótulo só.
3. **Generalização entre espécies** — agora com base sólida (rota miRNA consertada, fase
   calibrada). Cada espécie nova deve trazer classe ou mecanismo que as atuais não
   exercitam, não "mais do mesmo":
   - **Arroz** (`Oryza_sativa`): PHAS**24** (804 loci no sRNAanno) — fase de 24 nt meiótica,
     mecanismo distinto do 21 nt de Arabidopsis. O sRNAanno serve o próprio genoma
     (`http://121.37.229.61:84/sRNAannoDATA/genome/Oryza_sativa.genome.fasta`, 363 MB), então
     as coordenadas batem por construção. miRBase `osa` tem 595 precursores, **mas em MSU7** —
     conferir se casa com o genoma do sRNAanno antes de usar.
   - **Tribolium castaneum**: 148 clusters de piRNA no piRNAclusterDB
     (`https://www.pirnaclusterdb.org/data/FASTA/Tribolium_castaneum.piRNAclusters.gtf`),
     genoma ~165 MB (cabe na RAM), nomes `ChLG*`. Ressalva: são predições do proTRAC.
     miRBase `tca` deu HTTP 500 — verificar.
   - Falta em ambos: biblioteca de small RNA (SRA/ENA; não há sra-tools no ambiente, baixar
     fastq direto do ENA por HTTP) e índice bowtie novo. É meio dia de trabalho.

3b. ~~Controles negativos~~ **MEDIDOS (2026-08-07)** — `extras/validate_structural.py`,
   contra tRNA/rRNA/snoRNA/snRNA do GFF de cada organismo. **Decidir o que fazer** com o
   resultado (ver abaixo).

**A0b. O p-valor de fase é pouco sensível — a formulação, não a correção (2026-08-07)**

Medido sobre os 73 loci que caem em PHAS anotado (verdade conhecida): só **5,5% alcançam
p ≤ 0,01** e o p mediano é **0,368**. Com ≥200 reads sobe só para 14,8%. Em loci que a
referência diz serem fásicos, detectamos 1 em 20.

**Não é o Bonferroni.** Duas evidências: (i) ~28 janelas num locus de 600 bp, então p
corrigido 0,368 implica p **cru** ≈ 0,013 — não passaria em 0,01 nem sem correção; (ii) o
teste de permutação mostra taxa de falso positivo **no** nominal (0,25–2,75%), não abaixo,
o que aconteceria se estivéssemos sobre-corrigindo.

**O campo faz igual a nós — e sem corrigir.** Lido o código: o PhasiHunter
(`phaseFunction.py:472`) é hipergeométrico sobre **posições** (`total_n` = posições
ocupadas, `total_k` = em registro); a abundância só entra depois como filtro
`min_read_num`. O PHASIS guarda a mesma estrutura `n`/`k`. E **nenhum dos dois corrige a
varredura de janelas**, então os p deles são inflados do mesmo modo que os nossos eram.
Conclusão desconfortável mas honesta: **hoje somos mais corretos que o campo e menos
sensíveis.** Afrouxar a correção seria voltar ao erro deles; o ganho legítimo tem de vir de
estatística mais potente.

*(O sRNAminer, que gerou a anotação do sRNAanno, reporta `Phas_Ratio` — abundância — e p de
1e-9 a 1e-15. Não temos o código. Comparar nossa detecção com a anotação dele mede quão
liberal ele é, não quem está certo: por isso o conjunto experimental acima.)*

**A formulação ignora onde está o sinal.** O hipergeométrico conta *quantas posições* estão em registro e ignora
**quanta abundância** cada uma tem — uma posição em registro com 500 reads pesa igual a uma
com 2. O sRNAminer reporta 1e-9 a 1e-15 nos mesmos loci onde reportamos 1e-3 a 5e-2.

**Isto explica a discordância interna da tool:** o classificador usa `w_phase`, sigmoide
sobre `phase_score`, que **usa** abundância — por isso chama 633 loci. O `score` usa o
p-valor, que **não** usa — por isso quase nenhum passa. Os dois medem coisas diferentes.

**Consequência para o resto:** enquanto isto não for consertado, não dá para usar `score`
para negar rótulo (ver A1), porque 4 dos 13 loci em TAS têm score 0 — seriam descartados
loci reais.

**Resolvido (2026-08-08), spec em `docs/superpowers/specs/2026-08-08-phase-evidence-design.md`.**
A premissa acima — que faltava incorporar abundância — foi medida e rejeitada: o
`phase_ratio` cru (razão de abundância no registro) discrimina TAS/PHAS a AUC 0,415, pior
que moeda. O `phase_score` (que já usa abundância via `w_phase`) sai a 0,786 (TAS) e 0,562
(PHAS) contra o fundo não pareado — quase todo esse número já cai contra o casamento por
`distinct_seqs` sozinho: 0,754 (TAS, quase intacto) e 0,526 (PHAS, no acaso). Ou seja, casar
por profundidade **quase não move o TAS**; o que derruba o `phase_score` é o conjunto PHAS,
maior e mais heterogêneo, pareado ou não. (A cobertura em si é bem maior nos TAS — mediana
1.236 reads contra 85 do fundo não pareado, 14,6x, e 173 dos controles casados, 7,1x — mas
essa diferença de profundidade não é a causa da queda do PHAS.) O que faltava não era
abundância, era medir a *extensão* da corrida em registro em vez da surpresa estatística.
`e_phase` passou a ser `min(ciclos_em_registro, 11) / 11` com `p₀ = 1/phase_length` (antes
era `observados/sorteios` com `p₀ = sucessos/população`), mantendo o `phase_test` e o
p-valor intocados. AUC do `score` (`extras/validate_phase_evidence.py`, controles casados por
`distinct_seqs`): 0,618 → 0,861 em TAS (IC 95% bootstrap [0,738, 0,955], só 12 positivos),
0,639 → 0,810 em PHAS (IC 95% [0,736, 0,879]).

**A2. As quatro colunas órfãs — medido o que cada uma agrega (2026-08-08)**

`overlap_dup_frac`, `size_purity`, `dicer_call_2` e `unique_frac` são reportadas e não entram
em regra nenhuma. Medido nas corridas finais:

| | Drosophila | Arabidopsis |
|---|---|---|
| `overlap_dup_frac` | piRNA 0,0175 × siRNA 0,0168 — **não discrimina** | hcsiRNA 0,0286 × phasiRNA 0,0176 |
| `size_purity` | miRNA **0,598** × piRNA 0,190 | hcsiRNA 0,647 × phasiRNA 0,429 |
| `unique_frac` | miRNA **1,000** × piRNA **0,000** | hcsiRNA 0,845 × phasiRNA **0,093** |
| `dicer_call_2` | presente em **100%** | presente em **100%** |

**Destino de cada uma:**

1. ~~**`dicer_call_2` está quebrada como projetada.**~~ **CONSERTADO (2026-08-08, `c02c6e6`).**
   Disparava em **98,7%** dos loci em Drosophila e 94,0% em Arabidopsis, e metade desses
   segundos tamanhos (48,5% e 69,5%) era um comprimento **vizinho** do dominante — a
   dispersão da mesma população, não uma segunda via. A causa não era o limiar 0,8: era medir
   pureza sobre **um comprimento só**, quando o produto de Dicer se espalha ±1 nt, o que
   chamava de impuro um locus 21 nt limpo. Agora `size_purity` é a banda do pico ±1 nt e
   `dicer_call_2` só sai para um pico a ≥2 nt com ≥20% da massa na janela. Colunas
   diagnósticas; nenhuma regra as lê.
2. **`overlap_dup_frac` não separa o que se esperava.** A ideia era distinguir siRNA (dúplex
   de Dicer) de piRNA (ping-pong) em artrópode; dá 0,0175 × 0,0168. Coerente com o já medido:
   não há pico de dúplex em Drosophila. Documentar como medida que não discrimina aqui.
3. **`size_purity` e `unique_frac` discriminam forte, mas NÃO devem virar peso de classe.**
   Elas são **consequência** da classe, não evidência dela: piRNA vive em transpóson, logo
   `unique_frac` ≈ 0. Usar como peso ensinaria "repetitivo ⇒ piRNA", que é propriedade de
   **onde** o piRNA mora, não prova de que o locus **é** piRNA — a mesma armadilha da
   alocação guiada que foi descartada, e que penalizaria locus real por estar exatamente onde
   deveria estar.
4. **Onde elas servem é na REJEIÇÃO, não na decisão.** Um locus com `size_purity` 0,19 tem
   duas populações misturadas: o rótulo único já está errado, seja qual for a classe. Isso é
   evidência direta e não-circular para a categoria "não classificado" (A1) — não diz *qual*
   classe, diz que *nenhuma* descrição única serve. `unique_frac` fica como diagnóstico
   declarado: útil para o leitor saber que a evidência é emprestada de repetição.

**Desenho para A1 — TESTADO E REPROVADO (2026-08-08).** A proposta era rejeitar por `score`
baixo ou `size_purity` baixa. Medido em Arabidopsis, contra os 8 TAS experimentais (que a
rejeição deve **poupar**) e o RNA estrutural (que deve **pegar**):

| critério | perde dos TAS | pega do estrutural |
|---|---|---|
| `size_purity < 0,2` | 6,7% | 33,3% |
| `size_purity < 0,3` | 26,7% | 55,5% |
| `score < 0,05` | **80,0%** | 86,2% |
| `score < 0,2` | **93,3%** | 94,5% |
| `score < 0,2` **E** `pureza < 0,4` | 26,7% | 67,8% |

**Nenhum ponto é aceitável.** O melhor da combinação perde um quarto dos loci de verdade
conhecida. As distribuições se sobrepõem porque em planta o score mediano é 0,127 e os TAS,
sendo `phasiRNA` definido por `phase` — a evidência mais fraca da tool (1/8 a p ≤ 0,01) —
caem na mesma faixa do RNA estrutural.

**Conclusão:** consertar o `score` era necessário mas **não suficiente**. A rejeição por score
só fica viável quando a evidência de fase for forte o bastante para separar phasiRNA real de
ruído — e o levantamento já mostrou que nossa fase está em pé de igualdade com o campo, então
esse ganho não sai de ajuste, teria de vir de dado (biblioteca mais profunda, tecido onde os
PHAS são expressos) ou de uma evidência definidora adicional para a classe.

*Caminho alternativo que não depende disso:* rejeitar por **anotação** em vez de por score —
um `--exclude-loci` com tRNA/rRNA/snoRNA, análogo ao `--mirna-loci`. Pega o alvo real (RNA
estrutural) sem tocar em locus de verdade conhecida. Menos elegante, mas mensurável e sem
falso negativo.

**A2c. Organela — filtro por CONTIG, com flag explícita (desenhado 2026-08-08)**

Medido nas corridas finais: a organela dá **336 loci e 7,9% da massa de reads** em Arabidopsis
(cloroplasto 6,46%, mitocôndria 1,45%), quase toda chamada `hcsiRNA`; em Drosophila são 27
loci e 0,60%. Não é problema só de planta — mitocôndria existe em todo eucarioto — mas em
planta é dez vezes maior, porque o plastídio tem genoma grande e muito transcrito e a planta
não tem via de siRNA de organela: aquilo é degradação sendo rotulada.

*Corrige uma leitura anterior:* o ROADMAP registrava que a organela em planta tem score
mediano 0,95 contra 3,155 do nuclear, como evidência de que o score funcionava. Continua
verdade, mas **score baixo não é ausência** — os 336 loci estão lá, e quem somar massa por
classe em planta está incluindo 8% de organela.

**Natureza diferente do filtro de ncRNA, e por isso é outro mecanismo:**

| | ncRNA estrutural | organela |
|---|---|---|
| origem | Rfam, 829 famílias, 316 MB | o próprio genoma |
| filtra por | **sequência** (`bowtie --un`) | **contig** |
| onde | antes do mapeamento | pode ser no `classify`, o BAM já traz o nome |
| impacto medido | 36,5% / 25,7% da massa | 7,9% (planta) / 0,6% (mosca) |

**Descoberta automática foi considerada e DESCARTADA.** Os cabeçalhos do RefSeq trazem a
palavra (`>NC_000932.1 ... chloroplast, complete genome`), e daria para detectar por
palavra-chave. Mas: (a) o miRDeep2 recusa espaço no identificador, então qualquer preparo de
genoma para ele trunca no primeiro espaço e a palavra some; (b) genoma de montagem própria ou
de outra fonte não garante descrição nenhuma. Funcionaria nos nossos dois genomas e falharia
**em silêncio** no primeiro caso fora disso — o mesmo padrão do `-g` que não fazia nada e do
`dicer_call_2` que dispara sempre.

**Decisão: flag explícita** (`--organelle-contigs NC_000932.1,NC_037304.1`), sem adivinhação.
O log diz quais contigs foram excluídos e quanta massa saiu.

**FEITO E MEDIDO (2026-08-08, `eb793b1`).** No `classify`, antes da extração; contig que não
está no BAM **aborta** em vez de virar no-op silencioso (o erro a evitar é um typo analisando a
organela e reportando corrida limpa). O `lib_size` continua contando as reads excluídas, para
o RPM seguir comparável com corridas que as mantiveram.

**Os dois filtros são complementares, e a medição mostra por quê** (Arabidopsis, denominador
casado):

| | loci de organela | massa deles |
|---|---|---|
| sem filtro nenhum | 336 (1,12%) | **7,91%** |
| só `--ncrna-db` | 290 (1,00%) | **3,08%** |
| `--ncrna-db` + `--organelle-contigs` | **0** | **0%** |

O filtro de ncRNA sozinho tira **61% da massa** de organela — é rRNA de plastídio, que está no
Rfam — mas deixa **290 loci** de pé, porque o resto é degradação de mRNA do plastídio, que
banco de ncRNA nenhum contém. Só o corte por contig zera. Nada mais muda: a corrida perde
exatamente esses 290 loci (29.044 → 28.754) e o piso fica intacto (TAS 6/6, miRBase 242/242).

*Corolário para quem for ler massa por classe em planta:* rodar só o `--ncrna-db` ainda deixa
3,1% da massa em organela, quase toda chamada `hcsiRNA`.

**A2b-pre. O campo filtra RNA estrutural NA READ, não no locus (lido 2026-08-08)**

O sRNAminer **exige** um banco de ncRNA (`--nd`, marcado `[Required]`: rRNA, tRNA, snoRNA,
snRNA) e o usa em `Data_Preprocessing/Quality_control.py:19`:

```
bowtie -v 1 -f -p N <ncRNA_db> <reads> --un <saida>
```

Ou seja: alinha as reads contra o banco e **segue só com as que NÃO alinharam**. Mesma coisa
para organela (`:47`). O filtro acontece **antes** de qualquer clusterização.

**Isto é melhor que o `--exclude-loci` proposto acima**, por três motivos:
1. **Remove a causa, não o efeito** — não existe locus contaminado para depois rejeitar. E
   como medimos que RNA estrutural carrega **36,5% da massa de reads** em Arabidopsis e 25,7%
   em Drosophila, isso conserta junto toda estatística de abundância por classe.
2. **Não precisa de anotação do genoma**, só de um FASTA de sequências de ncRNA — que existe
   para muito mais espécie do que GFF com tRNA/rRNA anotado. Resolve o caso não-modelo que
   motivou este item, sem depender de predizer estrutura.
3. **Não tem falso negativo em locus real**, que foi o custo que reprovou a rejeição por score
   (perdia 26,7% dos TAS experimentais).

O `--nd` ser obrigatório lá diz algo sobre a expectativa do campo: rodar sem filtrar ncRNA
não é tratado como opção.

**FEITO E MEDIDO (2026-08-08).** `--ncrna-db` no `mapping` (`553cf28`, `bfa3a9e`), índice do
Rfam por `extras/fetch_ncrna_db.sh` (`d717a2b`), e as duas corridas refeitas ponta a ponta.

**O índice:** 829 famílias tipadas rRNA/tRNA/snRNA/snoRNA no `family.txt` do Rfam, 7,3 M
sequências, 1,8 GB de FASTA, `bowtie-build` de ~1h45 em 8 threads — cabe em 8 GB de RAM sem
deduplicar (37,7% são cópias exatas; a dedup encolheria, não foi necessária).

*A seleção estava errada e foi consertada antes de medir:* casar as palavras na **linha
inteira** do `family.txt` trazia 63 famílias que só as mencionam na descrição (RNase P,
tmRNA, telomerase, íntrons de grupo I/II, `Histone3`, `Tymo_tRNA-like`) e derrubava 5
snoRNAs de Giardia chamados `GlsR*_miR*`, pelo filtro anti-miRNA. Selecionar pelo campo de
tipo dá exatamente os 829 já registrados aqui.

**Reads removidas:** 41,5% em Arabidopsis (9,99 M de 24,1 M) e 31,0% em Drosophila
(3,41 M de 11,0 M). Da massa **mapeada**, saíram 39,2% e 32,9%.

**Denominador casado é obrigatório para ler isto.** `--mincov` é RPM, e a biblioteca encolhe
~35%, então a corrida nativa cai de profundidade 20→12 (planta) e 10→6 (mosca): 42.106 loci
contra 29.872, span mediano 594→556, reads/locus 121→76. **Nada disso é o filtro** — com
`--norm` na biblioteca original os loci voltam ao que eram (span 594→589, reads/locus
121→120, cobertura 6,13×→6,06×, 2,8% menos loci). O que muda é o que sobrou dentro deles.

| (denominador casado) | Arabidopsis antes → depois | Drosophila antes → depois |
|---|---|---|
| massa de reads em loci estruturais | 36,49% → **3,48%** | 25,73% → **4,39%** |
| loci sobre RNA estrutural | 676 (2,26%) → **332 (1,14%)** | 413 (1,56%) → **321 (1,25%)** |
| contaminação do phasiRNA | 16,11% → **8,96%** | — |
| contaminação do siRNA | 11,13% → **6,67%** | 3,88% → **3,27%** |
| review | 4,5% → **3,4%** | 13,9% → **12,9%** |
| **piso** | TAS 6/6, miRBase 244/244 → **6/6, 242/242** | piC 85/90, miRBase 134/134 → **85/90, 134/134** |

**A massa cai 10× mas a contaminação de classe só cai pela metade**, e isso é informativo: o
que sobrepõe RNA estrutural depois do filtro são 332 loci de massa pequena. Ou são reads que
não casam no Rfam a uma mismatch (variante específica da espécie), ou são a vizinhança
transcrita desses genes — que o próprio validador avisa não ser erro. O alvo declarado era a
massa, e essa foi.

**CUIDADO com o enriquecimento, de novo.** Nos piC de Drosophila ele sobe piRNA 5,1×→7,2× e
siRNA 4,2×→6,3×. **Os dois subindo juntos é assinatura de denominador**, exatamente como
quando tiramos o miRNA das classes: saiu um quarto da massa das classes, não melhorou a
discriminação. Ganho limpo é a linha do piso e a da massa estrutural.

**Nota operacional:** o índice em `/mnt/d` (drvfs) faz `bowtie-inspect` travar; copiado para
disco local o filtro leva ~6 min em 24 M de reads com 8 threads, das quais ~30 s são carga do
índice.

Com isso a predição de estrutura (A2b) perde a motivação: ela existia para o organismo sem
GFF, e o Rfam já atende esse caso sem depender de anotação nenhuma.

**A2b. Predição de RNA estrutural quando não há GFF anotado — ENCERRADO sem fazer (2026-08-08)**

O filtro de read pelo Rfam atende o organismo sem anotação, que era a única motivação. O que
segue fica como registro do que se avaliaria, caso a via de sequência falhe em algum clado.

O `--exclude-loci` precisa de um BED. Nos dois organismos de teste ele sai de graça do GFF do
RefSeq (Arabidopsis: 684 tRNA, 287 snoRNA, 82 snRNA, 14 rRNA), mas o alvo da tool inclui
**organismo não-modelo sem anotação** — arroz, Tribolium, as 51 espécies do pirnaclusterdb.

**Desenho pretendido:** o GFF é **opcional**; quando presente, a predição fica **desativada**
(não faz sentido predizer o que já está anotado). Só quando falta é que a ferramenta externa
roda. Mesmo padrão do `--mirna-loci`, que aceita anotação curada ou predição.

**Candidatas a avaliar** (nenhuma instalada no env hoje):

| ferramenta | cobre | nota |
|---|---|---|
| **Infernal `cmscan` + Rfam** | tRNA, rRNA, snoRNA, snRNA, miRNA, riboswitch — **tudo num passo** | modelos de covariância (sequência + estrutura); Rfam completo tem alguns GB e num genoma de planta leva horas |
| **tRNAscan-SE** | tRNA | padrão de fato; rápido |
| **barrnap** | rRNA | HMM, muito rápido |
| **ARAGORN** | tRNA, tmRNA | alternativa leve |

**O que medir antes de adotar qualquer uma:** (a) quanto do RNA estrutural anotado no GFF de
Arabidopsis/Drosophila cada ferramenta recupera — se não recuperar o que já sabemos estar lá,
não serve para onde não sabemos; (b) o custo em tempo e disco num genoma de planta; (c) se a
união tRNAscan+barrnap cobre o suficiente sem o peso do Rfam completo.

**Só adicionar dependência se o teste justificar.** O `--exclude-loci` funciona com BED de
qualquer origem, então a tool não precisa saber quem o produziu — e um extrator de GFF→BED em
`extras/` (três linhas de awk, sem dependência) atende todo caso anotado.

**A1. RNA estrutural contamina as classes, e falta a categoria "não classificado" (2026-08-07)**

Primeiro teste de **especificidade** do projeto — todos os outros validadores medem recall.

| | Arabidopsis | Drosophila |
|---|---|---|
| loci sobre RNA estrutural | 676 (2,26%) | 413 (1,56%) |
| **massa de reads deles** | **36,5%** | **25,7%** |
| contaminação de classe | phasiRNA **16,1%**, siRNA 11,1% | siRNA 3,9%, piRNA 1,1% |
| padrão dominante | 377 loci de tRNA → hcsiRNA | 153 loci de tRNA → **piRNA** |

Não é ruído aleatório: metades de tRNA (tRFs) têm 30–35 nt e caem dentro ou colado na janela
24–31 do piRNA; em planta, fragmentos de 24 nt caem no hcsiRNA. **tRF é classe biológica
real** que a tool não tem onde colocar, então ela cai na classe cujo tamanho bate.

Em planta a contaminação de classe é séria (16% do phasiRNA, logo depois de zerarmos a de
miRNA); em mosca é pequena (1–4%), mas a massa de reads é grande nos dois — o que distorce
qualquer estatística de abundância por classe.

**Falta uma categoria de rejeição.** Hoje a classificação é **fechada**: o `RuleScorer`
sempre escolhe a classe de maior score, e `review` significa "as classes competem", não
"nenhuma serve". Um tRF não é ambíguo entre siRNA e piRNA — não é nenhum dos dois, e receber
rótulo afirmado é dizer algo falso. São conceitos distintos e não devem ser fundidos.

O instrumento existiria (`score` mede se a evidência definidora está presente, e os loci
estruturais têm score mediano **0,00** contra 7,61/58,38 dos demais), **mas não pode ser
usado ainda**: com limiar 0,5 o Drosophila perde 3 piC e a Arabidopsis perde **4 dos 6 TAS**;
com 2, perde os 6. Causa em A0b — 4 dos 13 loci em TAS têm score 0 porque o p-valor de fase
é insensível. É a mesma armadilha do `min_evidence` (removido em 2026-08-02): limiar global
sobre métrica que não vale igual para todas as classes derruba a classe pior servida.

**Ordem correta:** consertar A0b → então `unclassified` por score vira viável.

**Opções para o RNA estrutural enquanto isso:**
- **Máscara** (`--exclude-loci BED`): remover as regiões do clustering, como o `--mirna-loci`
  faz o recorte. Simples e reusa `carve_out`. Custo: descarta sinal real que porventura
  exista ali.
- **Coluna diagnóstica** (`structural_frac`): reportar sem decidir, deixando o filtro para o
  usuário. Conservador e barato; não limpa as classes.
- **Classe tRF própria** no `rules.yaml`. Mais fiel à biologia, mais trabalho, e precisa de
  referência para validar.
- Nota: parte do overlap é esperada (esses genes ficam em vizinhança transcrita), então
  nenhum limiar aqui conta como "passar" — os números servem para comparar corridas.

**B. Features maiores, cada uma com spec própria**

4. **Módulo pós-miRDeep** (item 7 do Tier 3): curar `result_*.csv`, isomiRs, alimentar
   `--mirna-loci`.
5. **`MLDecider`** (item 9): RandomForest/XGBoost sobre o `.clusters.tsv`, SHAP. Precisa de
   loci rotulados.

**C. Encerrados nesta rodada** — alocação única de multi-read (descartada com evidência,
seção "Multimapeamento"); PHAS do sRNAanno (feito, `extras/validate_phas.py`); DBSCAN
(item 10) perdeu a motivação: existia para o caso de as ilhas fragmentarem clusters
grandes, que era o problema resolvido em 2026-08-06.

### Pendências menores

- `git push` (chave `~/.ssh/id_ed25519` tem passphrase; precisa de `ssh-add` num terminal seu).
- Minors diferidos no ledger: guarda de `evidence.density_score` (ordem dos checks com
  `lib_size=None`), e o teste de threads não empurra `genome_bp` não-zero pelo pool.
- **A maioria dos loci tem duas populações de tamanho.** O número "18%" registrado aqui era
  de `pad=75` e está obsoleto: com os loci atuais a `size_purity` mediana é **0,190** em
  Drosophila e **0,647** em Arabidopsis, ou seja, em mosca a maioria fica abaixo de 0,5. Isso
  não é degradação — loci maiores agregam mais populações — mas invalida qualquer limiar
  calibrado sobre o 18%. A tarefa 7
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
