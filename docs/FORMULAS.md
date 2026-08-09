# Fórmulas do `classify` — referência

**Autor:** Lucas
**Atualizado:** 2026-08-08

Toda quantidade que o `classify` reporta, com a fórmula, o nulo e onde ela é usada.
Comparação com o que o campo faz: `docs/superpowers/specs/2026-08-07-signature-literature-review.md`.

**Convenção que atravessa tudo:** cada read multimapeadora contribui `w = 1/NH`, e **todas**
as features são calculadas **só com as reads do próprio locus** (`features.py:48` e `:129`).
A biblioteca inteira entra apenas em `reads_rpm`.

---

## 1. Definição de locus (clustering)

| passo | fórmula | onde |
|---|---|---|
| profundidade | `depth[p] = Σ 1/NH` das reads que cobrem `p` | `features.chromosome_depth` |
| limiar | `min_depth = max(1, round(mincov_RPM × lib_size / 1e6))` | `classify.mincov_to_depth` |
| ilha | posições contíguas com `depth ≥ min_depth` | `clustering.call_islands` (pad 0) |
| extensão | cada ilha cresce `pad` **dos dois lados** → vão tolerado `2 × pad` | `clustering.merge_compatible` |
| merge | funde vizinhos só se o **perfil for compatível** | idem |
| disjunção | vizinhos que não fundem **dividem o vão no meio** | idem |
| recorte | intervalos de `--mirna-loci` viram loci próprios | `clustering.carve_out` |

**Compatibilidade de perfil** (`clustering._compatible`), duas condições:
- **tamanho:** a janela de classe que concentra **> 50% da massa** de reads da ilha tem que
  ser a mesma. Ilha sem janela dominante não bloqueia. *(A moda do tamanho foi abandonada: com
  ~5 reads por ilha ela vira ruído.)*
- **fita:** `|frac_top_a − frac_top_b| ≤ 0,5`.

`lib_size` = **Σ 1/NH** de toda a biblioteca (`samio.weighted_library_size`), não o número de
linhas de alinhamento.

---

## 2. Features por locus

```
reads        = Σ 1/NH                       (soma das duas fitas)
reads_rpm    = reads / lib_size × 1e6
distinct_seqs= nº de sequências distintas   ← é o "n" de todos os testes binomiais
frac_top     = Σ 1/NH(+) / reads
unique_frac  = Σ 1/NH das reads com NH=1 / reads
multimap_frac= nº de linhas com NH>1 / nº de linhas          (diagnóstico bruto)
```

**Tamanho** (reads orientadas: fita − é revcomp):
```
seq_lens[L]   = Σ 1/NH das reads de comprimento L, para dicermin ≤ L ≤ dicermax
frac_in_window= (reads − short − long) / reads
dicer_call    = argmax_L seq_lens[L]   se frac_in_window ≥ 0,8   senão "N"
L1            = argmax_L seq_lens[L]
banda(L)      = seq_lens[L−1] + seq_lens[L] + seq_lens[L+1]
size_purity   = banda(L1) / Σ seq_lens
L2            = argmax seq_lens[L] entre os L com |L − L1| ≥ 2
dicer_call_2  = L2, reportado só se banda(L2) / Σ seq_lens ≥ 0,2
```
**Por que banda e não comprimento único:** o produto de Dicer se espalha ±1 nt, então a
pureza sobre um comprimento só chamava de impuro um locus 21 nt limpo — mediana 0,196 em
Drosophila — e o `dicer_call_2` disparava em 98,7% dos loci, metade deles nomeando um
comprimento **vizinho**, que é a dispersão da mesma população e não uma segunda via.

**Base 5'** (`pct_1u`, `pct_5a`) e `pct_10a`: fração ponderada de reads cuja base 1 (ou 10)
é a indicada, lidas na orientação da read.

**Espectro de offsets 5'-5'** (`signals.offset_spectrum`), base do ping-pong:
```
counts[d] = Σ_p  plus5[p] × minus5[p + d − 1]      d = 1..dmax
```
`plus5[p]` = massa cujo 5' está em `p` na fita +; `minus5` idem na fita −.
```
pingpong_frac    = counts[10] / Σ_{d=1..20} counts[d]      (só se Σ ≥ 10 pares)
overlap_dup_frac = counts[dicer_call − 2] / Σ_{d=1..dmax} counts[d]
```
O offset do dúplex de Dicer é `L − 2` (overhang de 2 nt em 3'), **não** um valor fixo.

**Fase** (`signals.phase_test`), só com reads de comprimento `= dicer_call`:
```
registro(read) = (5' − início_do_locus) mod L,  com 5' da fita − corrigido em −2
phase_number   = nº de ciclos ocupados no registro dominante
phase_ratio    = massa no registro dominante / massa total
phase_score    = phase_ratio × phase_number × ln(massa no registro × 20e6/lib_size)
```
**A massa é normalizada para 20 milhões** (`signals.NORM_LIB_SIZE`), como o
`PHASMinerCLI --normLibSize 20000000` do sRNAminer. Sem isso o score cresce com a
profundidade da biblioteca — medido, correlação **0,572** com `log(reads)`: mediana 0,89 em
loci com <50 reads contra 10,88 acima de 5.000 — e o `phase_x0` do sigmoide vira um limiar da
biblioteca em que foi calibrado. O 20e6 é o valor deles e também é onde nossa calibração está
(a biblioteca de planta tem 19,93 M, fator 1,003), então os limiares publicados seguem valendo;
em Drosophila, que tem 9,66 M, o fator é 2,07.

---

## 3. Pontuações de evidência (`s_*`)

Todas são `score = −log10(p)`, limitadas a 300 (`evidence.score_from_p`).

| coluna | teste | nulo |
|---|---|---|
| `s_1u` | binomial unilateral, `k = pct_1u × distinct_seqs`, `n = distinct_seqs` | **composição do DNA do locus** (ver abaixo); 0,25 sem `-g` |
| `s_5a` | idem para 5'A | idem |
| `s_fractop` | binomial **bicaudal** de `frac_top` | 0,5 |
| `s_fractop_bid` | binomial unilateral da **fita minoritária** | **0,10** — medido: fração minoritária mediana em flamenco (uni-fita) = 0,099 |
| `s_fractop_bal` | TOST de equivalência (`0,5 ± 0,15`) | **obsoleto**, ver nota |
| `s_pingpong` | Poisson: `P(X ≥ counts[10])` | **taxa média dos offsets 1–20 do próprio locus, excluindo o 10** |
| `s_phase` | hipergeométrico, ver abaixo | — |

**Nulo da base 5'** (`samio.base_null`) — responde *"se o corte caísse em qualquer posição
deste locus, com que frequência o produto começaria com esta base?"*:
```
null(base) = frac_top × f(base no DNA) + (1 − frac_top) × f(complemento no DNA)
```
com `f` calculado sobre **todas** as posições do locus (N excluído do denominador). A
ponderação por `frac_top` existe porque uma read da fita − começa com T onde o DNA tem A.
Medido: 0,272 (Drosophila) e 0,332 (Arabidopsis) — não 0,25.

**Fase — hipergeométrico com correção de varredura:**
```
janela      = 11 ciclos, deslizando 1 ciclo por vez sobre o locus
população   = 2 × comprimento_da_janela − 1
sucessos    = 2 × ciclos − 1
sorteios    = posições distintas ocupadas na janela
observados  = phase_number do registro dominante
p_janela    = P(X ≥ observados)
p           = min(1, p_janela × nº de janelas testadas)      ← Bonferroni
```
A janela vencedora é escolhida por **evidência** (`number`, depois `score`), nunca pelo
p-valor — o hipergeométrico premia janela vazia. A correção é obrigatória: sem ela, no nulo
por permutação, **76%** dos loci sem fase alguma passavam a p ≤ 0,01 num locus de 120 ciclos.

**Forças de efeito (`e_*`)** — o que o `score` usa. Mesmos `k`, `n` e `p₀` das `s_*`
correspondentes, lidos como tamanho de efeito em vez de significância:

```
e_X = max(0, (L95(p̂, n) − p₀) / (1 − p₀))
```

onde `L95` é o limite inferior de Wilson a 95%. Fica em `[0,1]`: **0,61 significa que o
efeito está a 61% do caminho entre o nulo e o máximo, com 95% de confiança**. Não são
colunas do TSV — vivem no `ClusterRecord` e alimentam o `score`.

| força | p̂ | n | p₀ |
|---|---|---|---|
| `e_1u` | `pct_1u` | `distinct_seqs` | composição do DNA (T) |
| `e_5a` | `pct_5a` | `distinct_seqs` | composição do DNA (A) |
| `e_fractop_bid` | fita minoritária | `distinct_seqs` | 0,10 |
| `e_pingpong` | `counts[10] / Σcounts` | `Σcounts` | taxa de fundo / `Σcounts` |
| `e_phase` | `ciclos em registro / 11` | `11` (a janela) | `1 / phase_length` |

---

## 4. Classificação (`decider.RuleScorer`)

O score da classe **não** usa as colunas `s_*` — usa as features cruas:

```
_window(size, lo, hi, edge=3) = 1 dentro da janela; decai linearmente 1 − d/edge fora
_mid(frac_top, σ=0.15)        = exp(−(frac_top − 0,5)² / 2σ²)
_sigmoid(x, x0, scale)        = 1 / (1 + exp(−(x − x0)/scale))

S(classe) = w_size          × _window(dicer_call, lo, hi)
          + w_1u            × pct_1u
          + w_5a            × pct_5a
          + w_fractop_mid   × _mid(frac_top)
          + w_pingpong      × pingpong_frac × gate
          − w_pingpong_penalty × pingpong_frac × gate
          + w_phase         × _sigmoid(phase_score, phase_x0, phase_scale) × gate
```
`gate = _window(dicer_call, lo, hi)` da **própria classe** — ping-pong e fase só pontuam a
classe cujo tamanho já encaixa. É isso que permite um único código servir planta e artrópode.

```
cls        = argmax_classe S
confidence = softmax(S)[cls]
margin     = S(1º) − S(2º)
review     = confidence < review_confidence  OU  margin < min_margin
             OU  reads < min_reads  OU  dicer_call == "N"
score      = min( e_* das evidências em defining_evidence )     ← 0 se dicer_call == "N"
             ← NA em miRNA: a classe vem da anotação, não foi pontuada
e_X        = max(0, (L95(p̂, n) − p₀) / (1 − p₀))              ← L95 = limite inferior de Wilson
```
Em `defining_evidence`, uma lista aninhada é **OU** (usa o máximo); entradas no topo são
**E** (usa o mínimo). Ex.: `[[1u, pingpong]]` = "1U **ou** ping-pong".

---

## 5. Limitações conhecidas de cada fórmula

- **O `score` ordena confiabilidade desde 2026-08-08**, porque deixou de ser `−log10(p)`.
  Medido: correlação com o efeito 0,323 → **0,883**, com `distinct_seqs` 0,629 → **0,137**, e
  a evidência definidora sobe monotonicamente com a faixa de score nos dois reinos (1U
  0,38→0,82 em Drosophila; 5'A 0,43→0,86 em Arabidopsis).
  *Escalas diferem por classe e isso é honesto:* mediana 0,47 em Drosophila e **0,127** em
  Arabidopsis, porque o `hcsiRNA` é definido só por `5a` (0,43 contra um nulo de 0,332 deixa
  pouco espaço) enquanto o piRNA tem a rota `1u OU pingpong`. O score está dizendo que a
  evidência de hc-siRNA é mais fraca — o que é verdade, e a escala de p-valor escondia.
- **As `s_*` seguem sendo p-valores** e respondem outra pergunta ("isto surpreenderia sob o
  nulo?"). Não são o que o `score` usa.
- **`s_fractop_bal` está obsoleto.** Media simetria (50/50) quando a assinatura é
  bidirecionalidade, e não tinha poder: um locus perfeitamente equilibrado com n=50 tirava
  1,68, abaixo do 2 que deveria passar. Substituído por `s_fractop_bid` na
  `defining_evidence`; a coluna segue no TSV só por compatibilidade.
- **`s_fractop_bid` é ruidoso.** Metade dos loci de flamenco (uni-fita) passa do nulo de 0,10,
  e a separação uni-fita × bifita fica em 1,5–1,7× em qualquer limiar.
- **`s_phase` é pouco sensível.** 1 de 8 loci TAS experimentais a p ≤ 0,01. Medido contra a
  fórmula do sRNAminer nos **mesmos dados**: ele detecta 2 de 8, com 10× mais falso positivo
  (ele não corrige a seleção do registro). Não há caso para redesenhar *o teste* — a correção
  pela varredura de janelas está certa e é mais conservadora que o campo. O que não funciona
  é usar esse p-valor para **ranquear** loci entre si (AUC 0,560/0,481, no acaso — ver
  `e_phase` abaixo): esse trabalho passou para `e_phase`/`score`, que usa a mesma contagem
  de ciclos sem convertê-la em p-valor.
- **`s_pingpong` assume variância = média** (Poisson). Medido no fundo de offsets: mediana
  1,3, mas 22,5% dos loci acima de 5. Binomial negativa seria a síntese.
- **`phase_score` correlaciona 0,53 com `log(reads)`** — mede abundância junto com fase, então
  não é comparável entre loci. `phase_pvalue` também não serve para isso (AUC 0,560/0,481,
  no acaso — ver `e_phase` abaixo); para comparar loci entre si use `e_phase`/`score`.
- **`pingpong_z` é diagnóstico**, capado em 100 e instável em locus estreito e denso. Para
  classificar, use `pingpong_frac`.
- **`overlap_dup_frac`, `size_purity`, `dicer_call_2`, `unique_frac`** são reportados e **não
  entram em regra nenhuma**.
- **`s_density` foi removido** (2026-08-07): correlacionava 0,72 com `log(reads)` e era alto
  por construção, já que o locus só existe depois de passar no `mincov`, que é um corte de
  densidade. `evidence.density_score` continua no código, sem estar ligado.
- **A escala mudou**: um `score` de 50 na escala antiga vira ~0,7 na escala de effect strength,
  então corridas anteriores a 2026-08-08 não comparam com as novas.
- **`e_phase` mede a extensão da fase, não a surpresa.** É a fração dos 11 ciclos da janela que
  está no registro dominante — a corrida de cortes processivos é a assinatura do mecanismo.
  Medido contra controles casados por `distinct_seqs`: AUC 0,826 (TAS experimental) e 0,769
  (PHAS anotado), contra 0,560/0,481 do p-valor dos mesmos dados e 0,415/0,293 da razão de
  abundância. Ressalvas: a janela é escolhida por maximizar evidência, e o `p₀ = 1/phase_length`
  é anticonservador porque o registro reportado é o máximo entre 21.
- **`e_phase` satura com ocupação total, não com concentração num registro.** Ele só lê
  `phase_number` — quantos ciclos têm *alguma* leitura no registro dominante — e ignora
  `phase_ratio`. Um locus com uma leitura de 21 nt em **cada** posição de uma janela de 250 nt
  dá `phase_number = 11` (o máximo), mesmo com `phase_ratio = 0,048` (= 1/21, nenhum
  enriquecimento) e `phase_pvalue = 1,0`: `e_phase = 0,7282`, o teto. Em Arabidopsis, 8 loci
  caem nesse regime (`phase_number >= 11` e `phase_ratio < 0,15`), todos `hcsiRNA` de 24 nt
  definidos por `5a` — nada errado nesta corrida. Mas em `phasiRNA`, `score` **é** `e_phase`:
  um locus de 21 nt densamente sequenciado e não fásico, noutra biblioteca, carregaria
  evidência de fase máxima com p-valor 1. Teste:
  `tests/test_features_phase_evidence.py::test_full_occupancy_saturates_the_statistic_without_enrichment`.
- **`score == 0` num locus `phasiRNA` significa `phase_number <= 1` na melhor janela**, não
  "não foi procurada evidência de fase". `e_phase` é uma função-degrau de 12 níveis de
  `phase_number` (0 a 11): `0, 0,004, 0,052, 0,109, 0,173, 0,244, 0,322, 0,406, 0,499, 0,604,
  0,729`. Em Arabidopsis (`at_ephase.clusters.tsv`), 63% dos loci chamados `phasiRNA` têm
  `score == 0` (era 56% antes desta mudança).
- **`e_pingpong` é otimista**: o `n` são pares, que escalam com o quadrado da cobertura e não
  são ensaios independentes.
- **`e_fractop_bid` herda o ruído do `s_fractop_bid`**: a fração minoritária separa uni-fita de
  bifita por apenas 1,5–1,7× no nível do locus.
