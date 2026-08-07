# Levantamento: o que o campo usa para cada assinatura

**Autor:** Lucas
**Data:** 2026-08-07
**Estado:** parcial — 3 de 5 assinaturas fechadas com fonte primária. Ver "Pendente".

## Por que este documento existe

As pontuações de evidência do `classify` foram construídas seguindo uma fonte só
(proTRAC) e engenharia reversa de código que por acaso estava em disco. O resultado foi
descobrir problemas por tentativa e erro, e propagar pelo menos uma afirmação herdada sem
verificar. Este documento fixa, assinatura por assinatura: **qual estatística o campo usa,
contra qual nulo, e o que se sabe da limitação** — com fonte.

Regra adotada: nada entra aqui sem fonte primária ou código de referência lido. Onde não
houver, fica marcado como pendente, e não como "provavelmente é assim".

---

## 1. Ping-pong (sobreposição 5'-5' de 10 nt)

**Origem biológica.** A assinatura vem do ciclo de amplificação Aub/Ago3: o corte guiado
por uma piRNA gera a extremidade 5' da parceira, deixando sobreposição de exatamente 10 nt
entre as duas fitas (Brennecke et al., 2007; Gunawardane et al., 2007).

**O que o campo usa.** Z-score da contagem no offset 10 contra o **fundo empírico dos
demais offsets** — não contra um nulo teórico. Lido no código do PingPongPro
(`pingpongpro.cpp:695`, `calculateFDRs`): média e desvio são calculados sobre todos os
offsets **exceto o 10**, e o z é convertido em FDR por aproximação normal.

**O que ele faz e nós não:** estratifica o fundo por **altura da pilha, viés de base e
altura local** (`heightScoreBin`, `baseBiasBin`, `localHeightScoreBin`), ou seja, compara
cada sítio com o fundo de sítios **comparáveis em profundidade**. O nulo depende da
cobertura, e ele trata isso explicitamente.

**Nosso estado.**
- `signals.pingpong_zscore` — média/desvio dos offsets 1–20 exceto o 10. **Segue o padrão.**
- `evidence.pingpong_score` — `binomtest(k, n, 1/20)`, ou seja, **nulo uniforme sobre os
  offsets**. Isto é invenção nossa e contraria a escolha da referência, que usa fundo
  empírico justamente porque os offsets **não** são uniformes.
- Não estratificamos por profundidade.

**Ação:** substituir o nulo uniforme pelo fundo empírico do próprio locus (que já
calculamos para o z-score). Verificável contra os landmarks flamenco × 42AB.

---

## 2. Viés de base 5' (1U, 5'A)

**Origem biológica.** Piwi/Aub carregam preferencialmente piRNAs com 5'U; Ago3 carrega as
com 10A, consequência do corte (Brennecke et al., 2007). Em planta, a primeira base
classifica em Argonauta: 5'U → AGO1, 5'A → AGO4 (via RdDM).

**O que o campo usa.** O proTRAC (Rosenkranz & Zischler, 2012) testa a frequência da base
5' contra composição aleatória. **Não verifiquei no código dele qual valor exato usa** — o
código não está em `ref_tools/` e a afirmação "nulo 0,25, como o proTRAC" no nosso ROADMAP
é herdada e **não conferida**.

**Nosso estado.** `evidence.base_bias_score` usa `p0 = 0.25`, `n = distinct_seqs` do locus,
`k = pct_1u × distinct_seqs`. O escopo está certo (só reads do locus — `features.py:48`).

**O que foi medido (2026-08-07)** e o que **não** se conclui dele: a composição das reads é
T ≈ 0,307 (Drosophila) e 0,286 (Arabidopsis), não 0,25. **Isto não prova que o nulo deva ser
esse** — as reads são o produto da seleção que se quer detectar, então usar a composição
delas é circular. Um nulo não-circular seria a composição **genômica local** do locus, que
não é produto da seleção. Fica como hipótese a testar, não como conserto.

**Ação:** ler o proTRAC (obter o código ou o texto do método) antes de mexer.

---

## 3. Fase (phasing)

**O que o campo usa.** Hipergeométrico sobre **posições**, com a abundância entrando apenas
como filtro posterior. Lido no código:
- **PhasiHunter** (`bin/phaseFunction.py:472`): `total_n` = posições ocupadas na janela,
  `total_k` = posições em registro; abundância só no corte `min_read_num`.
- **PHASIS** (`phasclust.genome.v2.pl`): guarda por registro `(abun, n, k, ..., pvalue)` —
  mesma estrutura de contagem de posições.
- **Nenhum dos dois corrige a varredura de janelas.**

**Nosso estado.** Mesma família, **mais** a correção de Bonferroni sobre as janelas
varridas (`0aa8bae`). Ou seja: mais correto que a referência e menos sensível.

**Medido:** contra os 8 loci TAS experimentais, detectamos **1/8** a p ≤ 0,01. Sem
correção, o nulo por permutação dá 76% de falso positivo em loci de 120 ciclos.
A concentração de massa no registro dominante dos TAS é ~49% contra 4,8% do acaso — o sinal
é forte e o teste por posição não o enxerga.

**Divergência não resolvida:** o sRNAminer (que gerou a anotação do sRNAanno) reporta
`Phas_Ratio` (abundância) e p de 1e-9 a 1e-15 nos mesmos loci. Não temos o código nem o
método publicado em mãos. **Comparar nossa detecção com a anotação dele mede quão liberal
ele é, não quem está certo.**

**Ação:** obter o método do sRNAminer antes de redesenhar. O conjunto experimental
(`extras/phas_experimental_at_TAIR10.bed`) é a régua, não a anotação predita.

---

## 4. Viés de fita — PENDENTE

Nosso `s_fractop` (binomial bicaudal vs 0,5) e `s_fractop_bid` (fita minoritária vs nulo
0,10) não foram conferidos contra fonte. O 0,10 é convenção nossa. Falta levantar o que se
usa para separar cluster uni-fita (somático, tipo flamenco) de bifita (germinativo, tipo
42AB).

## 5. Densidade — PENDENTE

`evidence.density_score` (Poisson contra distribuição uniforme no genoma) já está marcada
diagnostic-only no README. Falta levantar se o campo usa algo equivalente e para quê.

---

## Problema transversal, independente de assinatura

Todas as pontuações são `log10(1/p)`, e **p confunde tamanho de efeito com tamanho
amostral**. Medido no piRNA: correlação do `score` com `distinct_seqs` = **0,598**, com o
efeito (`pct_1u`) = **0,398**. Consequência: os loci de score mais alto **não** são os mais
confiáveis — dentro de piC conhecido, a faixa 10–50 tem 38,1% e a faixa 50+ tem 25,5%,
abaixo do fundo de 28,3%.

Para que "maior = mais confiável" seja verdade, a pontuação precisa ser **tamanho de efeito
com incerteza** (p.ex. limite inferior do IC), não significância. Isso vale para todas as
assinaturas e é independente de consertar cada nulo.

---

## Ordem de trabalho decorrente

1. Fechar as duas assinaturas pendentes (fita, densidade) e o método do proTRAC e do
   sRNAminer.
2. Corrigir cada assinatura conforme o levantamento, uma a uma, validando contra as
   referências que já existem: `phas_experimental_at_TAIR10.bed` (8 TAS), `piC_ref_dm6.bed`
   (90 piC), miRBase (`dme`/`ath`), e o nulo por permutação.
3. Só então redesenhar a pontuação de classe sobre tamanho de efeito.

## Referências

Brennecke, J., Aravin, A. A., Stark, A., Dus, M., Kellis, M., Sachidanandam, R., & Hannon, G. J. (2007). Discrete small RNA-generating loci as master regulators of transposon activity in *Drosophila*. *Cell*, 128(6), 1089–1103. https://doi.org/10.1016/j.cell.2007.01.043

Chen, C., Li, J., Feng, J., Liu, B., Feng, L., Yu, X., Li, G., Zhai, J., Meyers, B. C., & Xia, R. (2021). sRNAanno — a database repository of uniformly annotated small RNAs in plants. *Horticulture Research*, 8, 45. https://doi.org/10.1038/s41438-021-00480-8

Rosenkranz, D., & Zischler, H. (2012). proTRAC — a software for probabilistic piRNA cluster detection, visualization and analysis. *BMC Bioinformatics*, 13, 5. https://doi.org/10.1186/1471-2105-13-5

Uhrig, S., & Klein, H. (2019). PingPongPro: a tool for the detection of piRNA-mediated transposon-silencing activity. *Bioinformatics*, 35(2), 335–336. https://doi.org/10.1093/bioinformatics/bty578
