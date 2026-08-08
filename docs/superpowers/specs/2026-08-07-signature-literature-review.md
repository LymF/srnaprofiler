# Levantamento: o que o campo usa para cada assinatura

**Autor:** Lucas
**Data:** 2026-08-07
**Estado:** completo — 5 de 5 assinaturas fechadas contra código de referência lido.

## Por que este documento existe

As pontuações de evidência do `classify` foram construídas seguindo uma fonte só
(proTRAC) e engenharia reversa de código que por acaso estava em disco. O resultado foi
descobrir problemas por tentativa e erro, e propagar pelo menos uma afirmação herdada sem
verificar. Este documento fixa, assinatura por assinatura: **qual estatística o campo usa,
contra qual nulo, e o que se sabe da limitação** — com fonte.

Regra adotada: nada entra aqui sem fonte primária ou código de referência lido. Onde não
houver, fica marcado como pendente, e não como "provavelmente é assim".

## Resumo — o que o campo faz, e o que é nosso

| assinatura | o campo | nós |
|---|---|---|
| ping-pong | z contra **fundo empírico** dos outros offsets, estratificado por profundidade (PingPongPro) | z igual ✓ **+** binomial com nulo uniforme 1/20 (**nosso**) |
| base 5' | **limiar de fração** 0,75 (ou 0,5 em ambos) — sem teste, sem nulo (proTRAC) | binomial contra **0,25** (**nosso**, atribuído falsamente ao proTRAC) |
| fase | hipergeométrico sobre posições, **sem** correção de varredura (PhasiHunter, PHASIS) | igual **+ correção de Bonferroni** (mais correto que o campo) |
| viés de fita | **limiar de fração** (proTRAC `-clstrand`, ShortStack `--strand_cutoff`) | limiar igual ✓ **+** binomial vs 0,5 e vs 0,10 (**nossos**) |
| densidade | percentil **empírico** das janelas da própria biblioteca (proTRAC `-pdens`) | Poisson contra **uniforme** (**nosso**) |

**Padrão que emerge:** o campo usa **limiares de fração e distribuições empíricas**; nós
trocamos quase tudo por **testes contra nulos teóricos**. Em três das cinco assinaturas o
nulo teórico é invenção nossa (0,25; 1/20; uniforme para densidade), e em duas delas a
referência **explicitamente evita** o nulo teórico porque os dados não são uniformes.

Isso não torna nossa escolha errada por si só — testar contra nulo é mais informativo que
cortar por limiar, quando o nulo está certo. Mas explica por que fomos descobrindo
miscalibração atrás de miscalibração: **estávamos assumindo uniformidade em dados que o
campo já sabe não serem uniformes.**

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

**O que o campo usa — VERIFICADO no código** (`proTRAC_2.3.1.pl:201-211`, obtido via
bioconda; o domínio uni-mainz saiu do ar):

```
-1Tor10A [0..1]   Fração das contagens (normalizadas) com 1T (1U) OU 10A  (default 0.75)
-1Tand10A [0..1]  Se ficar abaixo, aceita a janela se AMBOS 1U E 10A      (default 0.5)
```

**O proTRAC não testa a base 5' contra nulo nenhum.** Ele aplica **limiares de fração
fixos**. Não há binomial, não há 0,25, não há p-valor para base.

**Portanto a afirmação do nosso ROADMAP — "Nulo da base 5' = 0,25, como o proTRAC" — é
FALSA.** O 0,25 é invenção nossa, atribuída incorretamente. Isso foi propagado desde a spec
de 2026-08-02 e repetido por mim nesta sessão sem conferência.

**Nosso estado.** `evidence.base_bias_score` usa `p0 = 0.25`, `n = distinct_seqs` do locus,
`k = pct_1u × distinct_seqs`. O escopo está certo (só reads do locus — `features.py:48`).

**O que foi medido (2026-08-07)** e o que **não** se conclui dele: a composição das reads é
T ≈ 0,307 (Drosophila) e 0,286 (Arabidopsis), não 0,25. **Isto não prova que o nulo deva ser
esse** — as reads são o produto da seleção que se quer detectar, então usar a composição
delas é circular. Um nulo não-circular seria a composição **genômica local** do locus, que
não é produto da seleção. Fica como hipótese a testar, não como conserto.

**Ação:** decidir conscientemente qual nulo queremos, já que não há padrão a seguir. As
opções reais são (a) limiar de fração, como o proTRAC; (b) binomial contra composição
genômica **local** do locus (não circular, ao contrário da composição das reads); (c) manter
0,25 assumindo composição uniforme, mas documentando que é escolha nossa e não herança.

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

**sRNAminer — RESOLVIDO (2026-08-07, jar descompilado com CFR).**
`QuickStatPHASRegion.StatRegion:154`:

```java
pvalue = HyperGeometric(2*locusLen-1,            // N = população
                        (2*locusLen-1)/phasLen,  // M = sucessos
                        posStrandCount.size(),   // n = sorteios
                        k);                      // observados
```

- **Também conta posições, não abundância.** `k` vem de `phasSpeciesCount`, que só incrementa
  quando o par (fita, 5') é novo — deduplicado por `seenStringHits`. A abundância
  (`inPhasRatio`) é calculada e reportada, mas **fica fora do p-valor**; entra como limiar
  separado (`minPhasRatio`).
- **Não varre janelas:** usa o locus inteiro numa tacada, logo não há seleção de máximo e
  não precisa de correção. É a única diferença estrutural para nós.
- Soma os registros vizinhos (`maxPhase±1`) ao calcular a razão — tolerância de ±1 nt que
  não temos.

**A COMPARAÇÃO QUE MOTIVOU TUDO ISTO ERA INVÁLIDA.** Rodando a fórmula dele nos **nossos**
dados, nos 8 TAS: p de 1,8e-2 a 8,4e-4, **não** os 1e-9–1e-15 da anotação. Detecção a
p ≤ 0,01: **nós 1/8, ele 2/8**. Ou seja, os p minúsculos da anotação vêm das **bibliotecas
deles**, não de estatística superior. Eu estava comparando nossos p na nossa biblioteca com
os p deles nas bibliotecas deles — o que não mede nada sobre as estatísticas.

**Conclusão — MEDIDA (2026-08-07), três variantes, nos dois reinos:**

| variante | calibração no nulo (p ≤ 0,01; esperado 1%) | TAS detectados |
|---|---|---|
| **A — nossa** (varredura de janelas + Bonferroni) | 0,7–3,7% | **1/8** |
| B — sRNAminer (locus inteiro, sem correção) | **9,3–13,0%** | 2/8 |
| C — locus inteiro + correção de registro (×L) | 0–0,7% | 0/8 |

Trade-off monotônico, sem almoço grátis. A variante B é ~10× anticonservadora e o valor é
**plano no tamanho do locus** — logo a causa não é a varredura, é a **seleção do registro
dominante entre L**, que o sRNAminer não corrige. A sensibilidade extra dele (2/8) é comprada
com 10× mais falso positivo.

**Decisão: manter a nossa.** Está no melhor ponto dos três.

Detalhe que valida o desenho original da varredura: **TAS1c dá 2,0e-2 em A e 1,0 em C** — a
varredura acha uma janela localmente fásica que o cálculo sobre o locus inteiro dilui. A
proteção contra "fase local em locus grande" não é hipotética; ela opera num TAS real.

Padrão idêntico em 21 nt (planta) e 26 nt (Drosophila), então a escolha serve os dois reinos
com um só código — que é o requisito da tool.

---

## 4. Viés de fita

**O que o campo usa — VERIFICADO.** Também **limiar de fração**, não teste:
- proTRAC (`:239`): `-clstrand [0..1]` — fração das contagens normalizadas na fita
  principal.
- ShortStack (README:123): `--strand_cutoff`, ponto flutuante entre 0,5 e 1; a
  "strandedness" do locus é inferida por comparação de `FracTop` com esse corte.

Ou seja, **nenhuma das duas faz teste estatístico de fita** — ambas classificam por corte.

**Nosso estado.** Temos as duas coisas: o corte (`strand_cutoff: 0.8`, que gera
`strand_call` e segue o padrão do campo) **e** testes próprios (`s_fractop` binomial
bicaudal vs 0,5; `s_fractop_bid` fita minoritária vs 0,10). Os testes são nossos; o 0,10 é
convenção sem respaldo levantado.

**Nulo de 0,10 — MEDIDO (2026-08-07), e valida-se a posteriori.** A fração minoritária
mediana nos loci de **flamenco** (uni-fita por biologia conhecida, `reads > 20`, n=397) é
**0,099**. O 0,10 escolhido sem base é exatamente o que um locus genuinamente uni-fita vaza
para a outra fita. Agora tem justificativa empírica, e não herdada.

**Limitação medida, e é séria:** a quantidade é **ruidosa no nível do locus**. Metade dos
loci de flamenco passa de 0,10, e a separação uni-fita × bifita fica em **1,5–1,7× em
qualquer limiar** (0,05 / 0,10 / 0,20 / 0,25). A região é uni-fita no agregado (frac_top
0,88) mas os loci individuais variam muito — 42AB tem minoritária mediana 0,240 contra 0,099
de flamenco, faixas que se sobrepõem largamente.

*Ressalva de escopo:* o `s_fractop_bid` não existe para separar flamenco de 42AB (ambos são
piRNA), e sim para atestar "as duas fitas produzem", que define o siRNA. Ainda assim o ruído
medido significa que **essa evidência é fraca por natureza**, não por implementação — o que
limita quanto o `score` do siRNA pode valer, independente de qualquer conserto.

## 5. Densidade

**O que o campo usa — VERIFICADO** (`proTRAC_2.3.1.pl:171-190`): `-pdens 0.01` significa que
a contagem normalizada da janela deve estar **no top 1% de todas as janelas deslizantes** —
percentil **empírico** contra a distribuição observada naquela biblioteca, com opção `-est`
de estimar por amostragem de janelas aleatórias de 1 kb.

**Nosso estado.** `evidence.density_score` usa **Poisson contra distribuição uniforme no
genoma** — nulo teórico onde a referência usa distribuição empírica.

**Mas o problema maior é que ela é REDUNDANTE (medido 2026-08-07), não mal calibrada.**
Correlação com `log(reads)`: **0,718** (Drosophila) e **0,729** (Arabidopsis). Ela reafirma
a contagem de reads, que já está no TSV. E é alta por construção — 71,1% e 57,7% dos loci
acima de 2 — porque o locus só existe se passou no `mincov`, que **é** um corte de densidade.
Medir depois "este locus é mais denso que o uniforme" re-mede o critério de seleção.

**Isto também esclarece a comparação com o campo:** no proTRAC a densidade é o **critério de
chamada de locus** (`-pdens`, top 1% das janelas), não uma evidência de classe. Nós usamos o
`mincov` para o mesmo fim. Já usamos densidade, no lugar certo; o `s_density` é a mesma
informação contada duas vezes.

**Não vale consertar o nulo.** O único jeito de a densidade virar informação independente é
medi-la contra o fundo **local** (a vizinhança do locus) em vez do genoma inteiro — "este
locus se destaca do entorno?". Isso não é circular com o `mincov` e distinguiria cluster real
de região amplamente transcrita. É feature nova, não conserto, e de baixa prioridade
enquanto a densidade não decidir nada em regra alguma.

---

## 6. Chamada de loci e definição de classe — como cada ferramenta faz

Descompilado do `sRNAminer.jar` (CFR) e lido do código das demais:

| | chamada de locus | classe |
|---|---|---|
| **sRNAminer** | encadeia reads por vão (`minRegionGap` 300 bp) e **filtra por abundância TOTAL da região** (`minReadDepth` 1000) | fração de abundância na janela de tamanho (`--targetLens 23,24`, `minRepeatRatio`), ≥10 reads |
| **ShortStack** | profundidade ≥ `mincov` (1 RPM) → `pad` 200 bilateral → merge | `DicerCall` por tamanho dominante, `strand_cutoff` 0,8 |
| **proTRAC** | janela deslizante, densidade no **top 1%** (`-pdens`) | limiares de fração: `-pisize` 0,75, `-1Tor10A` 0,75, `-clstrand` |
| **nós** | profundidade ≥ `mincov` (RPM) → `pad` 500 bilateral → **merge condicionado por perfil** | pesos + janelas de tamanho no `rules.yaml` |

**O que aprendemos de útil aqui:**

1. **Massa total mínima por locus** — o sRNAminer exige `totalAbundance ≥ 1000` **na etapa de
   clustering**. Nós só temos profundidade por posição (`mincov`) e depois `min_reads: 10`
   como *flag de review*, o que não impede o locus de existir. Um mínimo de massa na
   clusterização garante poder estatístico por locus **por construção** — que foi o problema
   que perseguimos a sessão inteira.
   *Adaptação:* o valor tem de ser em RPM (o 1000 deles é absoluto e não transfere) e
   **por reino**: hc-siRNA de planta é difuso por natureza (espalhado por TEs, já medimos que
   subir corte perde sinal), enquanto cluster de piRNA em Drosophila é denso. Provavelmente
   um valor no `rules.yaml`, como já fizemos com `clustering.pad`.

2. **Definir classe por fração de tamanho é o padrão universal** — as três referências fazem
   isso, e nós também (`frac_in_window`, janelas no config). Aqui estamos alinhados.

3. **O merge condicionado por perfil é nosso** — ninguém mais faz. E foi medido: fundiu 26%
   dos loci sem perder cluster conhecido. Vale manter e documentar como diferencial, não como
   desvio.

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

1. ~~Fechar as assinaturas pendentes~~ **FEITO.** Falta só o método do sRNAminer, que
   importa apenas para a fase.
2. Corrigir cada assinatura conforme o levantamento, uma a uma, validando contra as
   referências que já existem: `phas_experimental_at_TAIR10.bed` (8 TAS), `piC_ref_dm6.bed`
   (90 piC), miRBase (`dme`/`ath`), e o nulo por permutação.
3. Só então redesenhar a pontuação de classe sobre tamanho de efeito.

## Referências

Brennecke, J., Aravin, A. A., Stark, A., Dus, M., Kellis, M., Sachidanandam, R., & Hannon, G. J. (2007). Discrete small RNA-generating loci as master regulators of transposon activity in *Drosophila*. *Cell*, 128(6), 1089–1103. https://doi.org/10.1016/j.cell.2007.01.043

Chen, C., Li, J., Feng, J., Liu, B., Feng, L., Yu, X., Li, G., Zhai, J., Meyers, B. C., & Xia, R. (2021). sRNAanno — a database repository of uniformly annotated small RNAs in plants. *Horticulture Research*, 8, 45. https://doi.org/10.1038/s41438-021-00480-8

Ma, X., Zhang, X., Zhao, K., Li, F., Li, K., Ning, L., He, J., Xin, Z., Yin, D., & Li, K. (2024). sRNAminer: A multifunctional toolkit for next-generation sequencing small RNA data mining in plants. *Science Bulletin*, 69(8), 1017–1021. https://doi.org/10.1016/j.scib.2023.12.049

Rosenkranz, D., & Zischler, H. (2012). proTRAC — a software for probabilistic piRNA cluster detection, visualization and analysis. *BMC Bioinformatics*, 13, 5. https://doi.org/10.1186/1471-2105-13-5

Uhrig, S., & Klein, H. (2019). PingPongPro: a tool for the detection of piRNA-mediated transposon-silencing activity. *Bioinformatics*, 35(2), 335–336. https://doi.org/10.1093/bioinformatics/bty578
