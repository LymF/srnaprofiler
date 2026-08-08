# Evidência de fase: medir o efeito, não o p-valor

**Autor:** Lucas
**Data:** 2026-08-08
**Estado:** desenho aprovado, aguardando plano de implementação

## Problema

O `score` de um locus `phasiRNA` é a força da sua evidência definidora, e a evidência
definidora dessa classe é `phase` (`rules.plant.yaml:61`). Hoje essa força sai de um
caminho que **não separa phasiRNA real de ruído**.

Medido em Arabidopsis, com os **8 intervalos TAS experimentais — cobertos por 12 loci da
corrida** — e os 68 loci sobre PHAS21 anotado como positivos, contra controles 21-nt casados por `distinct_seqs` (o casamento é obrigatório: o
`phase_ratio` mediano do fundo é 0,954 em loci com ≤20 sequências distintas, contra 1/21
esperado — com poucas posições, o **máximo entre 21 registros** é alto por acaso):

| quantidade | AUC (12 TAS) | AUC (68 PHAS) |
|---|---|---|
| `phase_number` cru — ciclos em registro | **0,826** | **0,769** |
| p de permutação sobre `phase_number` | 0,755 | 0,744 |
| `phase_score` (razão × ciclos × ln massa) | 0,754 | **0,526** |
| p de permutação sobre `phase_score` | 0,706 | 0,613 |
| `e_phase` atual, via `score` | 0,618 | 0,639 |
| `s_phase` = −log10(`phase_pvalue`) | 0,560 | **0,481** |
| `phase_ratio` cru (o que o sRNAminer reporta) | 0,415 | 0,293 |

Três leituras, todas medidas:

1. **O p-valor atual está no acaso.** 0,560 nos TAS e 0,481 no conjunto PHAS. Não é
   conservadorismo do Bonferroni: é a formulação. Isto explica o 5,5% de detecção a p ≤ 0,01
   registrado no ROADMAP.
2. **Abundância não é o que falta.** O `phase_ratio` cru discrimina **pior que moeda**
   (0,415). O `phase_score`, que usa abundância, parece bom nos TAS (0,754) e **desaba para
   0,526** no conjunto casado — os TAS têm 1.236 reads medianas contra 173 dos controles, e
   o que ele mede é profundidade, não fase.
3. **Toda transformação do `phase_number` em p-valor perde informação.** O p mistura o efeito
   com o número de posições, exatamente como o `score` fazia antes de 2026-08-08 (correlação
   com o efeito 0,323 → 0,883 quando passou a ser tamanho de efeito).

## Decisão

`e_phase` passa a medir **quantos ciclos da janela estão em registro**, e o `score` passa a
ler isso. O `phase_pvalue` e o `s_phase` ficam **exatamente como estão**.

```
e_phase = effect_strength(k = min(phase_number, CYCLES), n = CYCLES, p0 = 1/phase_length)
```

com `CYCLES = 11`, a mesma constante da janela de `phase_test`. Em palavras: *dos 11 ciclos
que o teste enxerga, quantos estão no registro dominante* — limitado a [0,1] por construção,
com o limite inferior de Wilson dando o conservadorismo para `k` pequeno, como em toda
`e_*`.

Resultado da variante escolhida:

| | AUC TAS | AUC PHAS | mediana TAS / PHAS / controle | controles ≥ 0,3 |
|---|---|---|---|---|
| `e_phase` nova | **0,826** | **0,769** | 0,177 / 0,052 / **0,000** | **0%** (42% dos TAS) |

### Denominador: por que a constante da janela, e não o locus

Duas alternativas foram medidas e **reprovadas**:

- **Normalizar pelos ciclos do próprio locus** (`phase_number / ciclos_ocupados`): AUC cai
  para 0,647 (TAS) e 0,504 (PHAS). Os TAS têm janela de 10 ciclos contra 5 dos controles, e
  dividir pelo tamanho da janela joga fora justamente o que discrimina.
- **`p0` empírico dos outros 20 registros**, análogo ao nulo empírico do ping-pong: sai
  **invertido**. O fundo mediano é 0,201 nos TAS contra 0,050 nos controles, porque locus
  denso tem alguma cobertura em todo registro; subtraí-lo penaliza o positivo.

O denominador constante tem uma consequência deliberada: um locus que abrange 4 ciclos nunca
passa de 4/11. Isso é correto — evidência de fase é uma **corrida** de ciclos consecutivos, e
um locus curto simplesmente tem menos evidência a exibir.

`p0 = 1/phase_length` é a chance de uma posição cair no registro por acaso. **Ressalva
conhecida:** o registro reportado é o **máximo entre 21**, então o nulo verdadeiro é maior e
este `p0` é anticonservador — do mesmo tipo da ressalva de superdispersão já registrada para
o Poisson do ping-pong. O efeito prático é pequeno (nenhum controle casado alcança 0,3) e a
afirmação calibrada continua sendo do `phase_pvalue`, que não muda.

## Por que isto é o que a biologia e o campo dizem

**Mecanismo.** phasiRNA vem de corte processivo do DCL a partir de um ponto de início
definido pela clivagem guiada por miRNA — o registro é herdado desse corte inicial e mantido
ao longo do precursor (Allen et al., 2005; Fei et al., 2013). A assinatura é, portanto, uma
sequência de unidades consecutivas em registro: o comprimento da corrida. Massa concentrada
num registro é outra coisa — uma única posição quente com milhares de reads produz `ratio`
alto sem processividade nenhuma, e é por isso que ela aparece no fundo (Axtell, 2013).

**Campo.** Todos contam posições ou espécies, não abundância, no p-valor — lido do código em
`docs/superpowers/specs/2026-08-07-signature-literature-review.md` §3:

- PhasiHunter (`phaseFunction.py:472`): `total_n` = posições ocupadas, `total_k` = em registro.
- PHASIS (`phasclust.genome.v2.pl`): mesma estrutura `n`/`k` por registro.
- sRNAminer (`QuickStatPHASRegion.StatRegion:154`): `k` vem de `phasSpeciesCount`, que conta
  espécies; `inPhasRatio` é reportada mas **fica fora** do p-valor, como limiar separado.

A mudança portanto não inventa estatística nova: alinha o `score` com a evidência que o campo
inteiro usa, e deixa a abundância no papel descritivo que o próprio sRNAminer lhe dá.

## Escopo

**Muda:**
- `features.build_record`: o cálculo de `r.e_phase`.
- `docs/FORMULAS.md`: a fórmula de `e_phase` e a ressalva do `p0`.

**Não muda:**
- `phase_test`, `phase_metrics`, `hypergeometric_phase_pvalue` — nenhuma linha.
- `phase_pvalue`, `s_phase`, `phase_score`, `phase_ratio`, `phase_number` no TSV.
- `w_phase` e a decisão de classe: o `score` **reporta**, não decide. Nenhum locus deve
  trocar de classe.
- O nulo por permutação (`tests/test_signals_phase_null.py`) continua valendo sem alteração,
  porque testa `phase_test`, que não é tocado.

**Fora de escopo, registrado para depois:** exigir que os ciclos em registro sejam
**consecutivos**. A biologia da processividade diz que deveriam ser, e `phase_number` hoje
aceita ciclos espalhados pela janela. É um refinamento medível sobre os mesmos conjuntos de
verdade, não um conserto.

## Critérios de aceite

1. **Discriminação** (Arabidopsis, contra controles casados por `distinct_seqs`): AUC do
   `score` nos 12 loci sobre TAS ≥ 0,80 (hoje 0,618) e no conjunto PHAS ≥ 0,75 (hoje 0,639).
2. **Piso intacto:** TAS 6/6, miRBase 242/242 (planta); piC 85/90, miRBase 134/134 (mosca).
3. **Zero mudança de classe** nos dois reinos — o `score` não decide.
4. **`phase_pvalue` byte a byte idêntico** antes e depois, nas duas corridas. É a prova de
   que só a evidência mudou.
5. O teste de calibração por permutação continua passando sem edição.

## Riscos

- **A escala do `score` de planta muda.** A mediana de `hcsiRNA` (definido por `5a`) fica
  onde está, mas a de `phasiRNA` se move. Qualquer limiar calibrado sobre a escala antiga
  precisa ser relido — e o ROADMAP já registra que não existe limiar em produção sobre o
  `score`, exatamente porque a rejeição por score foi reprovada em 2026-08-08.
- **Verdade de referência estreita.** São 12 loci sobre TAS experimental e 68 sobre PHAS
  anotado pelo sRNAminer, numa única biblioteca de flor. O conjunto PHAS herda o quanto o
  sRNAminer é liberal; o conjunto TAS é curado mas pequeno. Nenhum dos dois autoriza calibrar
  limiar — autorizam comparar formulações, que é o que foi feito.
- **`p0` anticonservador**, já descrito acima.

## Referências

Axtell, M. J. (2013). Classification and comparison of small RNAs from plants. *Annual Review
of Plant Biology*, 64, 137–159. https://doi.org/10.1146/annurev-arplant-050312-120043

Fei, Q., Xia, R., & Meyers, B. C. (2013). Phased, secondary, small interfering RNAs in
posttranscriptional regulatory networks. *The Plant Cell*, 25(7), 2400–2415.
https://doi.org/10.1105/tpc.113.114652

Allen, E., Xie, Z., Gustafson, A. M., & Carrington, J. C. (2005). microRNA-directed phasing
during trans-acting siRNA biogenesis in plants. *Cell*, 121(2), 207–221.
https://doi.org/10.1016/j.cell.2005.04.004
