# Score por tamanho de efeito, não por significância

**Autor:** Lucas
**Data:** 2026-08-07

## O problema

O `score` do `.clusters.tsv` deveria responder *"quão confiável é esta chamada?"*. Ele não
responde, e a falha é estrutural, não de calibração.

Hoje `score = min(s_* das evidências definidoras)`, e cada `s_* = −log10(p)`. Um p-valor
mistura **tamanho de efeito** com **tamanho amostral**: com `n` grande, um viés minúsculo
produz p astronômico. Medido no piRNA de Drosophila (corrida `dros_gn2`, já com o nulo
genômico):

- correlação do `score` com `distinct_seqs` (amostra): **0,629**
- correlação do `score` com `pct_1u` (efeito): **0,323**

Ou seja, o score mede mais o tamanho da amostra do que a força do sinal. A consequência
aparece no teste de ordenação contra os 90 piC de referência:

| faixa de `score` | loci | dentro de piC conhecido |
|---|---|---|
| 0–2 | 799 | 18,0% |
| 2–10 | 1.811 | 33,5% |
| 10–50 | 6.633 | 37,3% |
| **50+** | 12.773 | **24,9%** |
| fundo | — | 28,3% |

Os loci de **maior** score são indistinguíveis do acaso. Não é monotônico, logo o score não
serve para ranquear nem para filtrar.

**Isto não se conserta ajustando nulos.** Acabamos de corrigir o nulo mais consequente (base
5' passou de 0,25 fixo para a composição do DNA do locus) e a correlação com o efeito
**caiu** de 0,398 para 0,323. Enquanto a métrica for significância, o problema permanece.

## O que muda

`score` passa a ser o mínimo das **forças** das evidências definidoras, em `[0, 1]`. As
colunas `s_*` continuam como estão, medindo significância — que é uma pergunta legítima e
diferente ("isto seria surpreendente sob o nulo?").

Leitura da nova escala: **`score = 0,61` significa que a evidência mais fraca que define a
classe está a 61% do caminho entre o nulo e o máximo possível, com 95% de confiança.**

## A fórmula

Toda evidência nossa é uma proporção observada contra um nulo, então uma fórmula serve para
todas:

```
L      = limite inferior de Wilson da proporção k/n ao nível conf
força  = max(0, (L − p₀) / (1 − p₀))
```

Wilson e não a aproximação normal porque nossos loci têm proporções extremas (`pct_1u` de
0,75) e amostras pequenas (`distinct_seqs` mediano 365, mas cauda longa até dezenas), onde a
normal produz intervalos fora de `[0,1]`.

A normalização por `(1 − p₀)` põe todas as evidências na mesma escala, independente de o
nulo ser 0,10 ou 0,33: mede-se **a fração do espaço disponível acima do nulo** que o efeito
ocupa, não o excesso absoluto.

Propriedades que a fórmula garante, e que o p-valor não dá:
- **não cresce com `n`** — cresce com a *precisão*, e satura no efeito observado;
- **mesmo efeito, mesma força** (a menos da largura do intervalo), independente do tamanho
  do locus;
- **comparável entre evidências** — 0,6 de 1U e 0,6 de ping-pong significam a mesma coisa;
- **zero quando o efeito não supera o nulo**, inclusive quando o observado está *abaixo*
  dele.

## As cinco evidências

Só entram as que aparecem em `defining_evidence` nos dois configs (conferido: `s_fractop`
não é usada por classe nenhuma, então fica fora):

| evidência | p̂ | n | p₀ |
|---|---|---|---|
| `1u` | `pct_1u` | `distinct_seqs` | `base_null(seq, "T", frac_top)` |
| `5a` | `pct_5a` | `distinct_seqs` | `base_null(seq, "A", frac_top)` |
| `fractop_bid` | `min(frac_top, 1−frac_top)` | `distinct_seqs` | 0,10 |
| `pingpong` | `counts[10] / Σcounts` | `Σcounts` | `taxa_de_fundo / Σcounts` |
| `phase` | `observados / sorteios` | `sorteios` | `sucessos / população` |

`n = distinct_seqs` em três delas pela mesma razão já adotada no resto da tool: reads
idênticas empilhadas e duplicatas de PCR não são ensaios independentes.

## Arquitetura

Três peças, cada uma testável sozinha:

1. **`evidence.effect_strength(k, n, p0, conf=0.95)`** — pura, sem dependência do resto.
   Devolve `[0,1]`. Guardas: `n ≤ 0`, `k < 0`, `k > n`, `p0 ≥ 1` → `0.0`.
2. **`features.build_record`** — calcula `e_1u`, `e_5a`, `e_fractop_bid`, `e_pingpong`,
   `e_phase` ao lado das `s_*` correspondentes, reusando os mesmos `k`, `n`, `p0` que já
   estão ali. Não são colunas do TSV: ficam no `ClusterRecord` para o decider consumir.
3. **`decider.defining_scores`** — passa a ler `e_{nome}` em vez de `s_{nome}`. A regra de
   combinação (mínimo no topo, máximo dentro de lista) **não muda**, nem a sintaxe do
   `rules.yaml`.

O `reporter` não ganha colunas: `score` já existe e só muda de significado. As `s_*` seguem
como estão.

## Consumidores verificados (o que NÃO muda)

Conferido no código antes de escrever isto:

| | muda? | por quê |
|---|---|---|
| classe atribuída | **não** | `_score_class` soma pesos sobre features cruas (`pct_1u`, `frac_top`, `pingpong_frac`, `phase_score`); nunca lê `s_*` nem `score` |
| flag `review` | **não** | `decider.py:123-128` usa `confidence`, `margin`, `reads`, `dicer_call` — o `score` não entra |
| recall de piC / TAS / miRBase | **não** | dependem de classe e review |
| `extras/simulate_rules.py` | **não** | reconstrói o registro com 7 features cruas; `defining_scores` já usa `getattr(..., 0.0)` e o comportamento é idêntico com `e_*` |
| coluna `score` do TSV e `score=` do GFF3 | **sim** | é o objetivo |

Os únicos consumidores do `score` são o `reporter` (coluna do TSV, atributo do GFF3) e os
testes. Nada decide com ele.

## Endurecer `defining_scores` (entra no mesmo trabalho)

`defining_scores` faz `getattr(record, f"s_{name}", 0.0)` — **atributo ausente vira 0,0 em
silêncio**. Um erro de digitação em `defining_evidence` produz score zero sem avisar
ninguém, e a classe inteira passa a parecer sem evidência.

É o mesmo padrão que já custou caro aqui: o `extras/simulate_rules.py` **aborta** quando
falta coluna, com o comentário "defaultar em silêncio faz a simulação discordar da corrida
por motivos que não têm nada a ver com a mudança de regra" — porque isso de fato aconteceu
(`pct_5a` faltando inventou 3.080 mudanças de classe). O `decider` não tem essa guarda.

**Mudança:** nome de evidência desconhecido vira `KeyError` com a lista dos nomes válidos,
na construção do `RuleScorer` (falha no carregamento do config, não a cada locus).

**Exceção deliberada:** o registro reduzido do `simulate_rules.py` não tem as `e_*`, e isso
é legítimo. A validação é sobre o **nome estar no vocabulário**, não sobre o atributo existir
no objeto — então o simulador segue funcionando e um `defining_evidence: [fractopp]` passa a
estourar.

## Critério de aceite

**Obrigatório — ordenação monotônica.** A assinatura de piRNA (`pct_1u`) tem que **subir**
a cada faixa de `score`.

**CRITÉRIO CORRIGIDO (2026-08-08).** A primeira versão deste critério media a fração de loci
dentro dos 90 piC de referência, e estava **errada**. Ela pergunta *"os loci de score alto
estão nos clusters conservados?"* quando a pergunta é *"os loci de score alto têm assinatura
de piRNA?"*. Os piC de referência cobrem 5,66% do genoma e só os clusters conservados entre
cepas — exigir concentração ali **penaliza a tool por encontrar piRNA legítimo em outro
lugar**, que é justamente o trabalho dela.

Medido: os **4.542** loci de score ≥ 0,6 **fora** dos piC são indistinguíveis dos 1.689 de
dentro — 1U de 0,82 nos dois, tamanho 24–31 nt em ~85% nos dois, ping-pong 0,068 × 0,079. O
controle de score < 0,2 fora tem 1U de 0,38 e ping-pong **zero**. Ou seja, o score separa
piRNA de não-piRNA dentro e fora da referência; a fração dentro do piC caía nas faixas altas
apenas por diluição, porque há milhares de loci genuínos espalhados pelo genoma.

**Obrigatório — piso intacto.** piC 85/90, TAS 6/6, miRBase 134/134 e 244/244. Como o
classificador não usa `s_*` nem `e_*` (usa `pct_1u` cru, `decider.py:63`), **nenhum locus
pode trocar de classe**. Se trocar, há bug.

**Esperado — correlação invertida.** A correlação do `score` com o efeito tem que passar a
ser maior que com `distinct_seqs`. Hoje é 0,323 contra 0,629.

## Testes

Da função pura:
- efeito igual em `n` diferente dá força parecida (é a propriedade que hoje falha);
- `n` maior com o mesmo `p̂` dá força **maior ou igual**, nunca menor (o intervalo aperta);
- observado igual ao nulo dá 0; observado abaixo do nulo dá 0;
- observado em 1,0 com `n` grande dá perto de 1,0;
- guardas de entrada degenerada.

De integração:
- um locus com efeito forte e `n` pequeno **não** perde para um com efeito fraco e `n`
  grande — inversão que o `score` atual produz;
- `score` continua sendo o mínimo das definidoras e respeitando slots OU.

## Limitações que ficam, e vão documentadas

- **Ping-pong:** `n` são *pares*, que escalam com o quadrado da cobertura e não são ensaios
  independentes. A força será otimista ali — o p-valor já é, pelo mesmo motivo.
- **Fase:** a janela é escolhida por maximizar evidência, então o efeito estimado é enviesado
  para cima. O limite inferior atenua, não elimina. (O `s_phase` trata isso com Bonferroni,
  que é correção de *significância* e não se aplica a tamanho de efeito.)
- **`e_fractop_bid` herda o ruído do `s_fractop_bid`:** a fração minoritária separa uni-fita
  de bifita por apenas 1,5–1,7× no nível do locus. Trocar a métrica não conserta uma
  quantidade que é ruidosa na origem.
- **A escala muda.** Um `score` de 50 vira algo como 0,7. Corridas antigas não comparam com
  novas — mas já não comparavam entre si, porque a escala antiga não media o que dizia.

## Fora de escopo

- Rever `defining_evidence` de cada classe (quais evidências definem o quê).
- A categoria "não classificado" por score baixo — depende desta mudança, vem depois.
- `s_fractop`, `s_fractop_bal`, `overlap_dup_frac`, `size_purity`, `dicer_call_2`,
  `unique_frac`: seguem reportadas e fora de regra.
