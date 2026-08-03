# Design: score de evidência estatística por locus

**Autor:** Lucas
**Data:** 2026-08-02
**Contexto:** revisão de literatura em `2026-08-02-classification-confidence-review.md`

## Problema

O `confidence` atual é um softmax sobre a soma de pesos das classes. Ele é
**relativo**: compara as classes entre si e sempre elege um vencedor. Um locus
que não se parece com nada recebe confiança alta se as outras classes se
encaixarem ainda menos. Consequências medidas na rodada de Arabidopsis:

- 184 miRNAs anotados receberam classe errada **afirmada com confiança**;
- 495 loci de organela (11,7% da massa de reads) foram chamados `hcsiRNA`,
  numa região onde a via RdDM nem existe;
- 475 loci sobre tRNA/snoRNA receberam confiança 0,584, contra 0,597 da média —
  o score não os distingue em nada.

O objetivo é que o score seja **filtrável** e signifique certeza: acima de um
corte escolhido pelo usuário, a chamada é sustentada por evidência que não
seria esperada por acaso.

## Abordagem

Seguir o proTRAC (Rosenkranz & Zischler 2012): cada assinatura é testada contra
um nulo explícito e vira `score = log₁₀(1/p)`. Escala sempre positiva e
comparável entre features — score 2 é p = 0,01, score 3 é p = 0,001.

### Nulos por feature

| feature | teste | nulo |
|---|---|---|
| `pct_1u`, `pct_5a` | binomial unilateral | **EM ABERTO — ver "Decisão pendente" abaixo** |
| `frac_top` | binomial bilateral | 0,5 (fitas equiprováveis) |
| `pingpong_frac` | binomial unilateral | 1/20 (um par cair no offset 10 entre os 20 offsets testados) |
| `phase_score` | hipergeométrico | **já implementado** (`signals.phase_test`) |
| densidade de reads | Poisson | `lib_size` ponderado / bases das referências mapeadas × span do locus |

O `n` dos testes binomiais é a contagem ponderada por 1/NH arredondada, mantendo
a coerência com o resto do pipeline.

### Decisão pendente: o nulo do viés de base 5'

Nenhuma opção resolve os dois problemas ao mesmo tempo, e a escolha é
metodológica — não deve ser fixada sem fundamento.

| nulo candidato | cancela viés de ligase? | sofre auto-referência? |
|---|---|---|
| 0,25 fixo | não | não |
| frequência média da biblioteca | **sim** (o viés é constante dentro dela) | **sim** (a classe dominante define o fundo) |
| posições 3–8 da própria read | não (não passam por ligação) | não |

**O viés de ligase é real e documentado:** a T4 RNA ligase 1 liga preferindo
C > U ≥ A > G na extremidade 5' (Jayaprakash et al. 2011), e boa parte do viés
vem da estrutura secundária dos substratos (Sorefan et al. 2012). Logo, a base
5' observada mistura seleção por Argonauta com artefato de protocolo.

**Auto-referência medida nestes dados:** a frequência de 5'A na biblioteca de
Arabidopsis é ~0,49, e a média de `pct_5a` dos loci `hcsiRNA` é 0,51. Usar a
média da biblioteca como nulo tornaria a classe dominante indetectável.

Resolver antes de implementar: verificar que nulo o proTRAC usa para 1T/10A e
se alguma ferramenta trata o viés de ligase explicitamente. Enquanto não houver
fundamento, `0,25` é o default conservador — erra por excesso de significância,
mas não se auto-sabota nem depende de premissa não verificada.

### Score da classe

Cada classe declara no `rules.yaml` quais evidências a **definem**. O score da
classe é o **mínimo** entre os `log₁₀(1/p)` dessas evidências, com p corrigido
por Bonferroni pelo número de testes realizados.

Mínimo, não soma: somar deixaria uma evidência forte compensar uma ausente — um
locus com densidade altíssima e nenhum viés de base ainda somaria score alto,
que é exatamente o falso positivo que queremos eliminar. Com o mínimo, o score
é o elo mais fraco, e a leitura fica direta: **score 2 significa que todas as
evidências definidoras da classe têm p ≤ 0,01**. É a convergência de critérios
que o proTRAC exige, na forma mais simples.

Evidências não definidoras (ex.: ping-pong numa classe de planta) continuam
sendo reportadas em suas colunas, mas não entram no score.

### Saída

O `.clusters.tsv` ganha `score` (da classe) e uma coluna por evidência
(`s_1u`, `s_5a`, `s_fractop`, `s_pingpong`, `s_phase`, `s_density`).
`cls`, `confidence` e `support` permanecem, para não quebrar consumidores
existentes (`classplots`, `class_summary`, validadores).

**Não haverá coluna booleana de confirmação.** Um `confirmed=yes/no` embutiria
um limiar que é decisão do usuário; o score contínuo deixa o corte com quem
analisa.

## O que NÃO entra, e por quê

Cada item abaixo foi testado nesta sessão e descartado por evidência, não por
preferência:

- **Filtro de organelas.** Exigir 5'A significativo já rejeita **71%** dos loci
  de organela contra 40,6% dos nucleares, sem que a ferramenta precise saber o
  que é organela — o que preserva a premissa de funcionar sem anotação. Uma
  flag `--exclude-refs` pode ser oferecida depois, se o usuário quiser.
- **Heterogeneidade de tamanho ao longo do locus.** Separava perfeitamente em
  Arabidopsis (nuclear 0,00 vs organela 0,50) e **inverteu em Drosophila**
  (nuclear 0,64 vs rDNA 0,57), porque piRNA tem distribuição de tamanho larga.
  Media biologia de reino, não qualidade de locus.
- **Múltiplas classes por locus.** Loci nucleares são homogêneos ao longo de si
  mesmos (0–14% de variação); a mistura aparente vinha de organela. Só 39 loci
  (0,1%) passam de 5 kb.
- **natsiRNA.** Define-se pela origem (transcritos antisense sobrepostos), não
  pelo perfil de reads; exigiria anotação gênica.
- **ML supervisionado sobre os rótulos atuais.** Seria circular: aprenderia as
  regras que queremos corrigir.

## Limitação conhecida

O teste tem menos poder em loci rasos. Dos `hcsiRNA` nucleares, os 18.632 que
seriam rejeitados ainda mostram **2,71× de enriquecimento em transposons**
(contra 3,24× dos que passam), com 19 reads medianas contra 60. São
provavelmente reais, mas sem profundidade para prová-lo.

**Score baixo significa "não comprovado", não "falso".** Isso precisa estar no
README, para que ninguém filtre acreditando estar removendo lixo.

## Validação

A regra vale para os dois reinos ou não entra. Antes e depois, sobre as tabelas
de `docs/validation/`:

1. Recall dos 6 loci TAS (Arabidopsis) e dos 90 clusters de piRNA (Drosophila)
   não pode cair.
2. Enriquecimento em TE do `hcsiRNA` e nos piCs do `piRNA` deve subir ou ficar
   estável ao filtrar por score.
3. Os 495 loci de organela devem concentrar-se na faixa de score baixo.
4. Os 41 miRNAs do ShortStack (100% de precisão) servem de positivo de ouro.

## Fora de escopo (fica para depois)

Subclasses baratas — phasiRNA 21 vs 24nt e piRNA primário vs secundário — usam
features já presentes no TSV e entram num spec próprio, depois que o score
estiver de pé. Clustering não-supervisionado como diagnóstico das classes
também: ele fica melhor com as features já convertidas em score.

## Referências

Jayaprakash, A. D., Jabado, O., Brown, B. D., & Sachidanandam, R. (2011). Identification and remediation of biases in the activity of RNA ligases in small-RNA deep sequencing. *Nucleic Acids Research*, 39(21), e141. https://doi.org/10.1093/nar/gkr693

Rosenkranz, D., & Zischler, H. (2012). proTRAC — a software for probabilistic piRNA cluster detection, visualization and analysis. *BMC Bioinformatics*, 13, 5. https://doi.org/10.1186/1471-2105-13-5

Sorefan, K., Pais, H., Hall, A. E., Kozomara, A., Griffiths-Jones, S., Moulton, V., & Dalmay, T. (2012). Reducing ligation bias of small RNAs in libraries for next generation sequencing. *Silence*, 3(1), 4. https://doi.org/10.1186/1758-907X-3-4

Johnson, N. R., et al. (2025). Class-agnostic annotation of small RNAs balances sensitivity and specificity in diverse organisms. *Computational and Structural Biotechnology Journal*. https://doi.org/10.1016/j.csbj.2025.05.045
