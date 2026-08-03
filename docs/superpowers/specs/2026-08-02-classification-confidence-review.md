# Revisão de literatura: confiança na classificação de populações de small RNAs

**Autor:** Lucas
**Data:** 2026-08-02
**Objetivo:** embasar o redesenho do score do `classify` para que ele seja
filtrável e signifique certeza real, e avaliar a viabilidade de subclasses
(hc-siRNA vs natsiRNA, piRNA primário vs secundário, phasiRNA 21 vs 24nt).

Busca: Europe PMC (10 consultas temáticas, ordenadas por citação), leitura dos
textos completos dos quatro métodos centrais, DOIs verificados no CrossRef.

---

## 1. O problema é reconhecido pelo campo

A revisão mais recente do assunto (Johnson et al. 2025) afirma o que já
suspeitávamos, com todas as letras:

- **Viés de classe**: a maioria das ferramentas trata só miRNA, ignorando
  siRNA e piRNA.
- **Over-annotation**: anotam regiões de expressão trivialmente baixa como
  loci — é exatamente o nosso caso, com 70.239 loci e 33% em review.
- **Escopo de organismo**: "virtually no tools are designed for de novo sRNA
  annotation outside of plants and animals".
- **Não existe gold standard**: eles próprios usam reprodutibilidade entre
  réplicas (Jaccard) por falta de verdade de referência.

Ou seja, o gargalo que motivou esta revisão é real e está documentado.

## 2. Como cada ferramenta decide — e o que isso nos ensina

| ferramenta | o que afirma | mecanismo | limite declarado |
|---|---|---|---|
| **ShortStack** (Axtell 2013) | `DicerCall` + `MIRNA` Y/N | critério estrutural (hairpin, miRNA*) | recusa-se a nomear outras classes |
| **unitas** (Gebert et al. 2017) | miRNA, tRF, rRF, piRNA | 100% banco de dados, hierárquico | **não faz de novo**; piRNA é a classe residual |
| **proTRAC** (Rosenkranz et al. 2012) | cluster de piRNA | **Poisson + binomial, score = log₁₀(1/p)** | específico de piRNA cluster |
| **YASMA** (Johnson et al. 2025) | loci, sem classe | knee-finding, filtros de densidade | class-agnostic por opção |

**A lição central vem do proTRAC.** O score dele *é* um p-valor em escala log:
`score = log₁₀(1/p)`, com `score ≥ 2` equivalendo a `p ≤ 0,01`. Cada
característica (densidade de hits, assimetria de fita, viés 1T/10A,
distribuição de tamanho) é testada contra um nulo explícito — Poisson para
densidade sobre distribuição uniforme, binomial 50/50 para fita — e o
conjunto leva correção de Bonferroni. Exige **convergência de todos os
critérios**, não só densidade alta, justamente para que acúmulos de
transposon multi-cópia não virem cluster funcional.

Nosso score é uma soma de pesos arbitrários; o deles é interpretável. Essa é
a diferença entre reportar "0,73" e reportar "p ≤ 0,01".

**O unitas mostra o que não fazer**: hierarquia dependente de anotação, com
piRNA recebendo o que sobrou. Para organismo não-modelo — nosso alvo — isso
não funciona, e atribuir uma classe por eliminação é a mesma falha que
produziu nossos 184 miRNAs mal rotulados.

**O YASMA confirma a direção dos filtros**: mínimo de 50 reads, 100 reads/kb e
10 reads únicas/kb, explicitamente "favorecendo reduzir falsa descoberta ao
custo de alguma sensibilidade".

## 3. Viabilidade das subclasses pedidas

| distinção | assinatura discriminante | temos as features? |
|---|---|---|
| **phasiRNA 21 vs 24nt** | tamanho + phasing; os de 24nt são reprodutivos, dependentes de DCL5, presentes em angiospermas (Xia et al. 2019) | **sim, já medidas** |
| **piRNA primário vs secundário** | secundário = produto do ciclo ping-pong (10A, sobreposição 10nt); primário = fásico, dependente de Zucchini, 1U (Han et al. 2015; Mohn et al. 2015) | **sim, já medidas** |
| **hc-siRNA vs natsiRNA** | natsiRNA se define pela **origem** (pares de transcritos antisense sobrepostos), não pelo perfil de reads | **não** — exige anotação gênica |

As duas primeiras são subdivisões das classes atuais usando features que já
estão no `.clusters.tsv`. A terceira depende de anotação e é frágil
exatamente onde a ferramenta quer ser forte.

## 4. Métodos de fora da área

- **IDR** (Irreproducible Discovery Rate), padrão ENCODE em ChIP-seq/ATAC-seq:
  usa réplicas para achar onde o sinal vira ruído e produz medida
  interpretável. Não encontrei aplicação a loci de sRNA. Limite prático: exige
  réplicas, que raramente existem em sRNA-seq — serve para **calibrar e
  validar**, não como requisito de uso.
- **Positive-unlabeled learning** (Li et al. 2021): treina com poucos
  positivos confiáveis e muitos não rotulados, em vez de fingir que se conhece
  o rótulo de todos. É o enquadramento honesto para o nosso caso — temos
  poucos positivos de ouro (6 TAS, 90 piCs, 41 miRNAs do ShortStack) contra
  70 mil loci sem rótulo.
- **Alocação fracionada de multimapeadores** (Johnson et al. 2016; e a equação
  `c = Σ rᵢ/hᵢ` do unitas): mesma família do nosso `1/NH`, já implementado.

## 5. Implicações para o srnaprofiler

1. **Trocar o score por evidência estatística por feature.** Seguir o proTRAC:
   nulo explícito por feature (binomial para 1U/5'A/fita, Poisson para
   densidade), score em `log₁₀(1/p)`, correção para múltiplos testes. Já temos
   um caso funcionando — o `phase_pvalue` hipergeométrico — e ele mostrou que
   dá para fazer.
2. **Exigir convergência, não soma.** O proTRAC não deixa densidade alta
   sozinha nomear um cluster. Nosso `min_evidence` foi um primeiro passo nessa
   direção, mas ainda é soma ponderada.
3. **Filtros de densidade de locus** (reads/kb, reads únicas/kb), como YASMA,
   antes de classificar.
4. **Subclasses baratas primeiro**: phasiRNA 21/24 e piRNA primário/secundário
   saem das features existentes; natsiRNA fica fora até haver anotação.
5. **Réplicas para validar o p-valor, não para rodar.** Se os loci com p baixo
   forem os que se reproduzem entre rep1/2/3, o modelo nulo está calibrado — e
   aí o limiar vale para organismo sem réplica.
6. **Não treinar ML supervisionado nos nossos próprios rótulos.** Seria
   circular e reproduziria os erros com confiança maior. PU-learning ou
   não-supervisionado são os enquadramentos defensáveis.

## Referências

Aparicio-Puerta, E., Hirsch, P., Schmartz, G. P., Kern, F., Fehlmann, T., & Keller, A. (2022). sRNAbench and sRNAtoolbox 2022 update: accurate miRNA and sncRNA profiling for model and non-model organisms. *Nucleic Acids Research*, 50(W1), W710–W717. https://doi.org/10.1093/nar/gkac363

Axtell, M. J. (2013). ShortStack: Comprehensive annotation and quantification of small RNA genes. *RNA*, 19(6), 740–751. https://doi.org/10.1261/rna.035279.112

Gebert, D., Hewel, C., & Rosenkranz, D. (2017). unitas: the universal tool for annotation of small RNAs. *BMC Genomics*, 18(1), 644. https://doi.org/10.1186/s12864-017-4031-9

Han, B. W., Wang, W., Li, C., Weng, Z., & Zamore, P. D. (2015). piRNA-guided transposon cleavage initiates Zucchini-dependent, phased piRNA production. *Science*, 348(6236), 817–821. https://doi.org/10.1126/science.aaa1264

Johnson, N. R., Larson, N. B., & Axtell, M. J. (2016). Improved Placement of Multi-mapping Small RNAs. *G3: Genes|Genomes|Genetics*, 6(7), 2103–2111. https://doi.org/10.1534/g3.116.030452

Johnson, N. R., et al. (2025). Class-agnostic annotation of small RNAs balances sensitivity and specificity in diverse organisms. *Computational and Structural Biotechnology Journal*. https://doi.org/10.1016/j.csbj.2025.05.045

Li, F., Dong, S., Leier, A., Han, M., Guo, X., Xu, J., Wang, X., Pan, S., Jia, C., Zhang, Y., Webb, G. I., Coin, L. J. M., Li, C., & Song, J. (2021). Positive-unlabeled learning in bioinformatics and computational biology: a brief review. *Briefings in Bioinformatics*, 23(1), bbab461. https://doi.org/10.1093/bib/bbab461

Mohn, F., Handler, D., & Brennecke, J. (2015). piRNA-guided slicing specifies transcripts for Zucchini-dependent, phased piRNA production. *Science*, 348(6236), 812–817. https://doi.org/10.1126/science.aaa1039

Rosenkranz, D., & Zischler, H. (2012). proTRAC — a software for probabilistic piRNA cluster detection, visualization and analysis. *BMC Bioinformatics*, 13, 5. https://doi.org/10.1186/1471-2105-13-5

Srivastav, S. P., Feschotte, C., & Clark, A. G. (2024). Rapid evolution of piRNA clusters in the *Drosophila melanogaster* ovary. *Genome Research*, 34(5), 711–724. https://doi.org/10.1101/gr.278062.123

Xia, R., Chen, C., Pokhrel, S., Ma, W., Huang, K., Patel, P., Wang, F., Xu, J., Liu, Z., Li, J., & Meyers, B. C. (2019). 24-nt reproductive phasiRNAs are broadly present in angiosperms. *Nature Communications*, 10(1), 627. https://doi.org/10.1038/s41467-019-08543-0
