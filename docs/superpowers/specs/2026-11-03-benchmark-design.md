# Benchmark do sRNAProfiler — design

Data: 2026-11-03. Status: rascunho para aprovação.
Parte 1 de 5 do plano de evolução da ferramenta (ver "Contexto" abaixo).

## Intenção e por quê

O sRNAProfiler chama loci de sRNA a partir de cobertura, mede as assinaturas de cada
via de biogênese e atribui uma classe por regras editáveis (`rules.yaml`,
`rules.plant.yaml`). O objetivo do trabalho é um **artigo** que mostre duas coisas:

1. a ferramenta classifica melhor do que as alternativas existentes;
2. uma camada de aprendizado de máquina melhora a classificação por regras.

Nenhuma das duas afirmações é mensurável hoje. Falta um **benchmark**: um conjunto de
dados públicos com resposta conhecida e métricas definidas. Sem ele não se sabe quanto
falso positivo existe, em qual classe, nem se qualquer melhoria posterior melhorou algo.

Motivação concreta: a rodada existente `results-teste/dros_ref.clusters.tsv`
(*Drosophila*, dados públicos cuja origem não foi registrada) tem 25.738 loci, dos quais
18.349 (71%) chamados piRNA, 4.453 siRNA, 2.801 ambíguos e 135 miRNA. Há 3.614 loci com
100 pb ou menos, 987 com 50 pb ou menos, e 18.383 loci (71%) com menos de 100 reads. A
mediana de tamanho dos loci piRNA é 551 pb, enquanto clusters de piRNA reais costumam ter
dezenas a centenas de kb. Hipótese (não verificada): clusters verdadeiros estão sendo
fragmentados e pedaços pequenos ou pouco suportados entram na contagem como loci.

## Escopo

- **Táxons:** artrópodes e plantas.
- **Dados:** somente públicos. As bibliotecas de *Apis mellifera* do grupo
  (`~/Libs_smallRNA_LSV_eric` no litrp3) **não** entram.
- **Entrada da ferramenta no benchmark:** somente genoma e reads. Nenhuma anotação de genes,
  TE ou outra camada opcional. As camadas opcionais serão avaliadas depois, como ganho
  medido à parte (parte 3 do plano).

## O que é falso positivo aqui

Duas famílias, medidas separadamente:

1. **Classe errada:** o locus existe, mas recebeu a classe errada (ex.: siRNA chamado de
   piRNA), ou foi chamado numa região que o gabarito diz não ter essa classe.
2. **Locus mal delimitado:** o locus é minúsculo, fragmentado, sustentado por poucas
   reads, ou junta regiões biologicamente distintas.

## Definição de cluster (decisão do usuário)

Como no proTRAC: loci vizinhos da mesma natureza são **juntados** em regiões maiores, mas
podem ser **separados em loci distintos quando a biologia pedir**. Candidatos a critério de
separação (palpite, a definir na parte de delimitação): mudança de classe, de fita
(uni versus bifita), de registro de fase, ou queda de cobertura entre dois domínios. O
benchmark deve premiar as duas direções: recuperar um cluster real inteiro e não colar loci
que são distintos. A *implementação* da delimitação não faz parte desta parte; o benchmark
só a mede.

## Revisão de 2026-11-03 (decisões do usuário, prevalecem sobre o texto abaixo)

- **Mutantes saem do benchmark.** O gabarito passa a ser só anotações curadas e
  simulação; a métrica de falso positivo por mutante e os limiares de fold-change caem.
  Consequência: em dado real, falso positivo só é medido onde a anotação é ampla; locus
  fora da anotação continua "desconhecido". A simulação passa a carregar mais peso.
- **Download de SRA com iSeq.**
- **O módulo `mirna` roda em todas as bibliotecas**, mesmo sem referência de miRBase
  (`--mature` ausente), com `mirdeep2` para artrópodes e `mirdeep-p2` para plantas.

## Fontes de gabarito (do mais forte ao mais fraco)

1. **Mutantes como verdade biológica.** Um locus chamado numa classe que persiste num
   mutante em que aquela via está desligada é falso positivo comprovado, sem depender de
   anotação de terceiros. Exemplos de pares mutante/selvagem a procurar (a confirmar na
   literatura): *Drosophila*, mutantes das vias piRNA (*piwi*, *aub*, *ago3*) e siRNA
   (*dcr-2*); *Arabidopsis*, mutantes *dcl* e *rdr* para siRNA de 24 nt e ta-siRNA.
2. **Anotações curadas** como gabarito positivo: clusters de piRNA (ex.: piRBase,
   *flamenco*, 42AB), miRNAs (miRBase, MirGeneDB), loci PHAS e TAS da literatura.
3. **Simulação:** loci sintéticos de classe conhecida com ruído realista, para controlar a
   dificuldade e medir falso positivo sem depender das imperfeições das anotações.

## Organismos (palpite do assistente, a confirmar)

- *Drosophila melanogaster*: base. Melhor gabarito e mutantes bem estudados; já existe uma
  rodada de comparação.
- Mais um artrópode, para não ficar só em mosca: candidatos *Aedes aegypti* ou *Bombyx mori*.
- *Arabidopsis thaliana*: mutantes *dcl* e *rdr*, loci TAS.
- Arroz ou milho: phasiRNA reprodutivos.

A escolha final de táxons e dos acessos é do construtor, desde que cubra artrópode e
planta e inclua *Drosophila* e *Arabidopsis*.

## Métricas

- **Classificação:** precisão, recall e FDR por classe; calibração da confiança.
- **Delimitação:** fração de cada cluster verdadeiro recuperada; número de fragmentos por
  cluster verdadeiro; número de loci minúsculos ou com poucas reads que sobram; número de
  loci verdadeiros distintos que foram colados.
- **Validação cruzada por espécie:** quando houver ML (parte 4), treinar sem uma espécie e
  testar nela, para mostrar generalização.

## Dados e reprodutibilidade

- Lista de acessos (SRA e genomas) **versionada no git**, com metadados: organismo, tecido,
  genótipo, replicata, fonte do gabarito.
- Script que baixa tudo para uma pasta `data/` **fora do git** (adicionar ao `.gitignore`),
  no litrp3.
- Cada rodada registra a versão da ferramenta, os parâmetros e os acessos usados, para que o
  resultado seja reproduzível. A origem dos dados do `dros_ref` não foi registrada; essa
  rodada serve só como ponto de partida, não como resultado.

## Metas, não-metas, anti-metas

**Metas**
- Um conjunto de benchmark reproduzível e as métricas acima, rodando de ponta a ponta.
- Um primeiro número de falso positivo por classe e de qualidade de delimitação para a
  ferramenta atual.

**Não-metas (nesta parte)**
- Melhorar o classificador, a delimitação ou adicionar ML. Isso vem depois.
- Comparar com outras ferramentas (parte 5), embora o formato do benchmark deva permitir.
- Camadas opcionais de evidência (parte 3).

**Anti-metas** (o que seria fracasso mesmo que "funcione")
- Um gabarito que a ferramenta já tenha visto ou do qual dependa (circularidade).
- Resultados que só valem para *Drosophila*.
- Métricas que escondam o problema de loci minúsculos e fragmentados.

## Restrições

- Ambiente de execução: litrp3 (Debian 13, 56 cores, 1 TB RAM, ~2 TB livres em `/home`).
  Env conda `srnaprofiler` criado a partir de `environment.yml`, em
  `/media/nas1/LITRP.DBs/metagenDBs/conda_envs/srnaprofiler`.
- Sem sudo no litrp3.
- Python 3.13 é fixado pelo `environment.yml` (motivo documentado ali).

## Estado do repositório que o construtor precisa saber

O commit `ca6ddc2` ("chore: only what the pipeline needs to run") **removeu `tests/` e
`extras/` do repo**. Consequências:

- O README ainda cita 398 testes, que não estão mais no repo.
- `extras/fetch_ncrna_db.sh`, necessário para `setup --data`, não existe mais; o comando
  falha como está.
- Esse commit removeu também material útil como semente do benchmark, recuperável de
  `ca6ddc2^`: `extras/piC_ref_dm6.bed` (clusters de piRNA de referência, dm6),
  `extras/phas21_ref_at_TAIR10.gff3` e `extras/phas_experimental_at_TAIR10.bed`
  (loci PHAS de *Arabidopsis*), e os scripts `extras/validate_*.py`,
  `extras/compare_runs.py` e `extras/simulate_rules.py`.

Decisão pendente do usuário (fora do benchmark): se esses arquivos voltam ao repo, ficam
numa pasta separada, ou ficam só no histórico.

## O que fica para o construtor

Escolha dos acessos e táxons dentro do escopo; formato interno do gabarito; estrutura de
pastas do benchmark; como medir "cluster recuperado" (limiares de sobreposição); como gerar
a simulação.

## Decisões do assistente, para o usuário conferir

- Documento escrito em português (o README do repo está em inglês).
- Os táxons e exemplos de mutantes acima são sugestões, não verificadas na literatura.
- O benchmark é a primeira parte; a delimitação de loci é uma parte separada, logo depois.
- **Gabarito de mutantes:** avaliado por **região fixa**, com fold-change de cobertura
  (RPM) entre selvagem e mutante, e não por presença ou ausência estrita, porque mutantes
  costumam deixar piRNA ou siRNA residual. Reportar a taxa de falso positivo em vários
  limiares de fold-change em vez de fixar um.
- **Quais métricas vêm de qual gabarito:** classificação usa mutantes e anotações curadas;
  **delimitação só usa o subconjunto curado** (clusters de piRNA, loci PHAS e TAS) e a
  simulação. Mutantes não dão fronteiras. Quando duas anotações curadas discordam (ex.:
  tamanhos ou fusões diferentes), reportar contra cada fonte separadamente, sem fundi-las.
- **Classe "ambiguous" e loci sem anotação:** "ambiguous" conta como abstenção, com uma
  métrica de cobertura ao lado da precisão. Locus em região sem anotação conta como
  "desconhecido", fica fora da precisão, mas é sempre reportado à parte com sua distribuição
  de tamanho e de reads, para não esconder o problema dos loci minúsculos. Nos mutantes,
  esses loci são testáveis.
- **Loci minúsculos ou com poucas reads:** só reportados, em vários limiares (≤50 pb,
  ≤100 pb; <20, <50, <100 reads). Os limiares ficam num arquivo de configuração versionado.
  Sem meta de aprovação por enquanto.
- **Simulação:** loci sintéticos injetados em bibliotecas reais como fundo, com parâmetros
  gerados a partir de distribuições da literatura e não das regras da ferramenta, para
  reduzir circularidade. O risco de circularidade residual deve ser dito no artigo.
- **Orçamento de dados:** começar pequeno (por organismo, selvagem mais um ou dois mutantes,
  2 a 3 réplicas), *Drosophila* primeiro, e crescer depois de fechar o ciclo completo.
- **Arquivos recuperados de `ca6ddc2^`:** podem ser copiados para uma pasta própria do
  benchmark, com a proveniência anotada.

## Circularidade: espécies de ajuste e espécies de teste (resposta do usuário)

Segundo o usuário, as regras e limiares de `rules.yaml` e `rules.plant.yaml` foram ajustados
com *Drosophila*, *Arabidopsis* e **arroz**. Consequências para o desenho:

- *Drosophila*, *Arabidopsis* e arroz são **espécies de ajuste**. Números obtidos nelas
  medem o quanto a ferramenta se ajusta ao que já viu, não o quanto generaliza. Podem ser
  reportados, rotulados como "ajuste", mas não sustentam a afirmação de que a ferramenta
  classifica melhor.
- O **teste** precisa de espécies que não moldaram as regras: pelo menos um artrópode além
  de *Drosophila* (ex.: *Aedes aegypti* ou *Bombyx mori*) e pelo menos uma planta além de
  *Arabidopsis* e arroz (ex.: milho, para phasiRNA).
- Isso promove o "segundo artrópode" e o "milho" da lista de organismos de opcionais a
  **obrigatórios**. O arroz sai da lista de candidatos a teste.
- Dentro das espécies de ajuste, os arquivos `piC_ref_dm6.bed` e `phas21_ref_at_TAIR10.gff3`
  são material de ajuste e não podem ser usados como gabarito de teste.
- Validação cruzada por espécie na parte 4 (ML): treinar nas espécies de ajuste e testar nas
  espécies não vistas. Qualquer reajuste de regras após ver o teste invalida o teste.
- Mutantes, anotações curadas e simulação continuam como fontes de gabarito, mas a
  simulação deve ser gerada de forma independente das regras (ver decisões acima).
