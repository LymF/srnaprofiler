# Design: definição de loci (clustering)

**Autor:** Lucas
**Data:** 2026-08-04
**Antecede:** `2026-08-02-evidence-scores-design.md` — a branch `evidence-scores` está
parada até este trabalho entrar, porque o score mede propriedades do locus e o
locus está mal definido.

## Problema

O `classify` define loci como ilhas de cobertura fundidas quando o vão é ≤ `--pad`
(default 75). Isso nunca foi calibrado — o valor atravessou desde a substituição do
tiling de 2 kb. Medido:

- **Fragmentação.** No mesmo BAM, o ShortStack produz loci de 529 bp medianos; nós,
  88–140 bp. Metade dos nossos vãos entre loci vizinhos tem 75–200 bp, ou seja,
  cortamos por buracos pequenos demais para separar biologia.
- **Poder estatístico.** 34 reads por locus contra 444 do ShortStack. Os testes de
  evidência não têm `n` para provar nada.
- **Assinatura destruída.** A classe `siRNA` é definida por ser bifita, mas a
  distância mediana de `frac_top` a 0,5 é 0,262 — um cluster bifita picado em
  fragmentos vira uma coleção de pedaços uni-fita, e o teste de equivalência
  (corretamente) não encontra evidência.

## Mudanças

### 1. `mincov` em RPM

O `--mincov` passa a ser lido em **reads por milhão**, convertido para profundidade
absoluta uma vez por corrida, como o ShortStack faz:

```
min_depth = round(mincov_rpm * lib_size / 1e6)
```

`lib_size` é o ponderado por 1/NH que já calculamos. Hoje o valor é contagem crua e
por isso **não transfere entre bibliotecas** — o `mincov: 10` calibrado em
DRR351382 (30,8M alinhamentos) não significa nada em outra profundidade. Em RPM,
transfere.

Consequência: **invalida a calibração do item 1 do roadmap.** Recalibrar em RPM.

### 2. Extensão bilateral

Trocar "funde ilhas cujo vão ≤ pad" por "estende cada ilha em `pad` para os dois
lados, funde o que passar a se sobrepor" (o `bedtools slop` + `merge` do
ShortStack). O vão efetivo passa a ser `2 × pad`, e o pico ganha uma zona de
influência simétrica em vez de o critério olhar só o buraco entre dois picos.

### 3. `pad` por config

Os dois reinos discordam do ótimo, e isso é biologia:

| pad | Drosophila | Arabidopsis |
|---|---|---|
| 75 | recall 85/90, 28 reads/locus | TAS 6/6, 34 reads/locus |
| 200 | **recall 85/90**, 36 reads | **TAS 6/6**, 59 reads |
| 500 | recall 83/90, 61 reads | **TAS 6/6**, 131 reads |
| 1000 | recall 81/90 | — |

Em Drosophila, 500 já custa 2 clusters de piRNA; em Arabidopsis ainda é ganho puro.
Os piCs de Drosophila são regiões grandes com subestrutura que a fusão agressiva
apaga; os loci de planta são compactos. Portanto `pad` vai para o config, junto das
classes: **200** em `rules.yaml`, **500** em `rules.plant.yaml` (valores a
reconfirmar após as mudanças 1 e 2, que alteram as fronteiras).

### 4. Fusão condicionada ao perfil (YASMA)

Distância deixa de ser critério suficiente. Duas regiões só se fundem se também
compartilharem:

- **perfil de tamanho** — mesma faixa dominante (não fundir um vizinho de 24 nt com
  um de 21 nt, que são DCL3 e DCL4);
- **proporção de fita** semelhante (diferença < 0,5, como o YASMA).

Medição prévia: aumentar o `pad` de 75 para 500 **não** aumentou a mistura de
tamanhos (pureza do tamanho dominante fica em ~0,64 e a fração de loci bimodais cai
de 19,3% para 17,7%). Ou seja, este critério não está corrigindo um dano que a
fusão causa — ele previne um dano que ela **poderia** causar em outro organismo, e
torna a fusão defensável em vez de puramente geométrica.

### 5. miRNA dentro de cluster

Com loci maiores, um miRNA passa a caber dentro de um cluster que contém muito mais
coisa. Hoje `classify.mark_mirna` marca o locus inteiro como `miRNA` a qualquer
sobreposição — com `pad=75` isso afeta 1 locus acima de 2 kb, com `pad=200` já são 6,
e a curva sobe com o `pad`.

Regra nova: o locus só é rotulado `miRNA` se o miRNA **dominar** o locus (sobreposição
cobrindo a maior parte do span, limiar a definir na implementação). Caso contrário o
locus mantém sua classe e o miRNA é reportado como anotação sobreposta, sem sequestrar
o rótulo. Alternativa considerada e descartada: refinar as fronteiras do locus para o
miRNA (o que o ShortStack faz), porque descarta o restante do cluster.

### 6. Heterogeneidade de tamanho por locus

**18% dos loci têm duas populações de tamanho** (pureza do tamanho dominante < 0,5),
e isso já ocorre com `pad=75` — não é efeito de fusão. Hoje reportamos só o
`DicerCall` dominante e a segunda população desaparece do relatório.

Passa a ser reportada: `size_purity` (fração das reads no tamanho dominante) e
`dicer_call_2` (segundo tamanho mais abundante, vazio quando a pureza é alta). Sem
mudar a classe — o locus continua com um rótulo, o da população dominante.

Isso é deliberadamente conservador. Um locus com 62% de 24 nt e 31% de 21 nt tem
hc-siRNA e phasiRNA no mesmo espaço, e a classificação em duas classes exigiria mudar
o formato de saída, o `classplots` e os validadores. Reportar a composição primeiro
permite medir quanto disso existe antes de decidir se vale reestruturar.

*(Nota: em 2026-08-02 medimos heterogeneidade ESPACIAL — variação do tamanho dominante
ao longo do locus — e concluímos que loci nucleares eram homogêneos. Esta medição é de
mistura de tamanhos no MESMO espaço, que é outra coisa; as duas conclusões coexistem.)*

## Validação

Nenhuma mudança entra sem passar nos dois reinos:

1. Recall dos 90 piCs (Drosophila) e dos 6 TAS (Arabidopsis) não pode cair.
2. Reads por locus deve subir — é o objetivo.
3. `|frac_top − 0,5|` do `siRNA` deve cair (a assinatura bifita precisa reaparecer),
   e o `s_fractop_bal` da branch `evidence-scores` deve deixar de ser ~0.
4. Enriquecimento em piC (Drosophila) e em TE (Arabidopsis) deve subir ou ficar estável.
5. Pureza de tamanho não deve cair (a fusão não pode criar quimeras).

## Depois disto

- Revalidar a branch `evidence-scores` sobre os loci novos e mesclá-la.
- Testar o `pad` em mais organismos. `piRNAclusterDB` (pirnaclusterdb.org, 51 espécies,
  GTF) serve para artrópodes e vertebrados, com a ressalva de que são clusters
  **preditos por proTRAC**, não validados — servem para testar generalização entre
  clados, não para provar acerto.
