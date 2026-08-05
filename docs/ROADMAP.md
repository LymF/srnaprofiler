# srnaprofiler — Roadmap / Próximos passos

**Autor:** Lucas
**Atualizado:** 2026-08-04

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

**Antes → depois** (Drosophila DRR351382, 1 RPM, pad 500):

| | antes (pad 75) | depois |
|---|---|---|
| span mediano | 88 bp | 1.025 bp |
| reads por locus | 34 | 122 |
| taxa de review | 79% | 10,3% |
| recall de piC | 85/90 | 85/90 |
| massa de piRNA nos piCs | 3,0× | 5,8× |
| siRNA \|frac_top−0,5\| | **0,262** | **0,113** |

O último número é o que destravou a branch: os loci de siRNA passaram a ser bifita, então o
teste de equivalência do `s_fractop_bal` tem sobre o que se pronunciar.

Arabidopsis no mesmo ajuste: **TAS 6/6**, reads por locus 34 → 197, review 33,1% → 2,6%.

Varredura que fixou os valores (as duas em 2026-08-04, regra nova):

| mincov (RPM) | loci | reads/locus | review | piC |
|---|---|---|---|---|
| 0,5 | 228.329 | 41 | 21,5% | 87/90 |
| **1,0** | 171.322 | 55 | 15,0% | 85/90 |
| 2,0 | 116.010 | 75 | 12,2% | 78/90 |
| 5,0 | 52.565 | 116 | 11,1% | 65/90 |

| pad (a 1 RPM) | reads/locus | review | piC |
|---|---|---|---|
| 200 | 55 | 15,0% | 85/90 |
| **500** | 122 | 10,3% | 85/90 |
| 1000 | 228 | 8,4% | 81/90 |

**A preferência de `pad` inverteu-se**: na varredura antiga (regra de vão) o pad 500 custava
2 piCs em Drosophila; com extensão bilateral + gate de perfil ele não custa nenhum. Os dois
reinos convergiram para o mesmo valor, então o `pad` por organismo continua no config mas
hoje não separa nada — mantido porque a calibração é por organismo por natureza.

**Sobre `size_purity`:** mediana 0,19 em Drosophila e 0,63 em Arabidopsis, com o mesmo
código. Não é fusão indevida: dentro da corrida do Drosophila a pureza é plana por tamanho
de locus (0,21 em 150–300 bp, 0,15 acima de 1.200 bp) e por abundância. É a dispersão de
tamanho do piRNA (24–31 nt) contra o pico único de 24 nt da planta. O 0,64 citado antes foi
medido em loci de 41 bp medianos, onde pureza alta é trivial — não é comparável.

### O termo de fita mede o construto errado (aberto, 2026-08-04)

`s_fractop_bal` é um TOST de equivalência: afirma que `frac_top` está em 0,5 ± 0,15. Mas
**bifita não quer dizer 50/50** — um cluster genuinamente bidirecional pode ser 70/30 por
viés de transcrição ou de carregamento em AGO. Medido nos loci novos: 82,8% dos loci de
siRNA têm fita minoritária ≥ 0,2 (são bidirecionais), mas só 60,9% caem dentro da margem
de equivalência, e 26,4% dos bidirecionais são reprovados **só por assimetria**.

O construto certo seria unilateral sobre a fita minoritária (H₀: minoritária ≤ 0,10).
Ressalva medida: isso torna o termo MENOS discriminativo entre as classes (0,828 no siRNA
contra 0,625 no piRNA, contra 0,609 × 0,327 do TOST) — a assimetria separa melhor as
classes do que a bidirecionalidade. Ou seja, é preciso decidir antes se o termo existe
para **definir** o siRNA ou para **discriminar** dele o piRNA; hoje ele é usado como
definidor (`defining_evidence: [fractop_bal]`).

Não bloqueia: o `score` é reportado, **não decide** (`decider.py:116` — a classe sai dos
pesos, o review sai de confiança/margem/`min_reads`/`dicer_call`). O efeito é a coluna
`score` ficar pessimista em loci de siRNA assimétricos.

Relacionado: o piRNA não tem `w_fractop` porque `frac_top` é **bimodal** na classe
(flamenco 0,90 somático uni-fita, 42AB 0,50 germinativo bifita) — um peso único quebraria
um dos dois. A saída seria `defining_evidence` com duas rotas (bifita COM ping-pong, ou
uni-fita SEM), que a sintaxe de slots OR já suporta.

**A geometria de dúplex NÃO serve de saída** (testado 2026-08-04, e corrige uma leitura
anterior errada minha): com `overlap_dup_frac` corrigido, a fração no offset `L−2` fica em
0,020–0,023 em TODOS os tamanhos, **abaixo** do uniforme da janela (1/35 = 0,029),
inclusive nos 21 nt. Não há pico de dúplex em Drosophila. Falta checar em planta (21 nt de
DCL4 são dúplex de livro-texto) para separar "é biologia do tecido" de "o pico se dilui no
locus de 1 kb, como acontecia no phasing antes da varredura por janelas".

**Pendência achada na corrida:** `evidence.py:18` emite
`RuntimeWarning: overflow encountered in scalar divide` quando p ≈ 0 (o `1/p` estoura antes
do `min` com o teto). Não altera resultado; limpar.

---

## RETOMAR AQUI (sessão de 2026-08-04)

Branch **`evidence-scores`** — 17 commits à frente de `master`, 166 testes passando,
NÃO mesclada por decisão. Ledger completo em
`.superpowers/sdd/2026-08-02-evidence-scores/progress.md`.

### O que está pronto e esperando

1. **Score de evidência** (spec `2026-08-02-evidence-scores-design.md`, plano
   `2026-08-02-evidence-scores.md`): implementado, revisado, com os 5 achados da
   revisão final corrigidos. `evidence.py` + colunas `score` e `s_*` no TSV.
2. **Clustering**: **FEITO** — ver a seção acima.

### Por que a branch do score está parada

O teste de equivalência (TOST) do `s_fractop_bal` revelou que os loci de `siRNA`
**não são bifita** (`|frac_top−0,5|` mediano 0,262), então o score da classe fica ~0.
O teste está certo; os **loci** é que estão errados. Validar o score sobre loci
fragmentados não prova nada, então: **clustering primeiro, revalidar o score, mesclar.**

### A descoberta que reordenou tudo

Nossos loci têm 88–140 bp medianos; o ShortStack, **no mesmo BAM**, produz 529 bp.
Metade dos vãos entre loci vizinhos tem 75–200 bp — cortamos por buracos pequenos
demais. Consequência: **34 reads por locus contra 444 dele**, o que tira o poder de
todos os testes estatísticos, e destrói a assinatura de fita. O `--pad 75` nunca foi
calibrado; veio da substituição do tiling de 2 kb.

Varredura já feita (mas que **não transfere** depois da extensão bilateral, porque o
significado do `pad` muda):

| pad | Drosophila | Arabidopsis |
|---|---|---|
| 75 | recall 85/90, 28 reads/locus | TAS 6/6, 34 reads/locus |
| 200 | recall 85/90, 36 reads | TAS 6/6, 59 reads |
| 500 | recall 83/90, 61 reads | TAS 6/6, 131 reads |
| 1000 | recall 81/90 | — |

### Ordem da próxima sessão

1. ~~Executar o plano do clustering~~ e ~~recalibrar `mincov`/`pad`~~ — **FEITOS**, seção acima.
2. Revalidar a branch `evidence-scores` sobre os loci novos: rodar o TOST do
   `s_fractop_bal` agora que o `|frac_top−0,5|` do siRNA caiu para 0,113. Então mesclar.
3. Testar o `pad` em mais organismos: `pirnaclusterDB` mudou para **pirnaclusterdb.org**
   (o domínio uni-mainz saiu do ar), 51 espécies, clusters em GTF preditos com proTRAC
   2.4.4. Ressalva: são preditos, não validados — servem para testar generalização entre
   clados, não para provar acerto.

### Pendências menores

- `git push` (chave `~/.ssh/id_ed25519` tem passphrase; precisa de `ssh-add` num terminal seu).
- Minors diferidos no ledger: guarda de `evidence.density_score` (ordem dos checks com
  `lib_size=None`), e o teste de threads não empurra `genome_bp` não-zero pelo pool.
- **18% dos loci têm duas populações de tamanho** (pureza < 0,5), já em `pad=75`. A tarefa 7
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
