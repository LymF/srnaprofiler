# Saída do módulo `mirna`: uma tabela que junta os dois lados, e o layout dos outros subcomandos

**Autor:** Lucas
**Data:** 2026-08-10
**Estado:** desenho aprovado, aguardando plano de implementação

## Problema

O módulo entrega hoje **um arquivo**: o `<prefix>.mirna_loci.bed` que alimenta o
`classify --mirna-loci`. Tudo o mais que a corrida produziu — e ela produz muito — fica
espalhado em diretórios internos, com nomes e estrutura que dependem do motor.

Três defeitos concretos:

**1. Nenhum dos dois lados é utilizável sozinho, e ninguém os junta.** O
`miRScore_results.csv` traz 20 colunas de evidência (`mSeq, mLen, mStart, mStop, msSeq,
msLen, msStart, msStop, precSeq, precLen, result, flags, mismatches, bulges, mReads,
msReads, totReads, numLibraries, precision`), mas suas coordenadas são **dentro do grampo** e
ele não sabe cromossomo nem fita — trabalha sobre sequências nomeadas. Nossas predições têm
`chrom/start/end/strand` e nenhuma evidência. Quem quiser as duas coisas junta na mão.

**2. O layout é o único fora do padrão do projeto.** `qc`, `mapping` e `classify` escrevem
direto em `--outdir` com nomes `<prefix>.<coisa>` (`qc.py:223`, `mapping.py:401`,
`classify.py:207-210`). O `mirna` cria `outdir/mirna/` e, dentro, **duas árvores diferentes
conforme o motor**: `processed/` + `results/` para `mirdeep2`, `mirdp2/` para `mirdeep-p2`,
mais `validation/` nos dois. Quem troca de `--engine` não acha nada no mesmo lugar.

**3. O que a pessoa quer olhar está enterrado.** Os plots de estrutura secundária **já
existem para os dois motores** e não vêm do preditor: quem os gera é o miRScore, em
`RNAplots/{passed,failed}/` (RNAplot do ViennaRNA) e em `strucVis/{passed,failed}/`, este
ligado por padrão (`miRScore:1524`, desligado só com `-nostrucvis`). O strucVis sobrepõe a
**profundidade de reads** na estrutura, que é o que se olha para julgar um locus. Ficam
quatro níveis abaixo do `--outdir`. O miRDeep2 ainda produz um segundo conjunto em
`pdfs_<timestamp>/`, com timestamp no nome do diretório.

## Decisão

Três entregáveis novos, o layout alinhado com o resto da tool, e nada de lógica nova de
biologia.

### A tabela: `<prefix>.mirna.tsv`

**Uma linha por predição do motor**, não só pelas aprovadas. O BED continua levando apenas as
validadas; a tabela é o registro da corrida inteira, com o veredito e o motivo de cada
descarte. Isso não é preciosismo de auditoria: foi exatamente olhar o conjunto reprovado que
mostrou que o miRDeep2 encontra hc-siRNA em planta (41 das 94 predições com maduro de 24 nt),
e esse tipo de diagnóstico é invisível numa tabela só de aprovados.

Convenção de coordenadas: **0-based, semiaberto**, como o `.clusters.tsv` e o BED. Os motores
reportam 1-based inclusivo e a conversão é nossa, num lugar só. O critério é poder cruzar
`mirna.tsv` com `clusters.tsv` sem ajuste mental.

Ausente é `"NA"`, como no `reporter._cell`.

| coluna | origem | vale para |
|---|---|---|
| `name` | predição | os dois motores |
| `chrom`, `start`, `end`, `strand` | predição | os dois |
| `engine` | `--engine` | os dois |
| `mature_seq`, `mature_len` | predição | os dois |
| `star_seq`, `star_len` | miRScore (`msSeq`, `msLen`) | os dois, se validado |
| `precursor_seq`, `precursor_len` | miRScore (`precSeq`, `precLen`) | os dois, se validado |
| `mature_start`, `mature_end`, `star_start`, `star_end` | miRScore, posições **no precursor** | os dois, se validado |
| `result` | ver abaixo | sempre preenchida |
| `flags` | miRScore | os dois, se validado |
| `mismatches`, `bulges` | miRScore | os dois, se validado |
| `mature_reads`, `star_reads`, `total_reads` | miRScore (`mReads`, `msReads`, `totReads`) | os dois, se validado |
| `precision` | miRScore | os dois, se validado |
| `mirdeep_score`, `mirdeep_prob` | `result_*.csv` | só `mirdeep2` |
| `rfam_alert`, `randfold` | `result_*.csv` | só `mirdeep2` |
| `mirbase_hit`, `mirbase_same_seed` | `result_*.csv` | só `mirdeep2`, e só com `--mature`/`--hairpin` |

O `precursor_seq` publicado é o que **nós** recortamos do genoma com flanco de 25 nt
(`hairpin_with_flank`) e o miRScore validou — é a sequência que corresponde às coordenadas da
linha, não a que o motor emitiu rente ao maduro.

**`result` tem quatro valores, e nenhum deles é vazio.** O ROADMAP registra três ocorrências
do mesmo anti-padrão nesta base (o `0.0` do `pingpong_frac` sob a guarda de 10 pares, o score
de miRNA, o `NaN` na `auc()`): um valor que significa "não medido" sendo lido como medição.
Aqui isso se resolve na origem:

| valor | significa |
|---|---|
| `Pass` | o miRScore aprovou |
| `Fail` | o miRScore reprovou; `flags` diz por quê |
| `NotScored` | o miRScore rodou, mas esta predição não pôde ser submetida — o precursor não foi recortável do genoma (`hairpin_with_flank` devolveu `None`) |
| `NotValidated` | o miRScore não rodou nesta corrida (`--no-validate`, sem `-g/--genome`, ou binário ausente) |

`NotScored` existe porque hoje essas predições somem em silêncio entre `len(predictions)` e
`kept` (`mirna.py:1051-1065`). Sem esse valor elas apareceriam como `Fail`, o que é afirmar
uma reprovação que ninguém fez.

### Os três FASTA — só os validados

`<prefix>.mirna.mature.fa`, `<prefix>.mirna.star.fa`, `<prefix>.mirna.precursor.fa`, com o
`name` da tabela como identificador, os três com **exatamente o mesmo conjunto de nomes do
BED**.

Só as aprovadas, por dois motivos. O precursor recortado corresponde às coordenadas do BED,
então os quatro arquivos falam do mesmo conjunto e podem ser usados juntos sem filtro. E
FASTA de precursor reprovado é material que volta a circular como se fosse real — quem
precisar da sequência de um reprovado tem `precursor_seq` na linha da tabela.

### Os plots: `<prefix>.mirna_structures/{passed,failed}/`

Promovidos a entregável, os dois vereditos separados — auditar visualmente um descarte é a
única coisa que a tabela não faz. Fonte: o `strucVis/` do miRScore, que é o plot com
profundidade; quando ele não existir (binário ausente, ou `-nostrucvis`), cai para o
`RNAplots/`. Se nenhum dos dois existir, o diretório sai vazio e o log avisa em **WARN** —
não é erro fatal, mas também não é silêncio.

**Sem conversão de formato.** Cada ferramenta emite o que emite (`.ps` do strucVis, `.pdf` do
miRDeep2) e nós não convertemos: ganhar uniformidade cosmética não vale uma dependência de
ghostscript no `setup`, que é o módulo com o maior custo de fragilidade do projeto.

### O layout

```
outdir/
├── at.mirna.tsv                 # todas as predições + veredito + evidência
├── at.mirna_loci.bed            # só as validadas → --mirna-loci
├── at.mirna.mature.fa           # ┐
├── at.mirna.star.fa             # ├ só as validadas, mesmos nomes do BED
├── at.mirna.precursor.fa        # ┘
├── at.mirna.run.json            # motor, versões, contagens, o que rodou
├── at.mirna_structures/
│   ├── passed/
│   └── failed/
└── at.mirna_work/
    ├── engine/                  # mirdeep2 OU mirdp2, mesmo caminho
    └── validation/              # miRScore cru, incluindo RNAplots/
```

Entregáveis achatados no `--outdir` com o prefixo, como os outros três subcomandos. O bruto
vai para um diretório de trabalho só, e **`engine/` tem nome fixo** — a árvore deixa de
depender do `--engine`.

O `run.json` segue o `write_run_json` do `classify` e registra: motor e versão, se o miRScore
rodou e sua versão, reino passado a ele, e as contagens (`predicted`, `submitted`, `passed`,
`failed`, `not_scored`). É o que responde "esta tabela veio de uma corrida validada?" sem
abrir a tabela.

**Nada é apagado nesta versão.** A regra do ROADMAP para o SAM (só apagar depois de o
artefato final existir e ser verificado; nunca o que o usuário trouxe) continua valendo, mas
as cópias de input e o FASTA colapsado são pequenos perto do que o `mapping` já movimenta.
Limpeza entra depois, com `--keep-work`, se virar problema medido.

## Escopo

**Dentro:** a tabela, os três FASTA, o diretório de estruturas, o `run.json`, o layout novo, e
o `NotScored`/`NotValidated` no lugar de campo vazio.

**Fora, e cada um por um motivo:**

- **Perfil de isomiR** — decidido em 2026-08-10 que não entra: existe ferramenta dedicada
  (sRNAbench, QuagmiR, IsoMiRmap) e formato padrão (mirtop), e a pergunta que ele responde é
  de biologia de miRNA, não de classificação de locus.
- **MFE/AMFE/MFEI** — o miRScore calcula a MFE e descarta (`miRScore:679`, `:944`, `:1369`);
  os critérios dele são a geometria do dúplex, mais específica que o sumário energético do
  grampo inteiro.
- **Match nosso contra o miRBase** — o cálculo do miRDeep2 é lookup de seed 2–8 contra
  espécies relacionadas (`miRDeep2_core_algorithm.pl:382-400`), e seed é conceito animal;
  reproduzi-lo no motor de planta importaria a premissa errada. As colunas de miRBase são
  repasse, preenchidas só com `--engine mirdeep2` e referência informada.
- **Limpeza do `work/`** e conversão de `.ps` para PDF.

## Critérios de aceite

Medidos na mesma corrida da instalação limpa (Arabidopsis, reads filtradas de ncRNA), cujos
números o ROADMAP já fixa: **87 predições do miRDeep-P2, 56 validadas, 56/56 em miRBase**.

1. `at.mirna.tsv` tem **87 linhas** de dado, e `result == "Pass"` em **56**.
2. `at.mirna_loci.bed` tem **56** intervalos, e o conjunto de nomes é **idêntico** ao das
   linhas `Pass`.
3. Os três FASTA têm **56** registros cada, com o mesmo conjunto de nomes do BED.
4. `at.mirna_structures/` tem um plot para cada predição que o miRScore **pontuou** — as
   `Pass` em `passed/` e as `Fail` em `failed/`. As `NotScored` não têm plot, por construção:
   elas nunca chegaram ao miRScore. Se nenhuma fonte de plot existir, o WARN explica.
5. Rodando com `--no-validate`: 87 linhas, **todas** `result == "NotValidated"`, nenhum FASTA
   escrito, BED escrito com o WARN de não validado que já existe.
6. Os dois motores produzem **os mesmos nomes de arquivo** e a mesma árvore.
7. A suíte segue verde (455 no ambiente limpo), com os testes de caminho atualizados.

## Riscos

**Caminhos mudam, e há teste que depende deles.** O BED sai de `outdir/mirna/` para
`outdir/`. Mitigação: atualizar os testes junto e o log dizer onde cada coisa foi escrita. Não
há compatibilidade a manter — o BED é passado à mão para o `classify`.

**As colunas do miRScore podem mudar entre versões.** Mitigação: ler por **nome de coluna** no
cabeçalho e **falhar alto** se uma coluna necessária faltar, em vez de gravar `NA` e seguir.
É o padrão desta base: o `simulate_rules.py` aborta quando falta coluna justamente porque
assumir 0 em silêncio já inventou 3.080 mudanças de classe uma vez.

**`strucVis` ausente no PATH.** O miRScore o chama por padrão; se faltar, os plots com
profundidade não saem. Mitigação: cair para `RNAplots/` e avisar; o `setup --check` ganha uma
linha para o binário.

**Junção por nome.** Tabela e miRScore casam pelo `name` da predição. Se um motor emitir nomes
repetidos, a junção silenciosamente duplica. Mitigação: verificar unicidade dos nomes ao
montar a tabela e falhar alto se houver colisão.

## Testes

- **Unidade — junção:** predições + `miRScore_results.csv` sintético produzem uma linha por
  predição, com as colunas do miRScore nas linhas certas.
- **Unidade — os quatro `result`:** miRScore ausente ⇒ todas `NotValidated`; predição não
  submetida ⇒ `NotScored`; aprovada ⇒ `Pass`; reprovada ⇒ `Fail` com `flags` preenchida.
- **Unidade — coerência do conjunto:** os nomes do BED, dos três FASTA e das linhas `Pass` são
  o mesmo conjunto.
- **Unidade — coordenadas:** uma predição 1-based inclusiva conhecida vira o intervalo
  0-based semiaberto esperado na tabela e no BED.
- **Unidade — colisão de nomes:** predições com `name` repetido fazem a montagem falhar alto.
- **Unidade — coluna faltando** no CSV do miRScore faz a leitura falhar alto.
- **Integração:** a corrida de Arabidopsis reproduz os sete critérios de aceite.

## Referências

miRScore: a rapid and precise microRNA validation tool. *PLOS Computational Biology*.
https://doi.org/10.1371/journal.pcbi.1013663 — código em https://github.com/Aez35/miRScore
(autores a completar a partir do artigo; o README da ferramenta não os lista).
