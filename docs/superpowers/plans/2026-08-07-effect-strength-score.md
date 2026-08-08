# Score por Tamanho de Efeito — Plano de Implementação

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Trocar `score` de significância (`−log10(p)`) para força de efeito com incerteza em `[0,1]`, de modo que score maior signifique evidência mais forte.

**Architecture:** Uma função pura calcula o limite inferior de Wilson normalizado pelo nulo. `build_record` a aplica às cinco evidências que aparecem em `defining_evidence`, gravando `e_*` no `ClusterRecord` ao lado das `s_*` existentes. `defining_scores` passa a ler `e_*`. As `s_*` permanecem inalteradas como p-valores.

**Tech Stack:** Python, scipy (já dependência), pytest.

## Global Constraints

- Código, comentários, docstrings, mensagens de log e ajuda de CLI em **inglês**.
- Comentários curtos — só o necessário. Contexto vai na mensagem de commit.
- Nunca adicionar `Co-Authored-By` nem crédito de IA em commit.
- As colunas `s_*` **não mudam** de fórmula nem de significado.
- O `reporter` **não ganha colunas**: `e_*` ficam no `ClusterRecord`, fora do TSV.
- A sintaxe do `rules.yaml` (`defining_evidence`, listas aninhadas = OU) **não muda**.
- Nenhum locus pode trocar de classe. Se trocar, é bug.

---

## File Structure

| arquivo | responsabilidade |
|---|---|
| `evidence.py` | `effect_strength` — função pura do limite inferior normalizado |
| `tests/test_evidence_effect_strength.py` (novo) | testes da função pura |
| `signals.py` | `PhaseResult` ganha os componentes que a força da fase precisa |
| `features.py` | calcula `e_1u`, `e_5a`, `e_fractop_bid`, `e_pingpong`, `e_phase` |
| `records.py` | campos `e_*` no `ClusterRecord` |
| `decider.py` | `defining_scores` lê `e_*`; `RuleScorer.__init__` valida os nomes |
| `docs/FORMULAS.md`, `README.md`, `docs/ROADMAP.md` | documentação |

---

### Task 1: a função pura

**Files:**
- Modify: `evidence.py`
- Test: `tests/test_evidence_effect_strength.py` (novo)

**Interfaces:**
- Produces: `evidence.effect_strength(k, n, p0, conf=0.95) -> float` em `[0.0, 1.0]`

- [ ] **Step 1: Write the failing test**

```python
"""Strength of an effect, not its significance.

A p-value mixes effect size with sample size: with a big n, a tiny bias gives
an astronomical p. Measured on the piRNA loci, `score` correlated 0.629 with
distinct_seqs and only 0.323 with the effect, and loci scoring 50+ sat at 24.9%
inside known piRNA clusters, below the 28.3% background.

This returns the Wilson lower bound of k/n, rescaled by how much room the null
leaves above it, so 0.61 reads as "at least 61% of the way from the null to
certainty, at 95% confidence".
"""
import pytest

from evidence import effect_strength


def test_the_same_effect_scores_alike_at_different_n():
    """The property the p-value score does not have."""
    small = effect_strength(38, 50, 0.25)
    large = effect_strength(760, 1000, 0.25)
    assert abs(large - small) < 0.25


def test_more_data_never_lowers_the_strength():
    """A tighter interval can only move the lower bound up."""
    assert effect_strength(760, 1000, 0.25) >= effect_strength(38, 50, 0.25)


def test_effect_at_the_null_is_zero():
    assert effect_strength(250, 1000, 0.25) == 0.0


def test_effect_below_the_null_is_zero():
    assert effect_strength(100, 1000, 0.25) == 0.0


def test_a_perfect_effect_with_plenty_of_data_approaches_one():
    assert effect_strength(2000, 2000, 0.25) > 0.95


def test_a_perfect_effect_with_little_data_does_not():
    """Three out of three is not evidence of much."""
    assert effect_strength(3, 3, 0.25) < 0.7


def test_the_null_rescales_the_range():
    """Same observed proportion is worth more against a harder null."""
    assert effect_strength(600, 1000, 0.10) > effect_strength(600, 1000, 0.50)


def test_guards():
    assert effect_strength(10, 0, 0.25) == 0.0
    assert effect_strength(-1, 100, 0.25) == 0.0
    assert effect_strength(101, 100, 0.25) == 0.0
    assert effect_strength(50, 100, 1.0) == 0.0
    assert effect_strength(50, 100, None) == 0.0
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_evidence_effect_strength.py -q`
Expected: FAIL, `ImportError: cannot import name 'effect_strength'`

- [ ] **Step 3: Write the implementation**

Em `evidence.py`, após `score_from_p`:

```python
def effect_strength(k, n, p0, conf=0.95):
    """How far the effect clears its null, conservatively, in [0, 1].

    Wilson lower bound of k/n at `conf`, rescaled by the room the null leaves:
    (L - p0) / (1 - p0). Unlike -log10(p) this does not grow with n -- it grows
    with precision and saturates at the observed effect -- so the same biology
    scores the same in a small locus and a large one.

    Wilson rather than the normal approximation: our proportions run extreme and
    our samples are often small, where the normal interval leaves [0, 1].
    """
    if p0 is None or n is None:
        return 0.0
    n = int(round(n))
    k = int(round(k))
    if n <= 0 or k <= 0 or k > n or not (0.0 <= p0 < 1.0):
        return 0.0
    z = stats.norm.isf((1.0 - conf) / 2.0)
    phat = k / n
    denom = 1.0 + z * z / n
    centre = (phat + z * z / (2 * n)) / denom
    half = z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / denom
    lower = centre - half
    return round(max(0.0, (lower - p0) / (1.0 - p0)), 4)
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `python -m pytest tests/test_evidence_effect_strength.py -q`
Expected: PASS, 8 passed

- [ ] **Step 5: Commit**

```bash
git add evidence.py tests/test_evidence_effect_strength.py
git commit -m "feat(evidence): effect strength as a normalised Wilson lower bound

A p-value mixes effect size with sample size, so the class score correlated
0.629 with distinct_seqs and 0.323 with the effect, and its top band sat below
the background inside known piRNA clusters.

This grows with precision rather than with n and saturates at the observed
effect, so the same biology scores the same in a small locus and a large one.
Wilson rather than the normal approximation because our proportions run extreme
and our samples are often small."
```

---

### Task 2: expor os componentes da fase

**Files:**
- Modify: `signals.py`
- Test: `tests/test_signals_phase_test.py`

**Interfaces:**
- Consumes: nada
- Produces: `PhaseResult` com os campos novos `draws`, `observed`, `null_frac` (todos com default, para não quebrar construções existentes)

**Contexto:** `PhaseResult` hoje carrega só `score, number, ratio, pvalue`. A força da fase precisa de `observed / draws` (a proporção observada) e de `sucessos / população` (o nulo), que hoje ficam dentro de `_phase_window` e se perdem.

- [ ] **Step 1: Write the failing test**

Acrescentar a `tests/test_signals_phase_test.py`:

```python
def test_result_carries_what_the_effect_size_needs():
    """The hypergeometric inputs are discarded today, and e_phase needs them."""
    pos = [(1000 + c * 21, 4.0) for c in range(8)]
    r = phase_test(pos, region_start=1000, phase_length=21, cycles=11)
    assert r.draws > 0
    assert 0 < r.observed <= r.draws
    assert 0.0 < r.null_frac < 1.0
    # the null is the chance a random position lands in register, ~1/phase_length
    assert abs(r.null_frac - 1 / 21) < 0.02


def test_an_empty_locus_carries_neutral_components():
    r = phase_test([], region_start=0, phase_length=21, cycles=11)
    assert r.draws == 0 and r.observed == 0 and r.null_frac == 0.0
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_signals_phase_test.py -q`
Expected: FAIL, `AttributeError: 'PhaseResult' object has no attribute 'draws'`

- [ ] **Step 3: Write the implementation**

Em `signals.py`, trocar a definição:

```python
PhaseResult = namedtuple(
    "PhaseResult",
    ["score", "number", "ratio", "pvalue", "draws", "observed", "null_frac"],
    defaults=(0, 0, 0.0),
)
```

Em `_phase_window`, guardar os valores que já são calculados:

```python
def _phase_window(positions, region_start, phase_length, cycles):
    """Phasing metrics and p-value for a single window."""
    score, number, ratio = phase_metrics(positions, region_start, phase_length)
    window_len = cycles * phase_length
    population = max(2 * window_len - 1, 1)
    successes = max(2 * cycles - 1, 1)
    draws = len({p for p, _ in positions})
    observed = min(number, successes)
    pvalue = hypergeometric_phase_pvalue(
        population=population, phased_positions=successes,
        draws=draws, observed=observed)
    return PhaseResult(score, number, ratio, pvalue,
                       draws, observed, successes / population)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_signals_phase_test.py tests/test_signals_phase_null.py tests/test_signals_phasing.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add signals.py tests/test_signals_phase_test.py
git commit -m "refactor(signals): PhaseResult carries the hypergeometric inputs

The window's draws, observed count and null proportion are computed and then
discarded, and the effect size needs all three. Added with defaults so the
existing constructions keep working."
```

---

### Task 3: calcular as forças por locus

**Files:**
- Modify: `records.py`, `features.py`
- Test: `tests/test_features_effect_strength.py` (novo)

**Interfaces:**
- Consumes: `evidence.effect_strength` (Task 1), `PhaseResult.draws/observed/null_frac` (Task 2)
- Produces: `ClusterRecord.e_1u`, `.e_5a`, `.e_fractop_bid`, `.e_pingpong`, `.e_phase` (float, default 0.0)

- [ ] **Step 1: Write the failing test**

```python
"""Each defining evidence gets a strength alongside its p-value.

They are not TSV columns: the decider reads them to build `score`.
"""
from features import build_record
from samio import Aln

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _reads(seq, n, start=100, strand="+"):
    return [Aln("c", start + i * 3, strand, seq, 1) for i in range(n)]


def test_a_strong_1u_locus_gets_a_high_strength():
    r = build_record(("c", 0, 400), _reads("T" + "ACGA" * 6, 60), 1_000_000, CFG)
    assert r.e_1u > 0.5


def test_a_locus_without_1u_gets_zero():
    r = build_record(("c", 0, 400), _reads("G" + "ACGA" * 6, 60), 1_000_000, CFG)
    assert r.e_1u == 0.0


def test_strength_does_not_grow_with_the_locus():
    """The p-value score does; that is the bug this replaces."""
    small = build_record(("c", 0, 400), _reads("T" + "ACGA" * 6, 40), 1_000_000, CFG)
    big = build_record(("c", 0, 4000), _reads("T" + "ACGA" * 6, 400), 1_000_000, CFG)
    assert abs(big.e_1u - small.e_1u) < 0.25
    assert big.s_1u > small.s_1u * 3      # the p-value score, for contrast


def test_bidirectional_strength_needs_both_strands():
    uni = build_record(("c", 0, 400), _reads("T" + "ACGA" * 6, 60), 1_000_000, CFG)
    dual = build_record(("c", 0, 400),
                        _reads("T" + "ACGA" * 6, 30)
                        + _reads("T" + "ACGA" * 6, 30, start=200, strand="-"),
                        1_000_000, CFG)
    assert dual.e_fractop_bid > uni.e_fractop_bid
    assert uni.e_fractop_bid == 0.0


def test_an_empty_locus_has_zero_strengths():
    r = build_record(("c", 0, 400), [], 1_000_000, CFG)
    for attr in ("e_1u", "e_5a", "e_fractop_bid", "e_pingpong", "e_phase"):
        assert getattr(r, attr) == 0.0
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_features_effect_strength.py -q`
Expected: FAIL, `AttributeError: 'ClusterRecord' object has no attribute 'e_1u'`

- [ ] **Step 3: Write the implementation**

Em `records.py`, logo após `s_fractop_bid`:

```python
    # Effect strengths, in [0,1] -- what `score` is built from. Not TSV columns.
    e_1u: float = 0.0
    e_5a: float = 0.0
    e_fractop_bid: float = 0.0
    e_pingpong: float = 0.0
    e_phase: float = 0.0
```

Em `features.py`, importar a função:

```python
from evidence import (base_bias_score, strand_bias_score, strand_balance_score,
                      strand_bidirectional_score, effect_strength,
                      pingpong_score_empirical, score_from_p)
```

Logo após a linha `r.s_fractop_bid = ...`:

```python
    # Same k, n and null as the p-values above, read as an effect size instead.
    n = r.distinct_seqs
    r.e_1u = effect_strength(r.pct_1u * n, n, null_u)
    r.e_5a = effect_strength(r.pct_5a * n, n, null_a)
    r.e_fractop_bid = effect_strength(min(r.frac_top, 1 - r.frac_top) * n, n, 0.10)
```

Logo após `r.s_pingpong = pingpong_score_empirical(overlap_counts)`:

```python
    bg = [overlap_counts.get(o, 0.0) for o in range(1, 21) if o != 10]
    if total_pairs > 0 and bg and sum(bg) > 0:
        r.e_pingpong = effect_strength(overlap_counts.get(10, 0), total_pairs,
                                       (sum(bg) / len(bg)) / total_pairs)
```

Dentro do bloco `if r.dicer_call != "N":`, logo após `r.s_phase = score_from_p(...)`:

```python
        r.e_phase = effect_strength(res.observed, res.draws, res.null_frac)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_features_effect_strength.py tests/test_features_scores.py tests/test_features_core.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add records.py features.py tests/test_features_effect_strength.py
git commit -m "feat(features): effect strength for each defining evidence

Reuses the k, n and null already computed for the p-values, read as an effect
size instead. The strengths stay on the record rather than in the TSV: the
decider builds `score` from them and the s_* columns still answer the
significance question."
```

---

### Task 4: o decider passa a usar as forças, e valida os nomes

**Files:**
- Modify: `decider.py`
- Test: `tests/test_decider_class_score.py`

**Interfaces:**
- Consumes: `ClusterRecord.e_*` (Task 3)
- Produces: `record.score` em `[0,1]`; `RuleScorer(config)` levanta `KeyError` para nome de evidência desconhecido

**Contexto:** `defining_scores` faz `getattr(record, f"s_{name}", 0.0)` — nome errado no config vira score 0 em silêncio. O `extras/simulate_rules.py` reconstrói o registro com só 7 features cruas e **não** tem as `e_*`; o default 0.0 tem que continuar valendo para o *atributo*, e a validação recai sobre o **nome no vocabulário**, na construção do scorer.

- [ ] **Step 1: Write the failing test**

Acrescentar a `tests/test_decider_class_score.py`:

```python
def test_score_is_built_from_the_strengths_not_the_p_values():
    r = _rec(dicer_call="24", s_5a=99.0, s_fractop=99.0, e_5a=0.7, e_fractop=0.4)
    RuleScorer(CFG).decide(r)
    assert r.score == 0.4


def test_an_unknown_evidence_name_fails_loudly():
    bad = {"decision": CFG["decision"],
           "classes": {"hcsiRNA": {"size": {"lo": 23, "hi": 24}, "w_size": 1.0,
                                   "defining_evidence": ["5a", "fractopp"]}}}
    with pytest.raises(KeyError) as e:
        RuleScorer(bad)
    assert "fractopp" in str(e.value)


def test_a_record_without_the_strength_attribute_still_scores_zero():
    """extras/simulate_rules.py rebuilds records with raw features only."""
    import types
    r = types.SimpleNamespace(dicer_call="24", frac_top=0.5, pct_1u=0.1,
                              pct_5a=0.8, pingpong_frac=0.0, phase_score=0.0,
                              reads=500)
    RuleScorer(CFG).decide(r)
    assert r.score == 0.0
```

E no topo do arquivo, acrescentar `import pytest` (ele não está lá hoje — o arquivo
importa só `ClusterRecord` e `RuleScorer`). Atualizar também o docstring do arquivo, que
ainda diz *"Score 2 therefore reads as 'every defining evidence has p <= 0.01'"* — na escala
nova isso é falso.

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_decider_class_score.py -q`
Expected: FAIL — o score sai dos `s_*`, e `RuleScorer` não valida nome nenhum

- [ ] **Step 3: Write the implementation**

Em `decider.py`, no topo da classe:

```python
KNOWN_EVIDENCE = ("1u", "5a", "fractop", "fractop_bal", "fractop_bid",
                  "pingpong", "phase")
```

Em `RuleScorer.__init__`, após guardar o config:

```python
        for name, spec in (config.get("classes") or {}).items():
            for entry in spec.get("defining_evidence", []):
                for ev in (entry if isinstance(entry, list) else [entry]):
                    if ev not in KNOWN_EVIDENCE:
                        raise KeyError(
                            f"{name}: unknown evidence {ev!r} in "
                            f"defining_evidence; known: {', '.join(KNOWN_EVIDENCE)}")
```

Em `defining_scores`, trocar o prefixo:

```python
        scores = []
        for name in spec.get("defining_evidence", []):
            if isinstance(name, list):
                scores.append(max(getattr(record, f"e_{n}", 0.0) for n in name))
            else:
                scores.append(getattr(record, f"e_{name}", 0.0))
        return scores
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_decider_class_score.py tests/test_decider.py tests/test_config_defining_evidence.py -q`
Expected: PASS

- [ ] **Step 5: Rodar a suíte inteira**

Run: `python -m pytest -q`
Expected: PASS. Testes que afirmavam `score` na escala de p-valor precisam ser atualizados para a escala `[0,1]` — atualize a expectativa, não a fórmula.

- [ ] **Step 6: Commit**

```bash
git add decider.py tests/test_decider_class_score.py
git commit -m "feat(decider): score from effect strength, and reject unknown evidence names

The score is now the minimum of the defining evidences' strengths, in [0,1], so
a larger locus with a weak bias no longer outranks a small one with a strong
bias. The combination rule and the rules.yaml syntax are unchanged.

defining_scores turned an unknown evidence name into 0.0 in silence, which made
a whole class look evidence-free after a typo. The name is now validated when
the scorer is built. The attribute lookup keeps its default, because
extras/simulate_rules.py legitimately rebuilds records without the strengths."
```

---

### Task 5: validar nos dois reinos e documentar

**Files:**
- Modify: `docs/FORMULAS.md`, `README.md`, `docs/ROADMAP.md`

- [ ] **Step 1: Rodar Drosophila e Arabidopsis**

```bash
SC=/tmp/claude-1000/-home-lucas-srnaprofiler/50ec75a9-4868-4173-841f-ad43b4a8ce1f/scratchpad
OUT=/mnt/d/srnaprofiler-teste/dros/classify
python srnaprofiler.py classify --sam ~/dros_k20.sorted.bam -o $OUT -p dros_eff \
  --mincov 1.0 --threads 2 --config rules.yaml \
  -g /mnt/d/srnaprofiler-teste/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna \
  --mirna-loci $SC/dme_mirbase.bed

POUT=/mnt/d/srnaprofiler-teste/plant/classify
python srnaprofiler.py classify --sam /home/lucas/plant_data/at_flower.sorted.bam \
  -o $POUT -p at_eff --threads 1 --config rules.plant.yaml \
  -g /mnt/d/srnaprofiler-teste/plant/GCF_000001735.4_TAIR10.1_genomic.fna \
  --mirna-loci $SC/ath_mirbase.bed
```

- [ ] **Step 2: Verificar o critério de aceite obrigatório — ordenação**

```bash
python - <<'PY'
import pandas as pd, numpy as np, sys
sys.path.insert(0, "extras")
from validate_piclusters import load_reference
d = pd.read_csv("/mnt/d/srnaprofiler-teste/dros/classify/dros_eff.clusters.tsv",
                sep="\t", low_memory=False)
ref = load_reference("extras/piC_ref_dm6.bed")
ins = np.zeros(len(d), bool)
for r in ref.itertuples():
    ins |= ((d.chrom == r.chrom) & (d.start < r.end) & (d.end > r.start)).values
d["in_piC"] = ins
pi = d[d.cls == "piRNA"]
print("fracao dentro de piC por faixa de score (tem que SUBIR):")
for lo, hi in ((0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 1.01)):
    s = pi[(pi.score >= lo) & (pi.score < hi)]
    if len(s):
        print(f"  {lo}-{hi}: n={len(s):>6,} | {s.in_piC.mean():.1%}")
print(f"  (fundo {d.in_piC.mean():.1%})")
print(f"corr score x distinct_seqs: {np.corrcoef(pi.score, pi.distinct_seqs)[0,1]:.3f}")
print(f"corr score x pct_1u:        {np.corrcoef(pi.score, pi.pct_1u)[0,1]:.3f}")
PY
```

**Aceite:** a fração não pode cair na faixa mais alta (hoje cai de 37,3% para 24,9%), e a correlação com `pct_1u` tem que superar a com `distinct_seqs` (hoje 0,323 contra 0,629). Se falhar, **parar e reportar** — o desenho está errado e não entra.

- [ ] **Step 3: Verificar o piso — nada pode ter mudado de classe**

```bash
python extras/validate_piclusters.py /mnt/d/srnaprofiler-teste/dros/classify/dros_eff.clusters.tsv | sed -n '3p;6,7p'
python extras/validate_mirna.py /mnt/d/srnaprofiler-teste/dros/classify/dros_eff.clusters.tsv dme.gff3 | sed -n '3p;6p'
python extras/validate_plant_tas.py /mnt/d/srnaprofiler-teste/plant/classify/at_eff.clusters.tsv \
  /mnt/d/srnaprofiler-teste/plant/GCF_000001735.4_TAIR10.1_genomic.gff | sed -n '/TAS loci with/p'
python extras/validate_mirna.py /mnt/d/srnaprofiler-teste/plant/classify/at_eff.clusters.tsv ath.gff3 | sed -n '3p;6p'
python - <<'PY'
import pandas as pd
for a, b in (("dros/classify/dros_gn2", "dros/classify/dros_eff"),
             ("plant/classify/at_gn2", "plant/classify/at_eff")):
    x = pd.read_csv(f"/mnt/d/srnaprofiler-teste/{a}.clusters.tsv", sep="\t", low_memory=False)
    y = pd.read_csv(f"/mnt/d/srnaprofiler-teste/{b}.clusters.tsv", sep="\t", low_memory=False)
    m = x[["chrom","start","cls"]].merge(y[["chrom","start","cls"]],
                                         on=["chrom","start"], suffixes=("_a","_b"))
    print(f"{a}: classe mudou em {(m.cls_a != m.cls_b).mean():.2%}  (tem que ser 0.00%)")
PY
```

**Aceite:** piC 85/90, TAS 6/6, miRBase 134/134 e 244/244, e **0,00%** de mudança de classe nos dois.

- [ ] **Step 4: Atualizar a documentação**

Em `docs/FORMULAS.md`, na seção 4, trocar a linha do score por:

```
score = min( e_* das evidências em defining_evidence )      ← 0 se dicer_call == "N"
e_X   = max(0, (L95(p̂, n) − p₀) / (1 − p₀))                ← L95 = limite inferior de Wilson
```

e acrescentar à seção 3 uma tabela com o `p̂`, `n` e `p₀` de cada `e_*` (os mesmos das `s_*`).

Na seção 5 (limitações), **remover** o item "o `score` não ordena confiabilidade" e pôr no lugar as duas ressalvas que permanecem: `n` do ping-pong são pares (não independentes) e a janela da fase é escolhida por máximo (efeito enviesado para cima).

No `README.md`, substituir o aviso "`score` is not a reliability ranking" pela leitura nova da escala.

No `docs/ROADMAP.md`, registrar os números medidos no Step 2 e Step 3.

- [ ] **Step 5: Commit**

```bash
git add docs/FORMULAS.md README.md docs/ROADMAP.md
git commit -m "docs: score now ranks by effect strength

Records the measured ordering against the reference piRNA clusters, the
correlations with effect and with sample size, and that no locus changed class
in either kingdom."
```

---

## Self-Review

**Cobertura da spec:** fórmula e função pura — Task 1; componentes da fase que faltavam — Task 2; as cinco evidências — Task 3; troca no decider e combinação inalterada — Task 4; endurecimento do `defining_scores` — Task 4; critérios de aceite (ordenação, piso, correlação) — Task 5 Steps 2 e 3; documentação das limitações que ficam — Task 5 Step 4.

**Placeholders:** nenhum; todo passo de código traz o código.

**Consistência de tipos:** `effect_strength(k, n, p0, conf=0.95) -> float` é usada com essa assinatura nas Tasks 3; `PhaseResult.draws/observed/null_frac` (Task 2) são consumidos exatamente com esses nomes na Task 3; os cinco atributos `e_*` declarados em `records.py` (Task 3) são lidos por `defining_scores` via `f"e_{name}"` (Task 4), e os nomes batem com o vocabulário `KNOWN_EVIDENCE`.

**Risco a sinalizar para quem executar:** a Task 4 muda a escala do `score`, então testes existentes que afirmam valores na escala de p-valor (`tests/test_decider_class_score.py`, `tests/test_decider.py`) vão falhar. Atualize a **expectativa** para a escala `[0,1]` — se a tentação for mexer na fórmula para o teste antigo passar, o desenho está sendo revertido sem querer.

**Ponto de parada obrigatório:** se o Step 2 da Task 5 mostrar que a ordenação continua não-monotônica, **não seguir para o Step 4**. O desenho falhou o próprio critério e a decisão volta para o usuário.
