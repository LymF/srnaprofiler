# Motor de Classificação de sRNA (v1) — Plano de Implementação

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Substituir a heurística de classificação por tiles por um motor "medir→decidir" que define clusters por cobertura, extrai sinais robustos por locus e atribui classe + confiança.

**Architecture:** Módulos flat na raiz (padrão do projeto). Fluxo: SAM → `clustering` (ilhas de cobertura) → `features` (ClusterRecord por locus) → `decider` (RuleScorer plugável) → `reporter` (TSV/GFF/JSON). `signals` e `samio` são utilitários puros. O antigo `clusters.py`/`_classify_loci` fica intacto durante a transição.

**Tech Stack:** Python 3, pysam (SAM/BAM), Biopython (FASTA), PyYAML (config), statistics/math (stdlib), pytest.

## Global Constraints

- Linguagem: Python 3; módulos **flat na raiz** do repo (seguir padrão existente: `from utils import log`).
- Métodos reimplementados dos tools de referência — **nunca copiar código** (tinyRNA/pingpongpro/mirdeep2 são GPL).
- Sem numpy: usar `statistics` e `math` da stdlib.
- Campo de classe no dataclass chama-se `cls` (`class` é palavra reservada).
- Não tocar em `clusters.py` (mantido para A/B).
- Autor em qualquer doc = "Lucas" (nunca IA).
- Testes em `tests/`; rodar com `pytest`.

---

### Task 1: Scaffold de testes + `ClusterRecord`

**Files:**
- Create: `records.py`
- Create: `tests/__init__.py` (vazio)
- Test: `tests/test_records.py`

**Interfaces:**
- Produces: `ClusterRecord` dataclass com os campos abaixo e propriedade `length`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_records.py
from records import ClusterRecord

def test_clusterrecord_defaults_and_length():
    r = ClusterRecord(chrom="chr1", start=100, end=250)
    assert r.length == 150
    assert r.cls == ""
    assert r.review is False
    assert r.dicer_call == "N"
    assert r.size_counts == {}
    assert r.class_scores == {}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_records.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'records'`

- [ ] **Step 3: Write minimal implementation**

```python
# records.py
"""Contrato de dados do motor de classificação: um registro por locus."""
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class ClusterRecord:
    chrom: str
    start: int
    end: int
    name: str = ""
    # abundância
    reads: int = 0
    reads_rpm: float = 0.0
    distinct_seqs: int = 0
    major_rna: str = ""
    major_rna_reads: int = 0
    # tamanho
    dicer_call: str = "N"
    frac_in_window: float = 0.0
    short: int = 0
    long: int = 0
    size_counts: dict = field(default_factory=dict)
    # fita
    frac_top: Optional[float] = None
    strand_call: str = "."
    # composição
    pct_1u: float = 0.0
    pct_10a: float = 0.0
    gc: float = 0.0
    # biogênese
    pingpong_z: float = 0.0
    pingpong_pairs10: int = 0
    phase_score: float = 0.0
    phase_pvalue: float = 1.0
    phase_number: int = 0
    phase_ratio: float = 0.0
    overlap19_frac: float = 0.0
    # repetição
    multimap_frac: float = 0.0
    # classificação
    cls: str = ""
    confidence: float = 0.0
    review: bool = False
    class_scores: dict = field(default_factory=dict)
    evidence: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start
```

Também criar `tests/__init__.py` vazio.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_records.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add records.py tests/test_records.py tests/__init__.py
git commit -m "feat: add ClusterRecord data contract"
```

---

### Task 2: Ping-pong z-score (`signals.py`)

**Files:**
- Create: `signals.py`
- Test: `tests/test_signals_pingpong.py`

**Interfaces:**
- Produces: `pingpong_zscore(counts: dict[int,int], target: int = 10, lo: int = 1, hi: int = 20) -> tuple[float, int]` → `(z, pairs_at_target)`. Fundo = todos os offsets em [lo,hi] exceto `target`. Desvio populacional.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_signals_pingpong.py
from signals import pingpong_zscore

def test_pingpong_peak_at_10_gives_high_z():
    counts = {o: 2 for o in range(1, 21)}
    counts[10] = 100
    z, pairs = pingpong_zscore(counts)
    assert pairs == 100
    assert z > 3.0

def test_pingpong_flat_gives_zero_z():
    counts = {o: 5 for o in range(1, 21)}
    z, pairs = pingpong_zscore(counts)
    assert z == 0.0
    assert pairs == 5

def test_pingpong_insufficient_background():
    z, pairs = pingpong_zscore({10: 50})
    assert z == 0.0
    assert pairs == 50
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_signals_pingpong.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'signals'`

- [ ] **Step 3: Write minimal implementation**

```python
# signals.py
"""Sinais de biogênese: ping-pong z-score e phasing."""
import statistics


def pingpong_zscore(counts, target=10, lo=1, hi=20):
    """Z-score do overlap 5'-5' no `target` (default 10 nt) contra o fundo.

    counts: {overlap_length: n_pares}. Fundo = offsets em [lo,hi] menos o target.
    """
    background = [counts.get(o, 0) for o in range(lo, hi + 1) if o != target]
    peak = counts.get(target, 0)
    if len(background) < 2:
        return 0.0, peak
    sd = statistics.pstdev(background)
    if sd == 0:
        return 0.0, peak
    mean = statistics.fmean(background)
    return round((peak - mean) / sd, 3), peak
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_signals_pingpong.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add signals.py tests/test_signals_pingpong.py
git commit -m "feat: add ping-pong z-score signal"
```

---

### Task 3: Phase score + p hipergeométrico (`signals.py`)

**Files:**
- Modify: `signals.py`
- Test: `tests/test_signals_phasing.py`

**Interfaces:**
- Produces:
  - `phase_metrics(positions: list[tuple[int,int]], region_start: int, phase_length: int) -> tuple[float,int,float]` → `(phase_score, phase_number, phase_ratio)`. `positions` = lista de `(pos5, abundance)`.
  - `hypergeometric_phase_pvalue(population: int, phased_positions: int, draws: int, observed: int) -> float`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_signals_phasing.py
import math
from signals import phase_metrics, hypergeometric_phase_pvalue

def test_phase_metrics_perfectly_phased():
    # três reads em fase (registro 0), abundância 10 cada
    positions = [(100, 10), (121, 10), (142, 10)]
    score, number, ratio = phase_metrics(positions, region_start=100, phase_length=21)
    assert number == 3
    assert ratio == 1.0
    assert score == round(1.0 * 3 * math.log(30), 3)

def test_phase_metrics_no_reads():
    assert phase_metrics([], 0, 21) == (0.0, 0, 0.0)

def test_hypergeometric_pvalue_bounds():
    p = hypergeometric_phase_pvalue(population=41, phased_positions=21, draws=5, observed=5)
    assert 0.0 <= p <= 1.0
    # observar 0 em fase nunca é significativo
    assert hypergeometric_phase_pvalue(41, 21, 5, 0) == 1.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_signals_phasing.py -v`
Expected: FAIL with `ImportError: cannot import name 'phase_metrics'`

- [ ] **Step 3: Write minimal implementation**

Adicionar ao final de `signals.py`:

```python
import math
from collections import defaultdict
from math import comb


def phase_metrics(positions, region_start, phase_length):
    """Score de fase estilo PhasiHunter: ratio * number * log(abun) no registro dominante.

    positions: [(pos5, abundance)]. Registro = (pos5 - region_start) % phase_length.
    """
    reg_abun = defaultdict(int)
    reg_cycles = defaultdict(set)
    total = 0
    for pos5, ab in positions:
        r = (pos5 - region_start) % phase_length
        cyc = (pos5 - region_start) // phase_length
        reg_abun[r] += ab
        reg_cycles[r].add(cyc)
        total += ab
    if total == 0:
        return 0.0, 0, 0.0
    dom = max(reg_abun, key=reg_abun.get)
    phase_abun = reg_abun[dom]
    phase_number = len(reg_cycles[dom])
    phase_ratio = phase_abun / total
    score = phase_ratio * phase_number * math.log(phase_abun) if phase_abun > 0 else 0.0
    return round(score, 3), phase_number, round(phase_ratio, 3)


def hypergeometric_phase_pvalue(population, phased_positions, draws, observed):
    """P(X >= observed) sorteando `draws` de `population` com `phased_positions` sucessos.

    comb() devolve 0 para argumentos inválidos, então termos impossíveis somem.
    """
    if observed <= 0:
        return 1.0
    denom = comb(population, draws)
    if denom == 0:
        return 1.0
    p = 0.0
    for i in range(observed, min(phased_positions, draws) + 1):
        p += comb(phased_positions, i) * comb(population - phased_positions, draws - i) / denom
    return min(p, 1.0)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_signals_phasing.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add signals.py tests/test_signals_phasing.py
git commit -m "feat: add phase score and hypergeometric p-value"
```

---

### Task 4: Leitura de SAM (`samio.py`) + clustering por cobertura (`clustering.py`)

**Files:**
- Create: `samio.py`
- Create: `clustering.py`
- Test: `tests/test_clustering.py`

**Interfaces:**
- Produces (`samio.py`): `Aln` = namedtuple(`chrom, pos, strand, seq, nh`); `read_alignments(sam_path) -> Iterator[Aln]` (ignora secundários e não-mapeados; `pos` 0-based leftmost).
- Produces (`clustering.py`): `call_islands(depth_by_pos: dict[int,int], mincov: int, pad: int) -> list[tuple[int,int]]` (start 0-based, end exclusivo); `cluster_regions(sam_path: str, mincov: int, pad: int) -> list[tuple[str,int,int]]`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_clustering.py
from clustering import call_islands

def test_call_islands_merges_within_pad():
    depth = {10: 5, 11: 5, 12: 5,   15: 5, 16: 5}  # gap 12->15 = 2 vazios
    # pad=2 funde os dois blocos; pad=1 separa
    assert call_islands(depth, mincov=1, pad=2) == [(10, 17)]
    assert call_islands(depth, mincov=1, pad=1) == [(10, 13), (15, 17)]

def test_call_islands_respects_mincov():
    depth = {10: 1, 11: 3, 12: 1}
    assert call_islands(depth, mincov=2, pad=0) == [(11, 12)]

def test_call_islands_empty():
    assert call_islands({}, mincov=1, pad=0) == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_clustering.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'clustering'`

- [ ] **Step 3: Write minimal implementation**

```python
# samio.py
"""Leitura de alinhamentos (SAM/BAM) para o motor de classificação."""
from collections import namedtuple
import pysam

Aln = namedtuple("Aln", ["chrom", "pos", "strand", "seq", "nh"])


def read_alignments(sam_path):
    """Itera alinhamentos primários mapeados. Exclui secundários (flag 256)."""
    mode = "rb" if str(sam_path).endswith(".bam") else "r"
    with pysam.AlignmentFile(sam_path, mode) as fh:
        for rec in fh:
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            seq = rec.query_sequence or ""
            if not seq:
                continue
            nh = rec.get_tag("NH") if rec.has_tag("NH") else 1
            strand = "-" if rec.is_reverse else "+"
            yield Aln(rec.reference_name, rec.reference_start, strand, seq, nh)
```

```python
# clustering.py
"""Definição de clusters por ilhas de cobertura (estilo ShortStack)."""
from collections import defaultdict
from samio import read_alignments


def call_islands(depth_by_pos, mincov, pad):
    """Posições com depth>=mincov, fundindo ilhas cujo gap (vazios) <= pad.

    Retorna [(start, end)] com end exclusivo.
    """
    positions = sorted(p for p, d in depth_by_pos.items() if d >= mincov)
    if not positions:
        return []
    islands = []
    start = prev = positions[0]
    for p in positions[1:]:
        if p - prev - 1 > pad:
            islands.append((start, prev + 1))
            start = p
        prev = p
    islands.append((start, prev + 1))
    return islands


def cluster_regions(sam_path, mincov, pad):
    """Constrói cobertura por cromossomo e chama ilhas."""
    depth = defaultdict(lambda: defaultdict(int))
    for a in read_alignments(sam_path):
        for p in range(a.pos, a.pos + len(a.seq)):
            depth[a.chrom][p] += 1
    regions = []
    for chrom in sorted(depth):
        for start, end in call_islands(depth[chrom], mincov, pad):
            regions.append((chrom, start, end))
    return regions
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_clustering.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add samio.py clustering.py tests/test_clustering.py
git commit -m "feat: add SAM reader and coverage-island clustering"
```

---

### Task 5: FeatureExtractor — núcleo (contagem, DicerCall, FracTop, composição)

**Files:**
- Create: `features.py`
- Test: `tests/data/toy.sam` (fixture)
- Test: `tests/test_features_core.py`

**Interfaces:**
- Consumes: `samio.read_alignments`, `records.ClusterRecord`.
- Produces:
  - `revcomp(seq: str) -> str`
  - `build_record(region: tuple[str,int,int], alns: list[Aln], lib_size: Optional[int], config: dict) -> ClusterRecord` (preenche tudo MENOS os sinais de biogênese, que entram na Task 6).
  - `extract(sam_path: str, regions: list, lib_size, config) -> list[ClusterRecord]`.
- `config` usa chaves: `dicermin` (int), `dicermax` (int), `strand_cutoff` (float).

- [ ] **Step 1: Write the failing test**

Criar fixture `tests/data/toy.sam` (janela Dicer 20–24; 4 reads de 21nt na fita +, 1 read de 30nt "long"):

```
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	101	255	21M	*	0	0	TAGCTAGCTAGCTAGCTAGCT	*	NH:i:1
r2	0	chr1	101	255	21M	*	0	0	TAGCTAGCTAGCTAGCTAGCT	*	NH:i:1
r3	0	chr1	105	255	21M	*	0	0	TGGCTAGCTAGCTAGCTAGCT	*	NH:i:1
r4	0	chr1	108	255	21M	*	0	0	TAGCTAGCTAGCTAGCTAGCA	*	NH:i:2
r5	0	chr1	120	255	30M	*	0	0	GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	*	NH:i:1
```

```python
# tests/test_features_core.py
from samio import read_alignments
from features import build_record, revcomp

CONFIG = {"dicermin": 20, "dicermax": 24, "strand_cutoff": 0.8}

def test_revcomp():
    assert revcomp("TAGC") == "GCTA"

def test_build_record_core():
    alns = list(read_alignments("tests/data/toy.sam"))
    region = ("chr1", 100, 160)
    r = build_record(region, alns, lib_size=1_000_000, config=CONFIG)
    assert r.reads == 5
    assert r.frac_top == 1.0
    assert r.strand_call == "+"
    assert r.dicer_call == "21"          # 4/5 = 0.8 na janela
    assert r.short == 0 and r.long == 1  # o 30nt é "long"
    assert round(r.frac_in_window, 2) == 0.8
    assert r.pct_1u == round(4 / 5, 4)    # r1..r4 começam com T(=U); r5 com G
    assert r.reads_rpm == 5.0             # 5/1e6*1e6
    assert r.multimap_frac == round(1 / 5, 4)  # r4 tem NH=2
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_features_core.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'features'`

- [ ] **Step 3: Write minimal implementation**

```python
# features.py
"""Extração de features por cluster → ClusterRecord."""
from collections import Counter
from typing import Optional
from records import ClusterRecord
from samio import read_alignments

_COMP = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _region_alns(region, alns):
    chrom, start, end = region
    out = []
    for a in alns:
        if a.chrom == chrom and start <= a.pos < end:
            out.append(a)
    return out


def build_record(region, alns, lib_size, config):
    chrom, start, end = region
    dmin, dmax = config["dicermin"], config["dicermax"]
    scut = config["strand_cutoff"]
    r = ClusterRecord(chrom=chrom, start=start, end=end, name=f"{chrom}:{start}-{end}")

    local = _region_alns(region, alns)
    r.reads = len(local)
    if r.reads == 0:
        return r

    strand_count = Counter()
    seqs = Counter()
    seq_lens = Counter()
    gc_bases = 0
    total_bases = 0
    n_1u = n_10a = 0
    multimap = 0
    phase_ready = []  # (pos5, oriented_len) usado na Task 6

    for a in local:
        oriented = a.seq if a.strand == "+" else revcomp(a.seq)
        L = len(oriented)
        strand_count[a.strand] += 1
        seqs[oriented] += 1
        if L < dmin:
            seq_lens["short"] += 1
        elif L > dmax:
            seq_lens["long"] += 1
        else:
            seq_lens[str(L)] += 1
        gc_bases += sum(1 for b in oriented if b in "GCgc")
        total_bases += L
        if oriented[:1] in ("T", "U"):
            n_1u += 1
        if L >= 10 and oriented[9] == "A":
            n_10a += 1
        if a.nh > 1:
            multimap += 1

    r.distinct_seqs = len(seqs)
    r.major_rna, r.major_rna_reads = seqs.most_common(1)[0]
    r.frac_top = round(strand_count["+"] / r.reads, 4)
    if r.frac_top >= scut:
        r.strand_call = "+"
    elif r.frac_top <= 1 - scut:
        r.strand_call = "-"
    else:
        r.strand_call = "."
    r.short = seq_lens["short"]
    r.long = seq_lens["long"]
    predom_len, predom_reads = None, -1
    for L in range(dmin, dmax + 1):
        c = seq_lens.get(str(L), 0)
        r.size_counts[L] = c
        if c > predom_reads:
            predom_len, predom_reads = L, c
    in_window = r.reads - r.short - r.long
    r.frac_in_window = round(in_window / r.reads, 4)
    r.dicer_call = str(predom_len) if r.frac_in_window >= 0.8 else "N"
    r.pct_1u = round(n_1u / r.reads, 4)
    r.pct_10a = round(n_10a / r.reads, 4)
    r.gc = round(gc_bases / total_bases, 4) if total_bases else 0.0
    r.multimap_frac = round(multimap / r.reads, 4)
    r.reads_rpm = round(r.reads / lib_size * 1e6, 4) if lib_size else 0.0
    return r


def extract(sam_path, regions, lib_size, config):
    alns = list(read_alignments(sam_path))
    return [build_record(region, alns, lib_size, config) for region in regions]
```

Criar diretório `tests/data/` com o `toy.sam` acima.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_features_core.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add features.py tests/test_features_core.py tests/data/toy.sam
git commit -m "feat: FeatureExtractor core (counting, DicerCall, FracTop, composition)"
```

---

### Task 6: FeatureExtractor — sinais de biogênese (ping-pong, phasing, overlap19)

**Files:**
- Modify: `features.py`
- Test: `tests/data/pingpong.sam` (fixture)
- Test: `tests/test_features_signals.py`

**Interfaces:**
- Consumes: `signals.pingpong_zscore`, `signals.phase_metrics`, `signals.hypergeometric_phase_pvalue`.
- Produces: `build_record` agora também preenche `pingpong_z, pingpong_pairs10, overlap19_frac, phase_score, phase_number, phase_ratio, phase_pvalue`.

- [ ] **Step 1: Write the failing test**

Fixture `tests/data/pingpong.sam` — um par +/− formando overlap 5′-5′ de 10 nt. Geometria: o `+` tem 5′ em 100 (0-based); o `−` (10 nt, começando em 100) tem 5′ em `100+10-1 = 109`; overlap_length = `109-100+1 = 10` → assinatura ping-pong.

```
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
plus	0	chr1	101	255	10M	*	0	0	AAAAAAAAAA	*	NH:i:1
minus	16	chr1	101	255	10M	*	0	0	TTTTTTTTTT	*	NH:i:1
```

```python
# tests/test_features_signals.py
from samio import read_alignments
from features import build_record

CONFIG = {"dicermin": 8, "dicermax": 12, "strand_cutoff": 0.8}

def test_pingpong_pair_detected():
    alns = list(read_alignments("tests/data/pingpong.sam"))
    r = build_record(("chr1", 100, 160), alns, lib_size=None, config=CONFIG)
    assert r.pingpong_pairs10 == 1     # um par + / - com overlap 10
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_features_signals.py -v`
Expected: FAIL (`pingpong_pairs10 == 0`, pois os sinais ainda não são computados)

- [ ] **Step 3: Write minimal implementation**

No topo de `features.py`, adicionar import:

```python
from signals import pingpong_zscore, phase_metrics, hypergeometric_phase_pvalue
```

Antes do `return r` em `build_record`, inserir o bloco de sinais (usa `local`, `r`, `start`):

```python
    # --- sinais de biogênese ---
    plus5 = Counter()
    minus5 = Counter()
    phase_positions = []
    for a in local:
        L = len(a.seq)
        if a.strand == "+":
            five = a.pos
        else:
            five = a.pos + L - 1
        (plus5 if a.strand == "+" else minus5)[five] += 1
        phase_positions.append((five, 1))
    overlap_counts = Counter()
    for a5, ac in plus5.items():
        for d in range(0, 20):            # overlap_length = d+1 em 1..20
            bc = minus5.get(a5 + d, 0)
            if bc:
                overlap_counts[d + 1] += ac * bc
    r.pingpong_z, r.pingpong_pairs10 = pingpong_zscore(overlap_counts)
    total_pairs = sum(overlap_counts.values())
    r.overlap19_frac = round(overlap_counts.get(19, 0) / total_pairs, 4) if total_pairs else 0.0

    if r.dicer_call != "N":
        pl = int(r.dicer_call)
        r.phase_score, r.phase_number, r.phase_ratio = phase_metrics(
            phase_positions, region_start=start, phase_length=pl)
        n_pos = len({p for p, _ in phase_positions})
        r.phase_pvalue = hypergeometric_phase_pvalue(
            population=max((end - start) * 2 - 1, 1),
            phased_positions=pl, draws=n_pos, observed=r.phase_number)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_features_signals.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add features.py tests/test_features_signals.py tests/data/pingpong.sam
git commit -m "feat: wire ping-pong and phasing signals into FeatureExtractor"
```

---

### Task 7: Decisor (`RuleScorer`) + `rules.yaml`

**Files:**
- Create: `decider.py`
- Create: `rules.yaml`
- Test: `tests/test_decider.py`

**Interfaces:**
- Consumes: `records.ClusterRecord`.
- Produces: `load_config(path) -> dict`; `class RuleScorer` com `__init__(self, config: dict)` e `decide(self, record: ClusterRecord) -> ClusterRecord` (preenche `cls, confidence, review, class_scores, evidence`).
- `config` tem `classes` (dict por classe com pesos/janelas) e `decision` (`review_confidence`, `min_margin`).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_decider.py
from records import ClusterRecord
from decider import RuleScorer, load_config

def _cfg():
    return load_config("rules.yaml")

def test_pirna_strong_signal():
    r = ClusterRecord(chrom="c", start=0, end=100, reads=500,
                      dicer_call="27", frac_top=0.95, pct_1u=0.9, pingpong_z=8.0)
    RuleScorer(_cfg()).decide(r)
    assert r.cls == "piRNA"
    assert r.confidence > 0.5
    assert r.review is False

def test_sirna_dualstrand_21nt():
    r = ClusterRecord(chrom="c", start=0, end=100, reads=500,
                      dicer_call="21", frac_top=0.5, pct_1u=0.2, pingpong_z=0.0)
    RuleScorer(_cfg()).decide(r)
    assert r.cls == "siRNA"

def test_weak_locus_flagged_review():
    r = ClusterRecord(chrom="c", start=0, end=100, reads=5,
                      dicer_call="N", frac_top=0.5)
    RuleScorer(_cfg()).decide(r)
    assert r.review is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_decider.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'decider'`

- [ ] **Step 3: Write minimal implementation**

```yaml
# rules.yaml — config de score por organismo (default genérico)
decision:
  review_confidence: 0.6
  min_margin: 0.15
  min_reads: 10
classes:
  siRNA:
    size: {lo: 20, hi: 22}
    w_size: 1.0
    w_fractop_mid: 1.0        # premia FracTop ~0.5 (bifita)
    w_pingpong_penalty: 1.0   # penaliza ping-pong
  piRNA:
    size: {lo: 24, hi: 31}
    w_size: 1.0
    w_pingpong: 1.5
    w_1u: 1.0
```

```python
# decider.py
"""Decisor plugável. v1: RuleScorer (modelo de score transparente)."""
import math
import yaml


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def _sigmoid(x, x0=1.65, scale=1.0):
    return 1.0 / (1.0 + math.exp(-(x - x0) / scale))


def _window(size, lo, hi):
    if size is None:
        return 0.0
    if lo <= size <= hi:
        return 1.0
    d = lo - size if size < lo else size - hi
    return max(0.0, 1.0 - d / 3.0)


def _mid(frac_top):
    if frac_top is None:
        return 0.0
    return math.exp(-((frac_top - 0.5) ** 2) / (2 * 0.15 ** 2))


class RuleScorer:
    def __init__(self, config):
        self.cfg = config

    def _size(self, record):
        return int(record.dicer_call) if record.dicer_call != "N" else None

    def _score_class(self, name, spec, record):
        size = self._size(record)
        s = 0.0
        terms = {}
        if "w_size" in spec:
            t = spec["w_size"] * _window(size, spec["size"]["lo"], spec["size"]["hi"])
            s += t; terms["size"] = round(t, 3)
        if "w_pingpong" in spec:
            t = spec["w_pingpong"] * _sigmoid(record.pingpong_z)
            s += t; terms["pingpong"] = round(t, 3)
        if "w_1u" in spec:
            t = spec["w_1u"] * record.pct_1u
            s += t; terms["1u"] = round(t, 3)
        if "w_fractop_mid" in spec:
            t = spec["w_fractop_mid"] * _mid(record.frac_top)
            s += t; terms["fractop_mid"] = round(t, 3)
        if "w_pingpong_penalty" in spec:
            t = spec["w_pingpong_penalty"] * _sigmoid(record.pingpong_z)
            s -= t; terms["pingpong_penalty"] = round(-t, 3)
        return s, terms

    def decide(self, record):
        scores = {}
        evidence = {}
        for name, spec in self.cfg["classes"].items():
            scores[name], evidence[name] = self._score_class(name, spec, record)
        record.class_scores = {k: round(v, 3) for k, v in scores.items()}
        # softmax para confiança
        ordered = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
        top_name, top_score = ordered[0]
        exps = {k: math.exp(v) for k, v in scores.items()}
        z = sum(exps.values())
        conf = exps[top_name] / z if z else 0.0
        margin = top_score - (ordered[1][1] if len(ordered) > 1 else 0.0)
        dec = self.cfg["decision"]
        record.cls = top_name
        record.confidence = round(conf, 3)
        record.evidence = "; ".join(f"{k}={v}" for k, v in evidence[top_name].items())
        record.review = (
            conf < dec["review_confidence"]
            or margin < dec["min_margin"]
            or record.reads < dec["min_reads"]
            or record.dicer_call == "N"
        )
        return record
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_decider.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add decider.py rules.yaml tests/test_decider.py
git commit -m "feat: add RuleScorer decider and default rules.yaml"
```

---

### Task 8: Reporter (TSV / GFF3 / run.json)

**Files:**
- Create: `reporter.py`
- Test: `tests/test_reporter.py`

**Interfaces:**
- Consumes: `records.ClusterRecord`.
- Produces: `write_tsv(records: list, path)`; `write_gff3(records: list, path)`; `write_run_json(params: dict, path)`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_reporter.py
import json
from records import ClusterRecord
from reporter import write_tsv, write_gff3, write_run_json

def _rec():
    r = ClusterRecord(chrom="chr1", start=100, end=200, reads=50, dicer_call="21")
    r.cls, r.confidence, r.review = "siRNA", 0.8, False
    return r

def test_write_tsv(tmp_path):
    p = tmp_path / "out.tsv"
    write_tsv([_rec()], str(p))
    lines = p.read_text().splitlines()
    assert lines[0].startswith("chrom\tstart\tend")
    assert "siRNA" in lines[1]

def test_write_gff3(tmp_path):
    p = tmp_path / "out.gff3"
    write_gff3([_rec()], str(p))
    txt = p.read_text()
    assert txt.startswith("##gff-version 3")
    assert "chr1\tsrnaprofiler\tsiRNA\t101\t200" in txt  # GFF 1-based

def test_write_run_json(tmp_path):
    p = tmp_path / "run.json"
    write_run_json({"mincov": 5}, str(p))
    assert json.loads(p.read_text())["mincov"] == 5
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_reporter.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'reporter'`

- [ ] **Step 3: Write minimal implementation**

```python
# reporter.py
"""Saídas do motor: TSV (matriz de features), GFF3 e proveniência JSON."""
import json
from dataclasses import fields

_COLS = [
    "chrom", "start", "end", "name", "reads", "reads_rpm", "distinct_seqs",
    "major_rna", "dicer_call", "frac_in_window", "frac_top", "strand_call",
    "pct_1u", "pct_10a", "gc", "pingpong_z", "pingpong_pairs10", "overlap19_frac",
    "phase_score", "phase_pvalue", "phase_number", "phase_ratio", "multimap_frac",
    "cls", "confidence", "review", "evidence",
]


def write_tsv(records, path):
    with open(path, "w") as fh:
        fh.write("\t".join(_COLS) + "\n")
        for r in records:
            fh.write("\t".join(str(getattr(r, c)) for c in _COLS) + "\n")


def write_gff3(records, path):
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for r in records:
            attrs = f"ID={r.name};confidence={r.confidence};review={r.review}"
            fh.write("\t".join([
                r.chrom, "srnaprofiler", r.cls or "unclassified",
                str(r.start + 1), str(r.end), str(r.confidence),
                r.strand_call, ".", attrs,
            ]) + "\n")


def write_run_json(params, path):
    with open(path, "w") as fh:
        json.dump(params, fh, indent=2)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_reporter.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add reporter.py tests/test_reporter.py
git commit -m "feat: add reporter (TSV/GFF3/run.json)"
```

---

### Task 9: Orquestrador `classify.py` + subcomando + miRNA opcional

**Files:**
- Create: `classify.py`
- Modify: `srnaprofiler.py` (adicionar subparser `classify`)
- Test: `tests/test_classify_integration.py`

**Interfaces:**
- Consumes: `clustering.cluster_regions`, `features.extract`, `decider.RuleScorer/load_config`, `reporter.*`.
- Produces: `parse_mirna_bed(path) -> dict[str,list[tuple[int,int]]]`; `mark_mirna(records, mirna_intervals) -> None` (seta `cls="miRNA"` em quem sobrepõe); `run_classify(args)`.
- `args` (atributos): `sam, genome, prefix, outdir, mincov, pad, mirna_loci, config, norm`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_classify_integration.py
import types, os
from classify import run_classify, parse_mirna_bed, mark_mirna
from records import ClusterRecord

def test_mark_mirna_overlap():
    recs = [ClusterRecord("chr1", 100, 200), ClusterRecord("chr1", 500, 600)]
    mark_mirna(recs, {"chr1": [(150, 160)]})
    assert recs[0].cls == "miRNA" and recs[0].review is False
    assert recs[1].cls == ""

def test_run_classify_end_to_end(tmp_path):
    args = types.SimpleNamespace(
        sam="tests/data/toy.sam", genome=None, prefix="t", outdir=str(tmp_path),
        mincov=1, pad=5, mirna_loci=None, config="rules.yaml", norm=1_000_000)
    run_classify(args)
    assert os.path.exists(tmp_path / "t.clusters.tsv")
    assert os.path.exists(tmp_path / "t.clusters.gff3")
    assert os.path.exists(tmp_path / "t.run.json")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_classify_integration.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'classify'`

- [ ] **Step 3: Write minimal implementation**

```python
# classify.py
"""Orquestrador do subcomando `classify`."""
import os
from collections import defaultdict
from clustering import cluster_regions
from features import extract
from decider import RuleScorer, load_config
from reporter import write_tsv, write_gff3, write_run_json
from utils import log


def parse_mirna_bed(path):
    intervals = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.split("\t")
            intervals[f[0]].append((int(f[1]), int(f[2])))
    return intervals


def mark_mirna(records, mirna_intervals):
    for r in records:
        for (s, e) in mirna_intervals.get(r.chrom, []):
            if r.start < e and s < r.end:
                r.cls = "miRNA"
                r.confidence = 1.0
                r.review = False
                break


def run_classify(args):
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    config = load_config(args.config)
    ext_cfg = {"dicermin": 20, "dicermax": 24, "strand_cutoff": 0.8}
    ext_cfg.update(config.get("extractor", {}))

    log("Definindo clusters por cobertura...", "INFO")
    regions = cluster_regions(args.sam, args.mincov, args.pad)
    log(f"{len(regions)} clusters", "INFO")

    records = extract(args.sam, regions, args.norm, ext_cfg)

    mirna_intervals = parse_mirna_bed(args.mirna_loci) if args.mirna_loci else {}
    if mirna_intervals:
        mark_mirna(records, mirna_intervals)

    scorer = RuleScorer(config)
    for r in records:
        if r.cls != "miRNA":
            scorer.decide(r)

    base = os.path.join(outdir, args.prefix)
    write_tsv(records, base + ".clusters.tsv")
    write_gff3(records, base + ".clusters.gff3")
    write_run_json(vars(args), base + ".run.json")
    log(f"Escrito: {base}.clusters.tsv / .gff3 / .run.json", "INFO")
```

Em `srnaprofiler.py`, adicionar dentro de `build_parser()` (após o bloco do subparser `mirna`, antes de `return parser`):

```python
    # ==================== Classify Module ====================
    from classify import run_classify
    cl = subparsers.add_parser("classify", help="Classificar clusters de sRNA (motor medir->decidir)")
    cl.add_argument("-p", "--prefix", required=True)
    cl.add_argument("-o", "--outdir", required=True)
    cl.add_argument("--sam", required=True, help="SAM/BAM mapeado no genoma inteiro")
    cl.add_argument("-g", "--genome", help="FASTA do genoma (opcional no v1)")
    cl.add_argument("--mincov", type=int, default=2)
    cl.add_argument("--pad", type=int, default=75)
    cl.add_argument("--mirna-loci", dest="mirna_loci", help="BED de miRNA de alta confiança")
    cl.add_argument("--config", default="rules.yaml")
    cl.add_argument("--norm", type=int, help="Tamanho da biblioteca para RPM")
    cl.set_defaults(func=run_classify)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_classify_integration.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add classify.py srnaprofiler.py tests/test_classify_integration.py
git commit -m "feat: add classify orchestrator, subcommand and optional miRNA overlap"
```

---

### Task 10: Remover tiling do `qc`

**Files:**
- Modify: `qc.py` (função de preparação de genoma)
- Modify: `srnaprofiler.py` (args do subparser `qc`)
- Test: `tests/test_qc_no_tiling.py`

**Interfaces:**
- Produces: `prepare_genome(genome_path, outdir, prefix) -> str` que **padroniza** o FASTA (uppercase, quebra de linha consistente) e retorna o caminho, **sem tilear**.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_qc_no_tiling.py
from qc import prepare_genome

def test_prepare_genome_does_not_tile(tmp_path):
    g = tmp_path / "g.fa"
    g.write_text(">chr1\n" + "ACGT" * 1000 + "\n")   # 4000 nt, 1 sequência
    out = prepare_genome(str(g), str(tmp_path), "t")
    headers = [l for l in open(out) if l.startswith(">")]
    assert len(headers) == 1                      # continua 1 sequência (não tileou em janelas)
    assert headers[0].strip() == ">chr1"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_qc_no_tiling.py -v`
Expected: FAIL — a função atual tileia (gera múltiplos headers) ou tem outra assinatura.

- [ ] **Step 3: Write minimal implementation**

Em `qc.py`, localizar a função de preparação/tiling de genoma (a que responde ao `--kmer`/`--no-tiling`) e substituir por uma versão que só padroniza. Adicionar (ou substituir) esta função:

```python
def prepare_genome(genome_path, outdir, prefix):
    """Padroniza o FASTA (uppercase) sem tilear. Retorna o caminho de saída."""
    import os
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{prefix}.genome.fa")
    with open(genome_path) as fin, open(out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(line)
            else:
                fout.write(line.strip().upper() + "\n")
    return out
```

No `run_qc` (mesmo arquivo), trocar a chamada de tiling pela chamada a `prepare_genome`. Em `srnaprofiler.py`, no subparser `qc`, remover os argumentos `--kmer` e `--no-tiling` (o genoma não é mais tileado).

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_qc_no_tiling.py -v`
Expected: PASS

- [ ] **Step 5: Run full suite + commit**

Run: `pytest tests/ -v`
Expected: todos PASS

```bash
git add qc.py srnaprofiler.py tests/test_qc_no_tiling.py
git commit -m "refactor: remove genome tiling from qc (map to full genome)"
```

---

## Notas de execução

- **Dependências** (garantir no ambiente): `pip install pysam pyyaml biopython pytest`.
- **Validação de regressão (manual, pós-plano):** rodar `classify` em `files/finalPIRNA.sam` e comparar as chamadas com o `clusters` antigo; divergências são esperadas (bugs corrigidos) e devem ser anotadas.
- **Fora do v1 (specs seguintes):** pós-miRDeep (curar `result_*.csv` → alimentar `--mirna-loci`), `MLDecider` sobre `*.clusters.tsv`, clustering por DBSCAN, plots por classe no `reporter`.
