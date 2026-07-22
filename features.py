"""Extração de features por cluster → ClusterRecord."""
import bisect
from collections import Counter, defaultdict
from records import ClusterRecord
from samio import read_alignments
from signals import pingpong_zscore, phase_metrics, hypergeometric_phase_pvalue

_COMP = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def build_record(region, region_alns, lib_size, config):
    """Constrói um ClusterRecord a partir dos reads JÁ atribuídos à região.

    `region_alns` deve conter todos os reads que sobrepõem a região por span
    (atribuição feita em `extract`); build_record não refiltra por posição.
    """
    chrom, start, end = region
    dmin, dmax = config["dicermin"], config["dicermax"]
    scut = config["strand_cutoff"]
    r = ClusterRecord(chrom=chrom, start=start, end=end, name=f"{chrom}:{start}-{end}")

    local = region_alns
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

    # --- sinais de biogênese ---
    plus5 = Counter()
    minus5 = Counter()
    phase_positions = []
    for a in local:
        L = len(a.seq)
        if a.strand == "+":
            five = a.pos
            phase5 = five
        else:
            five = a.pos + L - 1
            # Correção do overhang 3' de 2 nt: os produtos fásicos de Dicer têm
            # overhang 3' de 2 nt, então o 5' da fita - fica +2 em relação ao
            # registro da fita + no mesmo corte. Subtrair 2 para juntar as duas
            # fitas no MESMO registro de fase. (sinal -2 fixado pela biologia/teste.)
            phase5 = five - 2
        (plus5 if a.strand == "+" else minus5)[five] += 1
        phase_positions.append((phase5, 1))
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
        # População do hipergeométrico limitada a ~10 ciclos de fase (Chen 2007 /
        # PhasiHunter), não à ilha inteira — evita população absurda em ilhas grandes.
        wlen = min(end - start, 10 * pl)
        population = max(2 * wlen - 1, 1)
        r.phase_pvalue = hypergeometric_phase_pvalue(
            population=population,
            phased_positions=pl, draws=n_pos, observed=r.phase_number)
    return r


def _assign_reads(regions, alns):
    """Atribui cada read às regiões que ele sobrepõe POR SPAN, em um passe.

    Overlap por span: `a.pos < region.end and a.pos + len(a.seq) > region.start`
    — consistente com a cobertura construída em `clustering.cluster_regions`.
    As regiões de `cluster_regions` são disjuntas e ordenadas por cromossomo, o
    que permite localizar o bloco de regiões sobrepostas por busca binária em vez
    de rescanear todos os reads por região (evita O(regiões × reads)).
    """
    buckets = [[] for _ in regions]
    per_chrom = defaultdict(list)   # chrom -> [(start, end, orig_index)] ordenado por start
    for i, (c, s, e) in enumerate(regions):
        per_chrom[c].append((s, e, i))
    for c in per_chrom:
        per_chrom[c].sort()
    starts = {c: [t[0] for t in rs] for c, rs in per_chrom.items()}

    for a in alns:
        rs = per_chrom.get(a.chrom)
        if not rs:
            continue
        aend = a.pos + len(a.seq)
        st = starts[a.chrom]
        # regiões com start < aend são as índices [0, j); percorra de trás para
        # frente enquanto end > a.pos (ends também são crescentes → bloco contíguo).
        j = bisect.bisect_left(st, aend)
        k = j - 1
        while k >= 0:
            s, e, oi = rs[k]
            if e <= a.pos:
                break
            buckets[oi].append(a)
            k -= 1
    return buckets


def extract(sam_path, regions, lib_size, config):
    alns = list(read_alignments(sam_path))
    buckets = _assign_reads(regions, alns)
    return [build_record(region, buckets[i], lib_size, config)
            for i, region in enumerate(regions)]
