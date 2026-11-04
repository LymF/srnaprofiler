import ast
import gzip
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

import benchmark.simulate as sim_module
from benchmark.simulate import SimSpec, simulate
from benchmark.truth import load_regions

COMP = str.maketrans("ACGTN", "TGCAN")
ALL = {"piRNA": 4, "siRNA": 4, "miRNA": 3, "phasiRNA": 3, "hcsiRNA": 4}
BG = [("bg1", "ACGTACGTACGTACGTACGTACGT"), ("bg2", "TTTTGGGGCCCCAAAATTTTGGGG")]


def revcomp(s):
    return s.translate(COMP)[::-1]


@pytest.fixture
def genome_file(tmp_path):
    rng = np.random.default_rng(7)
    seqs = {"chr1": list("".join(rng.choice(list("ACGT"), 9000))),
            "chr2": list("".join(rng.choice(list("ACGT"), 7000)))}
    seqs["chr1"][2000:2300] = "N" * 300
    p = tmp_path / "genome.fa"
    p.write_text("".join(f">{k}\n{''.join(v)}\n" for k, v in seqs.items()))
    return p


@pytest.fixture
def bg_file(tmp_path):
    p = tmp_path / "bg.fastq"
    p.write_text("".join(f"@{n}\n{s}\n+\n{'F' * len(s)}\n" for n, s in BG))
    return p


def read_genome(path):
    out, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            out[name] = ""
        else:
            out[name] += line
    return out


def read_fastq(path):
    lines = gzip.open(path, "rt").read().splitlines()
    return [(lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]


def locate(genome, seq):
    """(chrom, start, end, strand) of the unique genome match of a read."""
    for chrom, g in genome.items():
        for strand, q in (("+", seq), ("-", revcomp(seq))):
            i = g.find(q)
            if i >= 0:
                return chrom, i, i + len(seq), strand
    raise AssertionError(f"read {seq} not found in genome")


def run(tmp_path, genome_file, bg_file, classes=None, seed=1, name="out", **kw):
    return simulate(SimSpec(classes or dict(ALL), seed), bg_file, genome_file, tmp_path / name, **kw)


def sim_reads(result, genome, cls):
    """[(id, seq, chrom, start, end, strand, locus_name)] for one class."""
    out = []
    for rid, seq, _ in read_fastq(result.fastq):
        if rid.startswith(f"sim_{cls}_"):
            c, s, e, st = locate(genome, seq)
            out.append((rid, seq, c, s, e, st, f"{cls}_{rid.split('_')[2]}"))
    return out


def truth_rows(result):
    return [l.split("\t") for l in result.truth_bed.read_text().splitlines()]


def test_simulate_is_deterministic_for_same_seed(tmp_path, genome_file, bg_file):
    a = run(tmp_path, genome_file, bg_file, name="a")
    b = run(tmp_path, genome_file, bg_file, name="b")
    c = run(tmp_path, genome_file, bg_file, seed=2, name="c")
    assert a.fastq.read_bytes() == b.fastq.read_bytes()
    assert a.truth_bed.read_bytes() == b.truth_bed.read_bytes()
    assert a.truth_bed.read_bytes() != c.truth_bed.read_bytes()


def test_truth_coordinates_inside_chromosome(tmp_path, genome_file, bg_file):
    genome = read_genome(genome_file)
    rows = truth_rows(run(tmp_path, genome_file, bg_file))
    assert len(rows) == sum(ALL.values())
    spans = []
    for chrom, start, end, name, cls, strand in rows:
        start, end = int(start), int(end)
        assert 0 <= start < end <= len(genome[chrom])
        assert "N" not in genome[chrom][start:end]
        assert strand in "+-."
        spans.append((chrom, start, end))
    for i, (c1, s1, e1) in enumerate(spans):
        for c2, s2, e2 in spans[i + 1:]:
            assert c1 != c2 or e1 <= s2 or e2 <= s1


def test_pirna_reads_have_10nt_overlap_on_opposite_strands(tmp_path, genome_file, bg_file):
    genome = read_genome(genome_file)
    res = run(tmp_path, genome_file, bg_file, {"piRNA": 3}, reads_per_locus=60)
    reads = {r[0]: r for r in sim_reads(res, genome, "piRNA")}
    partners = [r for rid, r in reads.items() if rid.endswith("pp")]
    assert partners
    for p in partners:
        prim = reads[p[0][:-2]]
        assert p[5] != prim[5]
        assert min(p[4], prim[4]) - max(p[3], prim[3]) == 10
        assert p[1][9] == "A"


def test_simulate_does_not_import_rules():
    src = Path(sim_module.__file__).read_text()
    imported = set()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Import):
            imported |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            imported.add((node.module or "").split(".")[-1])
    assert not imported & {"classify", "decider", "features", "signals", "evidence"}
    assert "rules" not in src


def test_every_read_is_genome_substring_inside_its_locus(tmp_path, genome_file, bg_file):
    genome = read_genome(genome_file)
    res = run(tmp_path, genome_file, bg_file)
    loci = {r[3]: (r[0], int(r[1]), int(r[2])) for r in truth_rows(res)}
    for cls in ALL:
        reads = sim_reads(res, genome, cls)
        assert reads
        for rid, seq, chrom, s, e, strand, locus in reads:
            lc, ls, le = loci[locus]
            assert chrom == lc and ls <= s and e <= le
            assert seq == (genome[chrom][s:e] if strand == "+" else revcomp(genome[chrom][s:e]))
            assert "N" not in seq


def test_phasirna_starts_in_21nt_register(tmp_path, genome_file, bg_file):
    genome = read_genome(genome_file)
    res = run(tmp_path, genome_file, bg_file, {"phasiRNA": 3}, reads_per_locus=80)
    loci = {r[3]: int(r[1]) for r in truth_rows(res)}
    reads = sim_reads(res, genome, "phasiRNA")
    assert {r[5] for r in reads} == {"+", "-"}
    assert all(e - s == 21 for _, _, _, s, e, _, _ in reads)
    # plus reads start at 21k; minus reads sit 2 nt downstream-shifted (2-nt 3' overhang duplex)
    phase = [((s - loci[l]) % 21 == 0) if st == "+" else ((e - loci[l] - 19) % 21 == 0)
             for _, _, _, s, e, st, l in reads]
    assert sum(phase) / len(phase) > 0.8


def test_mirna_duplex_has_2nt_3prime_overhang(tmp_path, genome_file, bg_file):
    genome = read_genome(genome_file)
    res = run(tmp_path, genome_file, bg_file, {"miRNA": 3}, reads_per_locus=50)
    loci = {r[3]: (int(r[1]), int(r[2]), r[5]) for r in truth_rows(res)}
    reads = sim_reads(res, genome, "miRNA")
    for locus, (ls, le, ts) in loci.items():
        counts = Counter((s, e, st) for _, _, _, s, e, st, l in reads if l == locus)
        assert len(counts) == 2
        assert {st for _, _, st in counts} == {ts}
        (five, n5), (three, n3) = counts.most_common()
        assert n5 > n3 and 21 <= five[1] - five[0] <= 23
        if ts == "-":   # mirror into transcript orientation
            five, three = [(ls + le - e, ls + le - s) for s, e, _ in (five, three)]
        else:
            five, three = five[:2], three[:2]
        # base paired with the last base of the 3p arm lies 2 nt 5' of the 5p arm start
        assert (ls + le - 1) - (three[1] - 1) == five[0] - 2


def test_background_reads_present_unchanged(tmp_path, genome_file, bg_file):
    res = run(tmp_path, genome_file, bg_file)
    recs = read_fastq(res.fastq)
    assert recs[:2] == [(n, s, "F" * len(s)) for n, s in BG]
    assert all(r[0].startswith("sim_") and set(r[2]) == {"I"} for r in recs[2:])


def test_per_class_beds_load_with_load_regions(tmp_path, genome_file, bg_file):
    res = run(tmp_path, genome_file, bg_file)
    for cls, n in ALL.items():
        regs = load_regions(res.truth_bed.parent / f"sim_truth.{cls}.bed", cls, "sim")
        assert len(regs) == n and all(r.name.startswith(f"{cls}_") for r in regs)
    assert all(len(r) == 6 for r in truth_rows(res))


def test_unknown_class_raises(tmp_path, genome_file, bg_file):
    with pytest.raises(ValueError):
        run(tmp_path, genome_file, bg_file, {"lncRNA": 1})


def test_adapter_only_on_simulated_reads(tmp_path, genome_file, bg_file):
    res = run(tmp_path, genome_file, bg_file, {"siRNA": 1}, adapter="TGGAATTC", reads_per_locus=10)
    recs = read_fastq(res.fastq)
    assert len(recs) == 2 + 10
    assert recs[0][1] == BG[0][1]
    assert all(r[1].endswith("TGGAATTC") and len(r[2]) == len(r[1]) for r in recs[2:])


def test_class_signatures_5prime_and_length(tmp_path, genome_file, bg_file):
    genome = read_genome(genome_file)
    res = run(tmp_path, genome_file, bg_file, {"piRNA": 3, "siRNA": 3, "hcsiRNA": 3}, reads_per_locus=80)
    pi = [r[1] for r in sim_reads(res, genome, "piRNA") if not r[0].endswith("pp")]
    assert sum(s[0] == "T" for s in pi) / len(pi) > 0.7 and all(24 <= len(s) <= 30 for s in pi)
    si = sim_reads(res, genome, "siRNA")
    assert sum(len(r[1]) == 21 for r in si) / len(si) > 0.6
    assert 0.3 < sum(r[5] == "+" for r in si) / len(si) < 0.7
    hc = sim_reads(res, genome, "hcsiRNA")
    assert sum(len(r[1]) == 24 for r in hc) / len(hc) > 0.6
    assert sum(r[1][0] == "A" for r in hc) / len(hc) > 0.5
    assert {r[5] for r in hc} == {"+", "-"}
