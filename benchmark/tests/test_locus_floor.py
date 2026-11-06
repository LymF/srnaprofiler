"""Locus floor, read-bounded edges and fragments in process_chromosome."""
from collections import namedtuple
from types import SimpleNamespace

import features
from reporter import write_fragments_tsv

CAln = namedtuple("CAln", "chrom pos seq nh strand")


def _reads(pos, n, length=22, nh=1):
    return [CAln("chr1", pos, "A" * length, nh, "+") for _ in range(n)]


def _run(monkeypatch, alns, merged, cfg=None, carve=()):
    monkeypatch.setattr(features, "read_region", lambda *a, **k: iter(alns))
    monkeypatch.setattr(features, "merge_compatible", lambda *a, **k: merged)
    monkeypatch.setattr(features, "build_record",
                        lambda region, bucket, *a, **k: SimpleNamespace(
                            start=region[1], end=region[2],
                            reads=sum(1.0 / x.nh for x in bucket)))
    base = {"strand_cutoff": 0.8}
    base.update(cfg or {})
    return features.process_chromosome("x.bam", "chr1", 1, 50, 1e6, base,
                                       chrom_len=1000, carve=carve)


def test_floor_off_by_default(monkeypatch):
    recs, frags = _run(monkeypatch, _reads(200, 1), [(150, 310)])
    assert frags == []
    assert (recs[0].start, recs[0].end) == (150, 310)


def test_trim_to_reads_bounds_locus(monkeypatch):
    alns = _reads(200, 3) + _reads(238, 3)   # reads span 200..260
    recs, _ = _run(monkeypatch, alns, [(150, 310)], {"trim_to_reads": True})
    assert (recs[0].start, recs[0].end) == (200, 260)
    recs, _ = _run(monkeypatch, alns, [(150, 310)])
    assert (recs[0].start, recs[0].end) == (150, 310)


def test_trim_clipped_to_window(monkeypatch):
    # The read at 140-162 overlaps both windows, so it is assigned to both.
    # Unclipped, left would end at 162 and right would start at 140.
    alns = _reads(140, 3) + _reads(100, 2) + _reads(180, 2)
    recs, _ = _run(monkeypatch, alns, [(90, 150), (150, 210)],
                   {"trim_to_reads": True})
    assert [(r.start, r.end) for r in recs] == [(100, 150), (150, 202)]


def test_min_span_makes_fragment(monkeypatch):
    alns = [CAln("chr1", 200 + (i % 9), "A" * 22, 1, "+") for i in range(25)]
    alns[0] = CAln("chr1", 200, "A" * 22, 1, "+")
    alns[1] = CAln("chr1", 208, "A" * 22, 1, "+")   # span 200..230 = 30
    recs, frags = _run(monkeypatch, alns, [(150, 280)], {"locus_min_span": 50})
    assert recs == []
    assert len(frags) == 1
    f = frags[0]
    assert f["reason"] == "min_span" and f["reads"] == 25.0
    assert f["read_span"] == 30
    assert set(f) == {"chrom", "start", "end", "reads", "read_span", "reason"}


def test_min_reads_makes_fragment(monkeypatch):
    alns = [CAln("chr1", 200 + 14 * i, "A" * 22, 1, "+") for i in range(5)]
    _, frags = _run(monkeypatch, alns, [(150, 300)], {"locus_min_reads": 10})
    assert frags[0]["reason"] == "min_reads"
    _, frags = _run(monkeypatch, alns, [(150, 300)],
                    {"locus_min_reads": 10, "locus_min_span": 500})
    assert frags[0]["reason"] == "min_reads,min_span"


def test_carved_mirna_never_fragment_nor_trimmed(monkeypatch):
    alns = _reads(420, 3)
    recs, frags = _run(monkeypatch, alns, [(400, 480)],
                       {"locus_min_reads": 10, "locus_min_span": 50,
                        "trim_to_reads": True}, carve=[(400, 480)])
    assert frags == []
    assert (recs[0].start, recs[0].end) == (400, 480)


def test_all_fragments_returns_empty_records(monkeypatch):
    recs, frags = _run(monkeypatch, _reads(200, 2), [(150, 300)],
                       {"locus_min_reads": 10})
    assert recs == [] and len(frags) == 1


def test_empty_bucket_is_fragment_not_crash(monkeypatch):
    recs, frags = _run(monkeypatch, _reads(200, 2), [(150, 300), (600, 700)],
                       {"locus_min_reads": 1, "trim_to_reads": True})
    assert len(recs) == 1 and frags[0]["reads"] == 0 and frags[0]["read_span"] == 0


def test_write_fragments_tsv_header_only(tmp_path):
    p = tmp_path / "f.tsv"
    write_fragments_tsv([], str(p))
    assert p.read_text() == "chrom\tstart\tend\treads\tread_span\treason\n"


def test_write_fragments_tsv_rows(tmp_path):
    p = tmp_path / "f.tsv"
    write_fragments_tsv([{"chrom": "c", "start": 1, "end": 9, "reads": 2.5,
                          "read_span": 8, "reason": "min_reads"}], str(p))
    assert p.read_text().splitlines()[1] == "c\t1\t9\t2.50\t8\tmin_reads"


def test_run_classify_all_fragments(monkeypatch, tmp_path):
    """Every region a fragment: clusters.tsv is header-only, fragments has rows."""
    import classify
    from reporter import _COLS
    frags = [{"chrom": "chr1", "start": 10, "end": 90, "reads": 3.0,
              "read_span": 40, "reason": "min_reads"},
             {"chrom": "chr2", "start": 5, "end": 50, "reads": 1.5,
              "read_span": 20, "reason": "min_reads,min_span"}]
    monkeypatch.setattr(classify, "load_config", lambda p: {})
    monkeypatch.setattr(classify, "class_size_windows", lambda cfg: [])
    monkeypatch.setattr(classify, "RuleScorer", lambda cfg: None)  # no records to score
    monkeypatch.setattr(classify, "ensure_sorted_indexed_bam", lambda *a: "x.bam")
    monkeypatch.setattr(classify, "prefers_xm_as_nh", lambda b: False)
    monkeypatch.setattr(classify, "mapped_chromosomes", lambda b: ["chr1", "chr2"])
    monkeypatch.setattr(classify, "reference_alignment_counts", lambda b: [4, 2])
    monkeypatch.setattr(classify, "available_memory_mb", lambda: 10 ** 6)
    monkeypatch.setattr(classify, "extract_by_chrom",
                        lambda *a, **k: ([], frags))
    args = SimpleNamespace(outdir=str(tmp_path), config=None, prefix="t",
                           sam="x.sam", norm=1e6, genome="g.fa", mincov=1,
                           mirna_loci=None, organelle_contigs=None, threads=1,
                           pad=50)
    classify.run_classify(args)
    clusters = (tmp_path / "t.clusters.tsv").read_text().splitlines()
    assert clusters == ["\t".join(_COLS)]
    rows = (tmp_path / "t.fragments.tsv").read_text().splitlines()
    assert rows[0] == "chrom\tstart\tend\treads\tread_span\treason"
    assert rows[1:] == ["chr1\t10\t90\t3.00\t40\tmin_reads",
                        "chr2\t5\t50\t1.50\t20\tmin_reads,min_span"]
