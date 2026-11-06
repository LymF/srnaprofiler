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
    alns = _reads(140, 3) + _reads(100, 2) + _reads(180, 2)
    recs, _ = _run(monkeypatch, alns, [(90, 150), (150, 210)],
                   {"trim_to_reads": True})
    assert len(recs) == 2
    assert recs[0].end <= recs[1].start
    assert recs[0].start >= 90 and recs[1].end <= 210


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
