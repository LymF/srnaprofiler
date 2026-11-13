"""Task 3: support rule, capped confidence nudge and biology table."""
from pathlib import Path

import yaml

from annotation import Feature, IntervalIndex, Layer, apply_annotation, support_for
from records import ClusterRecord

REPO = Path(__file__).resolve().parents[2]
TABLE = {
    "piRNA": {"supports": ["te"], "conflicts": ["gene_exon"]},
    "miRNA": {"supports": ["intergenic", "gene_intron"], "conflicts": ["te"]},
}
CFG = {"nudge": 0.05, "min_overlap": 0.5, "classes": TABLE}


def _layer(kind, feats, chrom="chr1"):
    return Layer(kind, "x", {chrom: IntervalIndex(feats)}, len(feats))


def _rec(cls="piRNA", conf=0.80, review=False):
    return ClusterRecord(chrom="chr1", start=10, end=30, strand_call="+",
                         cls=cls, confidence=conf, review=review)


TE = [_layer("te", [Feature(0, 100, "+", "feature", "LTR")])]
EXON = [_layer("gene", [Feature(0, 100, "+", "exon", "g1")])]


def test_support_for_table():
    assert support_for("piRNA", {"te"}, TABLE) == "supports"
    assert support_for("piRNA", {"gene_exon"}, TABLE) == "conflicts"
    assert support_for("piRNA", {"te", "gene_exon"}, TABLE) == "supports"
    assert support_for("piRNA", set(), TABLE) == "neutral"
    assert support_for("foo", {"te"}, TABLE) == "neutral"


def test_apply_nudges_and_flags():
    a, b = _rec(), _rec()
    apply_annotation([a], TE, CFG)
    apply_annotation([b], EXON, CFG)
    assert (a.confidence, a.annotation_support, a.review) == (0.85, "supports", False)
    assert (b.confidence, b.annotation_support, b.review) == (0.75, "conflicts", True)
    assert a.cls == b.cls == "piRNA"
    assert a.annotation.startswith("te:")


def test_nudge_is_bounded():
    a, b = _rec(conf=0.98), _rec(conf=0.02)
    apply_annotation([a], TE, CFG)
    apply_annotation([b], EXON, CFG)
    assert a.confidence == 1.0 and a.annotation_support == "supports"
    assert b.confidence == 0.0 and b.annotation_support == "conflicts"


def test_review_never_cleared():
    r = _rec(review=True)
    apply_annotation([r], TE, CFG)
    assert r.review is True


def test_miRNA_and_ambiguous_untouched():
    for cls, want in (("miRNA", "conflicts"), ("ambiguous", "neutral")):
        r = _rec(cls=cls, conf=0.7)
        apply_annotation([r], TE, CFG)
        assert r.confidence == 0.7 and r.review is False
        assert r.annotation.startswith("te:")
        assert r.annotation_support == want


def test_no_layers_marks_none():
    recs = [_rec(), _rec(cls="miRNA")]
    assert apply_annotation(recs, [], CFG) == {
        "none": 2, "supports": 0, "conflicts": 0, "neutral": 0}
    for r in recs:
        assert r.annotation_support == "none" and r.annotation == ""


def test_shipped_configs_have_annotation_block():
    for name in ("rules.yaml", "rules.plant.yaml"):
        cfg = yaml.safe_load((REPO / name).read_text())
        ann = cfg["annotation"]
        assert ann["min_overlap"] == 0.5 and ann["nudge"] == 0.05
        for c in ann["classes"]:
            assert c == "miRNA" or c in cfg["classes"]
