"""Task 2: annotation layers and per-locus overlap summary."""
import pytest

import annotation
from annotation import (Feature, IntervalIndex, Layer, load_layers,
                        parse_layer_spec, summarize_locus)


def _layer(kind, feats, chrom="chr1"):
    return Layer(kind, "x", {chrom: IntervalIndex(feats)}, len(feats))


def test_parse_layer_spec_ok_and_errors():
    assert parse_layer_spec("te:a.out") == ("te", "a.out")
    with pytest.raises(ValueError, match="gene, te, ncrna, phas"):
        parse_layer_spec("foo:a.out")
    with pytest.raises(ValueError):
        parse_layer_spec("te:")


def test_te_context_and_text():
    L = _layer("te", [Feature(0, 100, "-", "feature", "LTR/Gypsy")])
    r = summarize_locus("chr1", 10, 30, "+", [L], 0.5)
    assert r["contexts"] == {"te"}
    assert r["text"] == "te:LTR/Gypsy:antisense:1.00"


def test_below_min_overlap_is_not_a_context():
    L = _layer("te", [Feature(0, 13, "+", "feature", "X")])
    r = summarize_locus("chr1", 10, 20, "+", [L], 0.5)
    assert r == {"contexts": set(), "text": ""}


def test_union_not_sum():
    L = _layer("te", [Feature(0, 14, "+", "feature", "A"),
                      Feature(10, 16, "+", "feature", "B")])
    r = summarize_locus("chr1", 10, 20, "+", [L], 0.5)
    assert r["text"].endswith(":0.60")


def test_gene_exon_vs_intron_vs_intergenic():
    G = _layer("gene", [Feature(0, 100, "+", "gene", "g1"),
                        Feature(0, 40, "+", "exon", "")])
    assert summarize_locus("chr1", 10, 30, "+", [G], 0.5)["contexts"] == {"gene_exon"}
    assert summarize_locus("chr1", 60, 80, "+", [G], 0.5)["contexts"] == {"gene_intron"}
    assert summarize_locus("chr1", 200, 220, "+", [G], 0.5)["contexts"] == {"intergenic"}
    T = _layer("te", [Feature(0, 5, "+", "feature", "X")])
    assert summarize_locus("chr1", 200, 220, "+", [T], 0.5)["contexts"] == set()


def test_bed_gene_layer_counts_as_exon_like():
    G = _layer("gene", [Feature(0, 100, ".", "feature", "")])
    r = summarize_locus("chr1", 10, 30, "+", [G], 0.5)
    assert r["contexts"] == {"gene_exon"}


def test_unknown_contig_locus_is_empty():
    L = _layer("te", [Feature(0, 100, "+", "feature", "X")])
    assert summarize_locus("chrX", 10, 30, "+", [L], 0.5) == {
        "contexts": set(), "text": ""}


def test_load_layers_skips_unknown_contigs_and_warns(tmp_path, monkeypatch):
    p = tmp_path / "x.bed"
    p.write_text("chr1\t0\t10\nscaf9\t0\t10\n")
    calls = []
    monkeypatch.setattr(annotation, "log", lambda m, lvl="INFO": calls.append((m, lvl)))
    layers = load_layers([f"te:{p}"], {"chr1"})
    assert list(layers[0].index) == ["chr1"]
    assert layers[0].n_features == 1
    assert len(calls) == 1 and calls[0][1] == "WARN"
    assert "1" in calls[0][0] and "scaf9" in calls[0][0]


def test_intergenic_not_claimed_on_unannotated_contig():
    G = _layer("gene", [Feature(0, 100, "+", "gene", "g1")])
    T = _layer("te", [Feature(0, 100, "+", "feature", "X")])
    r = summarize_locus("chrZ", 10, 30, "+", [G, T], 0.5)
    assert r == {"contexts": set(), "text": ""}


def test_intergenic_in_gap_on_annotated_contig():
    G = _layer("gene", [Feature(0, 50, "+", "gene", "g1"),
                        Feature(150, 200, "+", "gene", "g2")])
    r = summarize_locus("chr1", 80, 120, "+", [G], 0.5)
    assert r["contexts"] == {"intergenic"}


def test_relation_both_for_unstranded_feature_and_dot_for_unstranded_locus():
    L = _layer("te", [Feature(0, 100, ".", "feature", "X")])
    assert summarize_locus("chr1", 10, 30, "+", [L], 0.5)["text"] == "te:X:both:1.00"
    L2 = _layer("te", [Feature(0, 100, "+", "feature", "X")])
    assert summarize_locus("chr1", 10, 30, ".", [L2], 0.5)["text"] == "te:X:.:1.00"


def test_coverage_equal_to_min_overlap_counts():
    L = _layer("te", [Feature(0, 15, "+", "feature", "X")])
    r = summarize_locus("chr1", 10, 20, "+", [L], 0.5)
    assert r["contexts"] == {"te"} and r["text"].endswith(":0.50")


def test_load_layers_warning_details(tmp_path, monkeypatch):
    p = tmp_path / "y.bed"
    rows = ["chr1\t0\t10"] + [f"s{i}\t0\t10" for i in range(5)] + ["s0\t20\t30"]
    p.write_text("\n".join(rows) + "\n")
    calls = []
    monkeypatch.setattr(annotation, "log", lambda m, lvl="INFO": calls.append((m, lvl)))
    load_layers([f"te:{p}"], {"chr1"})
    assert len(calls) == 1
    msg = calls[0][0]
    assert "skipped 6 features" in msg
    assert sum(f"s{i}" in msg for i in range(5)) == 3
