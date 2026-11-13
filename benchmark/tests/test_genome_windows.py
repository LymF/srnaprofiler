"""Task 2: window features and gene models for the genome-browser view."""
from collections import namedtuple

import report_genome
from annotation import Feature, IntervalIndex, Layer
from report_genome import FEATURE_CAP, build_windows, window_features, window_for


def _layer(kind, feats, chrom="chr1"):
    return Layer(kind, "x", {chrom: IntervalIndex(feats)}, len(feats))


def test_window_for_clamps_and_zero_flank():
    assert window_for(100, 200, 50, 1000) == (50, 250)
    assert window_for(10, 200, 50, 1000) == (0, 250)
    assert window_for(900, 980, 50, 1000) == (850, 1000)
    assert window_for(100, 200, 0, 1000) == (100, 200)


def _ensembl():
    return [
        Feature(100, 1000, "+", "gene", "gene", "gene:G1", "", "Gname"),
        Feature(100, 900, "+", "mRNA", "mRNA", "transcript:T1", "gene:G1", "T1n"),
        Feature(100, 1000, "+", "mRNA", "mRNA", "transcript:T2", "gene:G1", "T2n"),
        # T1 (short): one exon
        Feature(100, 900, "+", "exon", "exon", "e1", "transcript:T1"),
        # T2 (long): UTR, CDS, exon, UTR
        Feature(100, 200, "+", "five_prime_UTR", "five_prime_UTR", "", "transcript:T2"),
        Feature(200, 400, "+", "CDS", "CDS", "", "transcript:T2"),
        Feature(100, 400, "+", "exon", "exon", "e2", "transcript:T2"),
        Feature(700, 1000, "+", "exon", "exon", "e3", "transcript:T2"),
        Feature(700, 800, "+", "CDS", "CDS", "", "transcript:T2"),
        Feature(800, 1000, "+", "three_prime_UTR", "three_prime_UTR", "", "transcript:T2"),
        Feature(150, 152, "+", "start_codon", "start_codon", "", "transcript:T2"),
    ]


def test_ensembl_gene_one_entry_longest_transcript():
    out = window_features("chr1", 0, 2000, (300, 350), [_layer("gene", _ensembl())])
    assert len(out["g"]) == 1
    label, strand, s, e, hl, blocks = out["g"][0]
    assert (label, strand, s, e, hl) == ("Gname", "+", 100, 1000, 1)
    assert sorted(blocks, key=lambda b: (b[1], b[2], b[0])) == blocks
    assert ["E", 100, 400] in blocks and ["U", 800, 1000] in blocks
    assert all(b[0] in "ECU" for b in blocks)  # start codon ignored
    assert len(blocks) == 6
    assert out["ly"] == ["gene"] and out["mo"]["g"] == 0


def test_gene_blocks_outside_window_still_complete():
    out = window_features("chr1", 0, 250, (0, 10), [_layer("gene", _ensembl())])
    assert len(out["g"]) == 1
    assert out["g"][0][3] == 1000 and len(out["g"][0][5]) == 6
    assert out["g"][0][4] == 0


def test_gene_without_overlap_dropped():
    out = window_features("chr1", 1500, 1600, (1500, 1600), [_layer("gene", _ensembl())])
    assert out["g"] == []


def test_bed_gene_layer_one_block_per_row():
    feats = [Feature(50, 80, "+", "feature", "feature", "", "", "b"),
             Feature(10, 30, "-", "feature", "feature", "", "", "a")]
    out = window_features("chr1", 0, 200, (60, 70), [_layer("gene", feats)])
    assert out["g"] == [["a", "-", 10, 30, 0, [["E", 10, 30]]],
                        ["b", "+", 50, 80, 1, [["E", 50, 80]]]]


def test_te_family_and_half_open_hl():
    feats = [Feature(0, 100, "+", "feature", "LTR/Gypsy", "", "", "Gyp-1"),
             Feature(100, 120, "-", "feature", "DNA/hAT", "", "", "hat"),
             Feature(150, 160, "+", "feature", "DNA/X", "", "", "x"),
             Feature(200, 300, "+", "feature", "DNA/Y", "", "", "y")]
    out = window_features("chr1", 0, 400, (100, 150), [_layer("te", feats)])
    byl = {r[0]: r for r in out["t"]}
    assert byl["Gyp-1"] == ["Gyp-1", "LTR/Gypsy", "+", 0, 100, 0]
    assert byl["hat"][5] == 1
    assert byl["x"][5] == 0  # starts exactly at locus end
    assert byl["y"][5] == 0
    assert out["ly"] == ["te"]


def test_ncrna_and_phas_use_ftype():
    nc = [Feature(10, 20, "+", "tRNA", "tRNA", "", "", "trn")]
    ph = [Feature(30, 40, "-", "PHAS", "PHAS", "", "", "ph1")]
    out = window_features("chr1", 0, 100, (15, 35), [_layer("ncrna", nc), _layer("phas", ph)])
    assert out["r"] == [["trn", "tRNA", "+", 10, 20, 1]]
    assert out["p"] == [["ph1", "PHAS", "-", 30, 40, 1]]
    assert out["ly"] == ["ncrna", "phas"]


def test_long_feature_keeps_real_coordinates():
    feats = [Feature(0, 5000, "+", "feature", "SINE", "", "", "big")]
    out = window_features("chr1", 1000, 2000, (1200, 1300), [_layer("te", feats)])
    assert out["t"][0][3:5] == [0, 5000]


def test_cap_highlighted_first_and_mo():
    feats = [Feature(i * 10, i * 10 + 5, "+", "feature", "F", "", "", f"t{i}")
             for i in range(100)]
    locus = (950, 995)  # overlaps the last features only
    out = window_features("chr1", 0, 1000, locus, [_layer("te", feats)])
    assert len(out["t"]) == FEATURE_CAP == 60
    assert out["mo"]["t"] == 40
    assert sum(r[5] for r in out["t"]) == 5
    starts = [r[3] for r in out["t"]]
    assert starts == sorted(starts)
    assert {"t95", "t99"} <= {r[0] for r in out["t"]}


def test_empty_contig_in_loaded_layer():
    out = window_features("chrZ", 0, 100, (10, 20),
                          [_layer("te", [Feature(0, 5, "+", "feature", "A")]),
                           _layer("gene", [Feature(0, 5, "+", "feature", "A")])])
    assert out["g"] == [] and out["t"] == [] and out["r"] == [] and out["p"] == []
    assert out["mo"] == {"g": 0, "t": 0, "r": 0, "p": 0}
    assert out["ly"] == ["gene", "te"]


def test_build_windows_keys_and_clamp():
    C = namedtuple("C", "name chrom start end")
    feats = [Feature(0, 5, "+", "feature", "A", "", "", "a")]
    out = build_windows([C("cl1", "chr1", 10, 20), C("cl2", "chr1", 480, 495)],
                        [_layer("te", feats)], 50, {"chr1": 500})
    assert set(out) == {"cl1", "cl2"}
    assert out["cl1"]["w"] == [0, 70] and out["cl1"]["L"] == [10, 20]
    assert out["cl2"]["w"] == [430, 500]
    assert len(out["cl1"]["t"]) == 1 and out["cl2"]["t"] == []
