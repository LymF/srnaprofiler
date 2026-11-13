"""Final-review fixes: gene-layer type whitelist, index worst case, parser, validation."""
import time

import pytest

from annotation import (Feature, IntervalIndex, Layer, load_layers,
                        parse_annotation, summarize_locus)

NCBI = (
    "##gff-version 3\n"
    "NC_1\tRefSeq\tregion\t1\t1000\t.\t+\t.\tID=NC_1:1..1000\n"
    "NC_1\tRefSeq\tgene\t101\t500\t.\t+\t.\tID=g1\n"
    "NC_1\tRefSeq\tmRNA\t101\t500\t.\t+\t.\tID=m1;Parent=g1\n"
    "NC_1\tRefSeq\texon\t101\t200\t.\t+\t.\tParent=m1\n"
    "NC_1\tRefSeq\tCDS\t121\t200\t.\t+\t0\tParent=m1\n"
    "NC_1\tRefSeq\texon\t401\t500\t.\t+\t.\tParent=m1\n")
FLYBASE = (
    "##gff-version 3\n"
    "2L\tFlyBase\tgolden_path_region\t1\t1000\t.\t.\t.\tID=2L\n"
    "2L\tFlyBase\tchromosome_arm\t1\t1000\t.\t.\t.\tID=2La\n"
    "2L\tFlyBase\ttransposable_element\t700\t800\t.\t+\t.\tID=te1\n"
    "2L\tFlyBase\tmatch\t700\t800\t.\t+\t.\tID=mt1\n"
    "2L\tFlyBase\tgene\t101\t500\t.\t+\t.\tID=g1\n"
    "2L\tFlyBase\tmRNA\t101\t500\t.\t+\t.\tID=m1\n"
    "2L\tFlyBase\texon\t101\t200\t.\t+\t.\tParent=m1\n"
    "2L\tFlyBase\texon\t401\t500\t.\t+\t.\tParent=m1\n")


def _gene_layer(tmp_path, text, name="g.gff3"):
    p = tmp_path / name
    p.write_text(text)
    chrom = text.split("\n")[1].split("\t")[0]
    return load_layers([f"gene:{p}"], {chrom})[0], chrom


@pytest.mark.parametrize("text", [NCBI, FLYBASE])
def test_stock_gene_file_contexts(tmp_path, text):
    L, c = _gene_layer(tmp_path, text)
    ctx = lambda s, e: summarize_locus(c, s, e, "+", [L], 0.5)["contexts"]
    assert ctx(900, 960) == {"intergenic"}
    assert ctx(250, 350) == {"gene_intron"}
    assert ctx(110, 190) == {"gene_exon"}


def test_gene_layer_type_rules():
    from annotation import _gene_class
    for t in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR", "start_codon",
              "stop_codon", "noncoding_exon", "pseudogenic_exon", "EXON", "feature"):
        assert _gene_class(t) == "exon", t
    for t in ("gene", "pseudogene", "ncRNA_gene", "protein_coding_gene", "mRNA", "lnc_RNA",
              "tRNA", "transcript", "primary_transcript", "pseudogenic_transcript"):
        assert _gene_class(t) == "span", t
    for t in ("region", "chromosome", "biological_region", "golden_path_region",
              "chromosome_arm", "transposable_element", "transposon_fragment", "match",
              "supercontig"):
        assert _gene_class(t) is None, t


def test_ignored_types_not_indexed(tmp_path):
    L, c = _gene_layer(tmp_path, NCBI)
    assert L.n_features == 5  # gene, mRNA, exon, CDS, exon (region ignored)


# ---- I1: index worst case ---------------------------------------------------

def test_index_long_early_and_late_features():
    fs = [Feature(0, 5_000_000, "+", "long", "e")]
    fs += [Feature(i * 100, i * 100 + 50, "+", "s", "") for i in range(1, 300)]
    fs += [Feature(9_000_000, 12_000_000, "+", "long", "l")]
    ix = IntervalIndex(fs)
    got = ix.overlaps(110, 120)
    assert [f.name for f in got] == ["e", ""]
    assert [f.name for f in ix.overlaps(10_000_000, 10_000_001)] == ["l"]
    assert ix.overlaps(6_000_000, 6_000_001) == []
    assert sorted((f.start, f.end) for f in ix.overlaps(0, 20_000_000)) == sorted(
        (f.start, f.end) for f in fs)


def test_index_timing_guard():
    fs = [Feature(i * 50, i * 50 + 30, "+", "s", "") for i in range(200_000)]
    fs.append(Feature(0, 10_000_000, "+", "chr", ""))
    ix = IntervalIndex(fs)
    t = time.time()
    for k in range(2000):
        ix.overlaps(9_000_000 + k * 3, 9_000_000 + k * 3 + 24)
    assert time.time() - t < 1.0


# ---- M4 / M5: parser ---------------------------------------------------------

def test_gff3_percent_decoded_and_id_kept(tmp_path):
    p = tmp_path / "t.gff3"
    p.write_text("c\tEDTA\tLTR_retrotransposon\t1\t10\t.\t+\t.\t"
                 "ID=TE_00001;Name=x;classification=LTR%2FGypsy\n")
    f = parse_annotation(str(p))["c"][0]
    assert f.name == "LTR/Gypsy" and f.fid == "TE_00001"


def test_gtf_id_kept(tmp_path):
    p = tmp_path / "a.gtf"
    p.write_text('chr1\ts\texon\t11\t20\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
    assert parse_annotation(str(p))["chr1"][0].fid == "t1"


def test_parse_known_contigs_counts_unknown_without_keeping(tmp_path):
    p = tmp_path / "a.gff3"
    p.write_text("chr1\ts\texon\t1\t10\t.\t+\t.\t.\n"
                 "scaf9\ts\texon\t1\t10\t.\t+\t.\t.\n"
                 "scaf9\ts\texon\t20\t30\t.\t+\t.\t.\n")
    unknown = {}
    d = parse_annotation(str(p), known_contigs={"chr1"}, unknown=unknown)
    assert list(d) == ["chr1"] and unknown == {"scaf9": 2}
    assert list(parse_annotation(str(p))) == ["chr1", "scaf9"]


def test_text_entry_fifth_field_and_sanitising():
    L = _L("te", [Feature(0, 100, "+", "feature", "LTR/Gypsy", "TE_1")])
    assert summarize_locus("chr1", 10, 30, "+", [L], 0.5)["text"] == \
        "te:LTR/Gypsy:sense:1.00:TE_1"
    L = _L("te", [Feature(0, 100, "+", "feature", "a|b:c", "bad:id")])
    assert summarize_locus("chr1", 10, 30, "+", [L], 0.5)["text"] == "te:a_b_c:sense:1.00"
    L = _L("te", [Feature(0, 100, "+", "feature", "X", "bad|id")])
    assert summarize_locus("chr1", 10, 30, "+", [L], 0.5)["text"] == "te:X:sense:1.00"


def _L(kind, feats, chrom="chr1"):
    return Layer(kind, "x", {chrom: IntervalIndex(feats)}, len(feats))


def test_intergenic_gets_text_entry():
    G = _L("gene", [Feature(0, 100, "+", "gene", "g1")])
    r = summarize_locus("chr1", 200, 220, "+", [G], 0.5)
    assert r["contexts"] == {"intergenic"} and r["text"] == "gene:intergenic:.:1.00"


def test_gff3_omits_empty_annotation_attr(tmp_path):
    from records import ClusterRecord
    from reporter import write_gff3
    r = ClusterRecord(chrom="chr1", start=10, end=30, strand_call="+", cls="piRNA",
                      confidence=0.8, review=False, annotation="",
                      annotation_support="neutral")
    p = tmp_path / "a.gff3"
    write_gff3([r], p, annotated=True)
    attrs = p.read_text().splitlines()[1].split("\t")[8]
    assert "annotation_support=neutral" in attrs and ";annotation=" not in attrs


# ---- I4: config validation ----------------------------------------------------

def _val(cfg):
    from annotation import validate_annotation_config
    return validate_annotation_config(cfg)


GOOD = {"nudge": 0.05, "min_overlap": 0.5,
        "classes": {"piRNA": {"supports": ["te"], "conflicts": ["gene_exon"]}}}


def test_validate_accepts_good_and_shipped():
    assert _val(GOOD)["nudge"] == 0.05
    import yaml
    from pathlib import Path
    for n in ("rules.yaml", "rules.plant.yaml"):
        cfg = yaml.safe_load((Path(__file__).resolve().parents[2] / n).read_text())
        _val(cfg["annotation"])


@pytest.mark.parametrize("bad,msg", [
    ({"nudge": 0.5}, "nudge"), ({"nudge": -0.1}, "nudge"), ({"nudge": "a"}, "nudge"),
    ({"min_overlap": 0}, "min_overlap"), ({"min_overlap": 1.5}, "min_overlap"),
    ({"classes": {"piRNA": {"supports": "te"}}}, "classes.piRNA.supports"),
    ({"classes": {"piRNA": {"conflicts": "gene_exon"}}}, "classes.piRNA.conflicts"),
    ({"classes": {"piRNA": {"supports": ["exon"]}}}, "exon"),
    ({"classes": ["piRNA"]}, "classes"), ({"classes": {"piRNA": ["te"]}}, "piRNA"),
])
def test_validate_rejects(bad, msg):
    with pytest.raises(ValueError, match=msg):
        _val(bad)


def test_validate_null_and_empty_warns_once(capsys):
    out = _val(None)
    assert out.get("classes") in (None, {}) or out["classes"] == {}
    err = capsys.readouterr().err
    assert err.count("[WARN]") == 1 and "no classes" in err
    _val(GOOD)
    assert capsys.readouterr().err.count("[WARN]") == 0


def test_apply_annotation_tolerates_none_cfg():
    from annotation import apply_annotation
    from records import ClusterRecord
    r = ClusterRecord(chrom="chr1", start=10, end=30, strand_call="+", cls="piRNA",
                      confidence=0.8, review=False)
    L = _L("te", [Feature(0, 100, "+", "feature", "X")])
    apply_annotation([r], [L], None)
    assert r.annotation_support == "neutral"


# ---- I3: bad layer files -----------------------------------------------------

def test_truncated_gz_is_clean_error(tmp_path):
    import gzip
    p = tmp_path / "a.gff3.gz"
    blob = gzip.compress(("chr1\ts\texon\t1\t10\t.\t+\t.\t.\n" * 5000).encode())
    p.write_bytes(blob[:len(blob) // 2])
    with pytest.raises(ValueError, match=r"a\.gff3\.gz"):
        load_layers([f"te:{p}"], {"chr1"})


def test_fake_gz_is_clean_error(tmp_path):
    p = tmp_path / "b.gff3.gz"
    p.write_text("chr1\ts\texon\t1\t10\t.\t+\t.\t.\n")
    with pytest.raises(ValueError, match=r"b\.gff3\.gz"):
        load_layers([f"te:{p}"], {"chr1"})
