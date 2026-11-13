"""Final-review fixes: gene name by parent chain, TE trees, _clean."""
import pytest

from annotation import (Feature, IntervalIndex, Layer, load_layers, summarize_locus)

NEST = (
    "##gff-version 3\n"
    "chr1\tt\tgene\t1\t50000\t.\t+\t.\tID=gene:HOST;Name=Host\n"
    "chr1\tt\tmRNA\t1\t50000\t.\t+\t.\tID=transcript:HOST-RA;Parent=gene:HOST\n"
    "chr1\tt\texon\t1\t500\t.\t+\t.\tID=e1;Parent=transcript:HOST-RA\n"
    "chr1\tt\texon\t49500\t50000\t.\t+\t.\tID=e2;Parent=transcript:HOST-RA\n"
    "chr1\tt\tgene\t20001\t22000\t.\t-\t.\tID=gene:NEST;Name=Nested\n"
    "chr1\tt\tmRNA\t20001\t22000\t.\t-\t.\tID=transcript:NEST-RA;Parent=gene:NEST\n"
    "chr1\tt\texon\t20001\t22000\t.\t-\t.\tID=nx1;Parent=transcript:NEST-RA\n"
    "chr1\tt\ttransposable_element_gene\t30001\t32000\t.\t+\t.\tID=gene:TE1;Name=ATTE1\n"
    "chr1\tt\tmRNA\t30001\t32000\t.\t+\t.\tID=transcript:TE1.1;Parent=gene:TE1;Name=ATTE1.1\n"
    "chr1\tt\texon\t30001\t32000\t.\t+\t.\tID=tx1;Parent=transcript:TE1.1\n"
    "chr1\tt\tCDS\t30001\t32000\t.\t+\t0\tID=tc1;Parent=transcript:TE1.1\n")


def _layer(tmp_path, text, name="g.gff3", chrom="chr1"):
    p = tmp_path / name
    p.write_text(text)
    return load_layers([f"gene:{p}"], {chrom})[0]


def test_nested_gene_exon_names_its_own_gene(tmp_path):
    L = _layer(tmp_path, NEST)
    assert summarize_locus("chr1", 20100, 20200, "-", [L], 0.5)["text"] == "gene:Nested:sense:1.00:nx1"


def test_host_exon_names_host(tmp_path):
    L = _layer(tmp_path, NEST)
    assert summarize_locus("chr1", 100, 200, "+", [L], 0.5)["text"] == "gene:Host:sense:1.00:e1"


def test_host_intron_names_host_not_nested(tmp_path):
    L = _layer(tmp_path, NEST)
    r = summarize_locus("chr1", 5000, 5100, "+", [L], 0.5)
    assert r["text"].startswith("gene:Host:sense:1.00")


GTF = (
    'chr1\tt\tgene\t1\t50000\t.\t+\t.\tgene_id "H"; gene_name "Host";\n'
    'chr1\tt\ttranscript\t1\t50000\t.\t+\t.\tgene_id "H"; transcript_id "HT";\n'
    'chr1\tt\texon\t1\t500\t.\t+\t.\tgene_id "H"; transcript_id "HT";\n'
    'chr1\tt\tgene\t20001\t22000\t.\t-\t.\tgene_id "N"; gene_name "Nested";\n'
    'chr1\tt\ttranscript\t20001\t22000\t.\t-\t.\tgene_id "N"; transcript_id "NT";\n'
    'chr1\tt\texon\t20001\t22000\t.\t-\t.\tgene_id "N"; transcript_id "NT";\n')


def test_gtf_chain(tmp_path):
    L = _layer(tmp_path, GTF, "g.gtf")
    r = summarize_locus("chr1", 20100, 20200, "-", [L], 0.5)["text"]
    assert r == "gene:Nested:sense:1.00:NT"


def test_no_parent_chain_falls_back_to_largest_overlap_deterministically():
    fs = [Feature(0, 1000, "+", "gene", "gene", "g1", label="Zed"),
          Feature(0, 1000, "+", "gene", "gene", "g2", label="Alpha"),
          Feature(10, 20, "+", "feature", "feature", "", label="")]
    L = Layer("gene", "x", {"c": IntervalIndex(fs)}, 3)
    t = summarize_locus("c", 12, 18, "+", [L], 0.5)["text"]
    assert t.startswith("gene:Alpha:sense:1.00")


def test_bed_fallback(tmp_path):
    p = tmp_path / "a.bed"
    p.write_text("c\t0\t100\tfoo\t0\t+\n")
    L = load_layers([f"gene:{p}"], {"c"})[0]
    assert summarize_locus("c", 10, 50, "+", [L], 0.5)["text"] == "gene:foo:sense:1.00"


def test_te_gene_children_are_dropped(tmp_path):
    L = _layer(tmp_path, NEST)
    kept = sorted(f.fid for f in L.index["chr1"].overlaps(0, 100000) if f.fid)
    assert not any("TE1" in k or k in ("tx1", "tc1") for k in kept)
    assert L.n_features == 7  # host: 4, nested: 3
    r = summarize_locus("chr1", 30100, 30900, "+", [L], 0.5)
    assert r["contexts"] == {"gene_intron"}  # the host's span only; was gene_exon
    assert "ATTE" not in r["text"]


def test_te_tree_order_independent(tmp_path):
    lines = NEST.split("\n")
    body = [l for l in lines[1:] if l]
    text = lines[0] + "\n" + "\n".join(reversed(body)) + "\n"
    L = _layer(tmp_path, text)
    assert L.n_features == 7


def test_te_children_not_in_windows(tmp_path):
    from collections import namedtuple
    from report_genome import build_windows
    C = namedtuple("C", "name chrom start end")
    L = _layer(tmp_path, NEST)
    out = build_windows([C("cl", "chr1", 30100, 30900)], [L], 100, {"chr1": 60000})
    labels = [r[0] for r in out["cl"]["g"]]
    assert "ATTE1" not in labels and "ATTE1.1" not in labels
    assert labels == ["Host"]  # the TE gene is gone; only the host spans this window


def test_clean_strips_whitespace_controls():
    from annotation import _clean
    assert _clean("a\tb\nc\rd|e:f") == "a b c d_e_f"
    fs = [Feature(0, 100, "+", "feature", "LTR\t/Gy\npsy", "t1")]
    L = Layer("te", "x", {"c": IntervalIndex(fs)}, 1)
    t = summarize_locus("c", 10, 30, "+", [L], 0.5)["text"]
    assert "\t" not in t and "\n" not in t and t.startswith("te:LTR /Gy psy:")
