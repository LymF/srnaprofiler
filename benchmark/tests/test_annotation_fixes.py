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
