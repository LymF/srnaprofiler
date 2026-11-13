"""Task 7: predicted ORFs (EXPERIMENTAL): find_orfs, the --predict-orfs flags and the gm "o" key."""
import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
import report_cli  # noqa: E402
import report_orfs  # noqa: E402
import srnaprofiler  # noqa: E402
import test_html_report_cli as base  # noqa: E402
from report_orfs import find_orfs  # noqa: E402
from test_genome_report import GFF, payload  # noqa: E402

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(COMP[c] for c in reversed(s))


def plus_orf(k):
    """A stop, then ATG + k GCT + stop: k + 1 aa. The ORF is the last 3 * (k + 2) bases."""
    return "TAA" + "ATG" + "GCT" * k + "TAA"


def test_plus_orf_exact_coordinates():
    seq = "C" * 50 + plus_orf(20) + "C" * 40
    assert find_orfs(seq, 0, min_aa=10) == [(53, 53 + 3 * 22, "+", (53 % 3))]
    assert find_orfs(seq, 1000, min_aa=10) == [(1053, 1053 + 66, "+", 53 % 3)]


def test_minus_orf_exact_coordinates():
    seq = "C" * 50 + rc(plus_orf(20)) + "C" * 40
    # the leading TTA (stop of the minus ORF) sits at 50; the ORF is 3 * 22 bases from there
    got = find_orfs(seq, 0, min_aa=10)
    assert [(s, e, st) for s, e, st, _ in got] == [(50, 50 + 66, "-")]
    assert find_orfs(seq, 500, min_aa=10)[0][:3] == (550, 616, "-")


def test_frame_is_the_codon_phase_on_the_scanned_strand():
    for pad in (0, 1, 2):
        seq = "C" * (50 + pad) + plus_orf(20) + "C" * 40
        got = find_orfs(seq, 0, min_aa=10)
        assert len(got) == 1 and got[0][3] == (53 + pad) % 3 and got[0][3] in (0, 1, 2)


def test_min_aa_counts_amino_acids_without_the_stop():
    seq = "C" * 30 + plus_orf(20) + "C" * 30     # 21 aa
    assert len(find_orfs(seq, 0, min_aa=21)) == 1
    assert find_orfs(seq, 0, min_aa=22) == []


def test_edge_touching_orfs_are_dropped():
    no_stop = "C" * 30 + "TAA" + "ATG" + "GCT" * 30            # runs off the end
    assert find_orfs(no_stop, 0, min_aa=10) == []
    no_upstream = "ATG" + "GCT" * 30 + "TAA" + "C" * 30         # the start may be further up
    assert find_orfs(no_upstream, 0, min_aa=10) == []
    reverse_no_stop = rc(no_stop)
    assert find_orfs(reverse_no_stop, 0, min_aa=10) == []


def test_n_breaks_an_orf():
    seq = "C" * 30 + "TAA" + "ATG" + "GCT" * 15 + "NNN" + "GCT" * 15 + "TAA" + "C" * 30
    assert find_orfs(seq, 0, min_aa=10) == []
    # the codon after the N is not a stop: the second ORF needs its own stop and ATG
    seq3 = "C" * 30 + "TAA" + "ATG" + "GCT" * 15 + "NNN" + "TAA" + "ATG" + "GCT" * 15 + "TAA" + "C" * 9
    assert [(s, e) for s, e, *_ in find_orfs(seq3, 0, min_aa=10)] == [(87, 87 + 3 * 17)]


def test_lowercase_is_read():
    seq = ("C" * 30 + plus_orf(20) + "C" * 30).lower()
    assert len(find_orfs(seq, 0, min_aa=10)) == 1


def test_nested_orfs_sharing_a_stop_collapse_to_the_longest():
    seq = "C" * 30 + "TAA" + "ATG" + "GCT" * 5 + "ATG" + "GCT" * 20 + "TAA" + "C" * 30
    got = find_orfs(seq, 0, min_aa=10)
    assert len(got) == 1
    s, e, st, _ = got[0]
    assert (s, e, st) == (33, 33 + 3 + 15 + 3 + 60 + 3, "+")      # starts at the first ATG


def test_back_to_back_orfs_in_one_frame_are_both_kept():
    seq = "C" * 30 + "TAA" + "ATG" + "GCT" * 20 + "TAA" + "ATG" + "GCT" * 20 + "TAA" + "C" * 30
    got = find_orfs(seq, 0, min_aa=10)
    assert [(s, e) for s, e, *_ in got] == [(33, 33 + 66), (99, 99 + 66)]


# ---- the payload, end to end ----

def fasta(tmp_path, seq, index=True, name="g.fa"):
    p = tmp_path / name
    p.write_text(">c1\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    if index:
        pysam.faidx(str(p))
    return str(p)


def genome_seq():
    # contig c1 (5000 bp, like the test BAM): plus ORF at 300 (21 aa), minus ORF at 1000
    parts = ["C" * 300 + plus_orf(20)]
    s = "".join(parts)
    s += "C" * (1000 - len(s)) + rc(plus_orf(20))
    return s + "C" * (5000 - len(s))


def gff_cds(tmp_path, a, b):
    p = tmp_path / "cds.gff3"
    p.write_text("##gff-version 3\n"
                 f"c1\ts\tgene\t{a}\t{b}\t.\t+\t.\tID=gene:G1;Name=G1\n"
                 f"c1\ts\tmRNA\t{a}\t{b}\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
                 f"c1\ts\texon\t{a}\t{b}\t.\t+\t.\tID=e1;Parent=transcript:T1\n"
                 f"c1\ts\tCDS\t{a}\t{b}\t.\t+\t0\tID=c1;Parent=transcript:T1\n")
    return str(p)


def run(tmp_path, ann, **kw):
    a = base.args_for(tmp_path, annotation=["gene:" + ann], flank=2000, predict_orfs=True,
                      genome=fasta(tmp_path, genome_seq()), orf_min_aa=kw.pop("orf_min_aa", 20), **kw)
    report_cli.run_report(a)
    p = payload(tmp_path)
    idx = [i for i, r in enumerate(p["rows"]) if r["n"] == "cl1"][0]
    return p, p["details"][str(idx)]["gm"]


def test_orfs_in_gm_and_meta(tmp_path):
    g = tmp_path / "far.gff3"
    g.write_text(GFF)      # one gene at 60-400 without CDS lines
    p, gm = run(tmp_path, str(g))
    assert len(gm["o"]) == 2
    assert p["meta"]["orf_min_aa"] == 20 and p["meta"]["gm"] == 1


def test_orf_coordinates_in_the_payload(tmp_path):
    g = tmp_path / "far.gff3"
    g.write_text(GFF)
    _, gm = run(tmp_path, str(g))
    # plus_orf: leading stop at 300..303, ATG at 303; the ORF (ATG..stop) is 66 bases
    assert gm["o"] == [[303, 369, "+", 21], [1000, 1066, "-", 21]]


def test_orf_overlapping_a_cds_is_excluded(tmp_path):
    _, gm = run(tmp_path, gff_cds(tmp_path, 350, 500))     # 1-based 350..500 overlaps the plus ORF (303..369)
    assert gm["o"] == [[1000, 1066, "-", 21]]
    _, gm = run(tmp_path, gff_cds(tmp_path, 370, 500))     # 0-based 369..500: just past the ORF (half-open)
    assert [r[2] for r in gm["o"]] == ["+", "-"]


def test_orfs_without_overlap_when_everything_is_coding_leaves_no_o(tmp_path):
    _, gm = run(tmp_path, gff_cds(tmp_path, 250, 1100))
    assert "o" not in gm


def test_min_aa_filters_in_the_payload(tmp_path):
    _, gm = run(tmp_path, gff_cds(tmp_path, 9000, 9100), orf_min_aa=22)
    assert "o" not in gm


def test_orf_cap_longest_first_and_mo(tmp_path):
    k = [10 + i for i in range(70)]               # 70 ORFs of 11..80 aa
    seq = "".join("C" * 6 + plus_orf(x) for x in k)
    p = tmp_path / "big.fa"
    p.write_text(">c1\n" + seq + "\n")
    pysam.faidx(str(p))
    from types import SimpleNamespace as NS
    c = NS(name="x", chrom="c1", start=0, end=10)
    gmd = {"x": {"w": [0, len(seq)], "L": [0, 10], "mo": {}}}
    report_orfs.add_orfs(gmd, [c], [], str(p), 11)
    g = gmd["x"]
    assert len(g["o"]) == 60 and g["mo"]["o"] == 10
    assert min(r[3] for r in g["o"]) == 21           # the ten shortest (11..20 aa) dropped
    assert [r[0] for r in g["o"]] == sorted(r[0] for r in g["o"])


def test_contig_missing_from_the_genome_gives_no_orfs(tmp_path):
    from types import SimpleNamespace as NS
    p = fasta(tmp_path, genome_seq())
    gmd = {"x": {"w": [0, 100], "L": [0, 10], "mo": {}}}
    report_orfs.add_orfs(gmd, [NS(name="x", chrom="zzz", start=0, end=10)], [], p, 10)
    assert "o" not in gmd["x"]


# ---- flags ----

def test_parser_flags():
    a = srnaprofiler.build_parser().parse_args(["report", "--clusters", "a", "--bam", "b", "-o", "o"])
    assert a.predict_orfs is False and a.genome is None and a.orf_min_aa == 100
    a = srnaprofiler.build_parser().parse_args(
        ["report", "--clusters", "a", "--bam", "b", "-o", "o", "--predict-orfs", "--genome", "g.fa",
         "--orf-min-aa", "50"])
    assert a.predict_orfs is True and a.genome == "g.fa" and a.orf_min_aa == 50


def test_help_says_experimental():
    rp = [a for a in srnaprofiler.build_parser()._subparsers._group_actions[0].choices["report"]._actions
          if "--predict-orfs" in a.option_strings][0]
    assert "EXPERIMENTAL" in rp.help


def _spy(monkeypatch):
    calls = []
    monkeypatch.setattr(report_cli.report_signals, "signals_for", lambda *a, **k: calls.append(1) or {})
    return calls


def test_predict_orfs_needs_genome(tmp_path, monkeypatch):
    calls = _spy(monkeypatch)
    g = tmp_path / "g.gff3"
    g.write_text(GFF)
    a = base.args_for(tmp_path, annotation=["gene:" + str(g)], predict_orfs=True, genome=None)
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(a)
    assert "--genome" in str(e.value) and not calls


def test_predict_orfs_needs_annotation(tmp_path, monkeypatch):
    calls = _spy(monkeypatch)
    a = base.args_for(tmp_path, annotation=None, predict_orfs=True, genome=fasta(tmp_path, genome_seq()))
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(a)
    assert "--annotation" in str(e.value) and not calls


def test_unindexed_or_missing_genome(tmp_path, monkeypatch):
    calls = _spy(monkeypatch)
    g = tmp_path / "g.gff3"
    g.write_text(GFF)
    ann = ["gene:" + str(g)]
    a = base.args_for(tmp_path, annotation=ann, predict_orfs=True, genome=fasta(tmp_path, genome_seq(), index=False))
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(a)
    assert ".fai" in str(e.value) and "samtools faidx" in str(e.value)
    a = base.args_for(tmp_path, annotation=ann, predict_orfs=True, genome=str(tmp_path / "nope.fa"))
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(a)
    assert "not found" in str(e.value)
    assert not calls


def test_min_aa_must_be_positive(tmp_path, monkeypatch):
    calls = _spy(monkeypatch)
    g = tmp_path / "g.gff3"
    g.write_text(GFF)
    a = base.args_for(tmp_path, annotation=["gene:" + str(g)], predict_orfs=True,
                      genome=fasta(tmp_path, genome_seq()), orf_min_aa=0)
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(a)
    assert "--orf-min-aa" in str(e.value) and not calls


def test_default_report_has_no_orf_trace(tmp_path):
    g = tmp_path / "g.gff3"
    g.write_text(GFF)
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + str(g)], flank=2000))
    p = payload(tmp_path)
    assert "orf_min_aa" not in p["meta"]
    assert all("o" not in d.get("gm", {}) for d in p["details"].values())
    assert "orf" not in str(p["details"]).lower().replace("orfs", "")


def test_genome_without_the_flag_is_inert(tmp_path):
    g = tmp_path / "g.gff3"
    g.write_text(GFF)
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + str(g)], flank=2000))
    one = payload(tmp_path)
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + str(g)], flank=2000,
                                        predict_orfs=False, genome=fasta(tmp_path, genome_seq()), orf_min_aa=100))
    assert payload(tmp_path) == one


def test_contract_documents_o():
    import report_html
    doc = report_html.__doc__
    assert "orf_min_aa" in doc and " o " in doc
