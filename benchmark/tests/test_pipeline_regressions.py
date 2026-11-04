"""Regressions found by the first real benchmark run (SRR30544370)."""
import subprocess
from pathlib import Path

import mirna
import qc


def test_filter_fasta_min_length_drops_short_reads(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(
        ">a\n" + "A" * 15 + "\n"
        ">b\n" + "C" * 17 + "\n"
        ">c\n" + "G" * 18 + "\n"
        ">d\n" + "T" * 22 + "\n"
    )
    out = tmp_path / "out.fa"
    kept, dropped = mirna.filter_fasta_min_length(src, out, 18)
    assert (kept, dropped) == (2, 2)
    assert out.read_text() == ">c\n" + "G" * 18 + "\n>d\n" + "T" * 22 + "\n"


def test_filter_fasta_min_length_joins_wrapped_sequences(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(">w\n" + "A" * 10 + "\n" + "C" * 10 + "\n>s\nACGT\n")
    out = tmp_path / "out.fa"
    kept, dropped = mirna.filter_fasta_min_length(src, out, 18)
    assert (kept, dropped) == (1, 1)
    assert out.read_text() == ">w\n" + "A" * 10 + "C" * 10 + "\n"


def test_fastqc_uses_outdir_flag(monkeypatch, tmp_path):
    seen = []
    monkeypatch.setattr(qc, "which", lambda name: "/bin/fastqc")
    monkeypatch.setattr(qc, "run_command", lambda cmd, *a, **k: seen.append(cmd))
    qc.fastqc("reads.fq", tmp_path, 2)
    cmd = seen[0]
    assert "--outdir" in cmd and cmd[cmd.index("--outdir") + 1] == str(tmp_path)
    assert "-o" not in cmd


# --- adapter detection, compound read extensions, quantifier without references ---

import gzip

import pytest


def _write_fastq(path, seqs):
    text = "".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(seqs))
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        Path(path).write_text(text)


BIO = "ACGTTGCAGGCTAACCGTTAGC"


def test_detect_adapter_truseq(tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "AGATCGGAAGAGCACAC"] * 10)
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS) == "AGATCGGAAGAGC"


def test_detect_adapter_illumina_srna(tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "TGGAATTCTCGGGTGCC"] * 10)
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS) == qc.FASTQ_ADAPTERS["illumina_srna"]


def test_detect_adapter_none_when_absent(tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "CCCCCCCCCCCC"] * 10)
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS) is None


def test_detect_adapter_plain_file(tmp_path):
    fq = tmp_path / "t.fastq"
    _write_fastq(fq, [BIO + "AGATCGGAAGAGCAC"] * 10)
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS) == "AGATCGGAAGAGC"


def test_detect_adapter_respects_read_cap(tmp_path):
    fq = tmp_path / "t.fastq"
    _write_fastq(fq, [BIO + "CCCCCCCCCCCC"] * 20 + [BIO + "AGATCGGAAGAGCAC"] * 20)
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS, n_reads=20) is None
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS, n_reads=40) == "AGATCGGAAGAGC"


def _run_fastp(monkeypatch, tmp_path, reads, adapter, in2=None):
    seen = []
    monkeypatch.setattr(qc, "which", lambda name: "/bin/fastp")
    monkeypatch.setattr(qc, "run_command", lambda cmd, *a, **k: seen.append(cmd))
    out2 = tmp_path / "o2.fq" if in2 else None
    qc.trim_with_fastp(str(reads), tmp_path / "o1.fq", in2, out2, 18, 20,
                       adapter, 1, 33, "lib")
    return seen[0]


def test_fastp_detects_truseq_when_no_adapter(monkeypatch, tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "AGATCGGAAGAGCAC"] * 10)
    cmd = _run_fastp(monkeypatch, tmp_path, fq, None)
    assert cmd[cmd.index("--adapter_sequence") + 1] == "AGATCGGAAGAGC"


def test_fastp_single_end_without_adapter_passes_none(monkeypatch, tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "CCCCCCCCCCCC"] * 10)
    cmd = _run_fastp(monkeypatch, tmp_path, fq, None)
    assert "--adapter_sequence" not in cmd


def test_fastp_paired_end_keeps_detect_flag(monkeypatch, tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "CCCCCCCCCCCC"] * 10)
    cmd = _run_fastp(monkeypatch, tmp_path, fq, None, in2=str(fq))
    assert "--detect_adapter_for_pe" in cmd and "--adapter_sequence" not in cmd


def test_fastp_explicit_adapter_skips_detection(monkeypatch, tmp_path):
    def boom(*a, **k):
        raise AssertionError("detect_adapter must not be called")
    monkeypatch.setattr(qc, "detect_adapter", boom)
    cmd = _run_fastp(monkeypatch, tmp_path, tmp_path / "x.fq", "ACGTACGT")
    assert cmd[cmd.index("--adapter_sequence") + 1] == "ACGTACGT"


@pytest.mark.parametrize("name,expected", [
    ("x.input_reads.fastq.gz", "x_input_reads.fastq.gz"),
    ("x.ncrna_filtered.fastq", "x_ncrna_filtered.fastq"),
    ("reads.fq.gz", "reads.fq.gz"),
    ("reads", "reads"),
    ("a.b.fa.bz2", "a_b.fa.bz2"),
])
def test_staged_reads_name_keeps_compound_extension(name, expected):
    assert mirna.staged_reads_name(name) == expected


def test_quantifier_skipped_without_references(monkeypatch, tmp_path):
    def boom(*a, **k):
        raise AssertionError("quantifier.pl must not run")
    monkeypatch.setattr(mirna.subprocess, "run", boom)
    assert mirna.run_quantifier(tmp_path / "r.fa", tmp_path, "Drosophila melanogaster") is True


def _mdp2_file(tmp_path, rows):
    """A miRDeep-P2 `*_filter_P_prediction` file: (read id, name, mature coords, precursor coords, mature)."""
    path = tmp_path / "filter_P_prediction"
    lines = [f"1\t-\t{read}\t{name}\t{mature_at}\t{pre_at}\t{seq}\t....."
             for read, name, mature_at, pre_at, seq in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_mirdeep_p2_repeated_precursor_is_renamed_not_dropped(tmp_path):
    path = _mdp2_file(tmp_path, [
        ("read02059488_x61", "1_139", "32686180..32686200", "32686095..32686200", "TGAAGCTGCCAGCATGATCTG"),
        ("read02059489_x666", "1_139", "32686179..32686200", "32686095..32686200", "TGAAGCTGCCAGCATGATCTGA"),
        ("read00000001_x7", "1_200", "100..120", "90..130", "GGGGCCCCAAAATTTTGGGG"),
    ])
    preds = mirna.parse_mirdeep_p2_predictions(path)
    assert [p["name"] for p in preds] == ["1_139", "1_139_iso2", "1_200"]
    assert preds[0]["p2_name"] == "1_139" and preds[1]["p2_name"] == "1_139"
    assert preds[1]["mature"] == "TGAAGCTGCCAGCATGATCTGA"
    assert preds[0]["start"] == 32686095


def test_mirdeep_p2_three_duplicates_get_iso2_and_iso3(tmp_path):
    path = _mdp2_file(tmp_path, [
        ("r1_x1", "4_1", "10..30", "1..60", "AAAAAAAAAAAAAAAAAAAAA"),
        ("r2_x1", "4_1", "11..31", "1..60", "CCCCCCCCCCCCCCCCCCCCC"),
        ("r3_x1", "4_1", "12..32", "1..60", "GGGGGGGGGGGGGGGGGGGGG"),
    ])
    preds = mirna.parse_mirdeep_p2_predictions(path)
    assert [p["name"] for p in preds] == ["4_1", "4_1_iso2", "4_1_iso3"]
    assert {p["p2_name"] for p in preds} == {"4_1"}


def test_mirdeep_p2_suffix_already_in_file_is_not_reused(tmp_path):
    path = _mdp2_file(tmp_path, [
        ("a_x1", "5_1", "10..30", "1..60", "AAAAAAAAAAAAAAAAAAAAA"),
        ("b_x1", "5_1_iso2", "10..30", "1..60", "CCCCCCCCCCCCCCCCCCCCC"),
        ("c_x1", "5_1", "10..30", "1..60", "GGGGGGGGGGGGGGGGGGGGG"),
    ])
    preds = mirna.parse_mirdeep_p2_predictions(path)
    assert [p["name"] for p in preds] == ["5_1", "5_1_iso2", "5_1_iso3"]
    assert [p["p2_name"] for p in preds] == ["5_1", "5_1_iso2", "5_1"]


def test_validate_writes_unique_names_to_mature_and_hairpin(monkeypatch, tmp_path):
    path = _mdp2_file(tmp_path, [
        ("r1_x1", "1_139", "10..30", "1..60", "AAAAAAAAAAAAAAAAAAAAA"),
        ("r2_x9", "1_139", "11..31", "1..60", "CCCCCCCCCCCCCCCCCCCCC"),
        ("r3_x1", "1_200", "10..30", "1..60", "GGGGGGGGGGGGGGGGGGGGG"),
    ])
    preds = mirna.parse_mirdeep_p2_predictions(path)
    workdir = tmp_path / "ms"

    def fake_run(cmd, **kwargs):
        out = workdir / "mirscore"
        out.mkdir(parents=True, exist_ok=True)
        (out / "miRScore_results.csv").write_text("")

    monkeypatch.setattr(mirna.shutil, "which", lambda name: "/bin/miRScore")
    monkeypatch.setattr(mirna, "open_genome", lambda path: object())
    monkeypatch.setattr(mirna, "hairpin_with_flank",
                        lambda genome, chrom, start, end, mature, flank=25: "ACGT" * 10)
    monkeypatch.setattr(mirna, "stage_reads_for_mirscore", lambda reads, workdir: reads)
    monkeypatch.setattr(mirna.subprocess, "run", fake_run)
    monkeypatch.setattr(mirna, "read_mirscore_rows", lambda p: {})
    monkeypatch.setattr(mirna, "no_reads_counted", lambda p: False)

    result = mirna.validate_with_mirscore(preds, "genome.fa", tmp_path / "reads.fa",
                                          "Plant", workdir)
    expected = ["1_139", "1_139_iso2", "1_200"]
    for fa in ("mature.fa", "hairpin.fa"):
        headers = [line[1:] for line in (workdir / fa).read_text().splitlines()
                   if line.startswith(">")]
        assert headers == expected
    assert result.submitted == set(expected)


# --- randomized insert ends (NEXTflex-style 4N) ---

import argparse

import mapping


@pytest.mark.parametrize("rate0,rate4,expected", [
    (0.0014, 0.84, 4),
    (0.93, 0.95, 0),
    (0.10, 0.25, 0),
    (0.0, 0.0, 0),
])
def test_choose_random_ends(rate0, rate4, expected):
    assert qc.choose_random_ends(rate0, rate4) == expected


def test_trim_random_ends_removes_ends_and_drops_short_reads(tmp_path):
    fq = tmp_path / "t.fastq.gz"
    long_seq = "AAAA" + "ACGTACGTACGTACGTAC" + "TTTT"   # 26 nt -> 18 nt
    _write_fastq(fq, [long_seq, "AAAA" + "ACGTACG" + "TTTT", long_seq])  # middle: 15 nt < 18
    kept, dropped = qc.trim_random_ends([fq], 18)
    assert (kept, dropped) == (2, 1)
    with gzip.open(fq, "rt") as fh:
        lines = fh.read().splitlines()
    assert lines[1] == "ACGTACGTACGTACGTAC" and lines[3] == "I" * 18
    assert len(lines) == 8
    assert not (tmp_path / "tmp.t.fastq.gz").exists()


def test_trim_random_ends_keeps_mates_in_sync(tmp_path):
    r1, r2 = tmp_path / "a_1.fq", tmp_path / "a_2.fq"
    _write_fastq(r1, ["N" * 4 + "A" * 20 + "N" * 4, "N" * 4 + "C" * 20 + "N" * 4])
    _write_fastq(r2, ["N" * 4 + "G" * 20 + "N" * 4, "N" * 4 + "T" * 5 + "N" * 4])
    assert qc.trim_random_ends([r1, r2], 18) == (1, 1)
    assert r1.read_text().splitlines()[1] == "A" * 20
    assert r2.read_text().splitlines()[1] == "G" * 20


def _qc_args(tmp_path, **kw):
    d = dict(random_ends="auto", index_prefix=None, prefix="lib", threads=2, min_len=18)
    d.update(kw)
    return argparse.Namespace(**d)


def _no_bowtie(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("bowtie must not be called")
    monkeypatch.setattr(qc, "measure_random_end_rates", boom)
    monkeypatch.setattr(mapping, "bowtie_build", boom)


def test_random_ends_zero_never_detects(monkeypatch, tmp_path):
    _no_bowtie(monkeypatch)
    fq = tmp_path / "t.fastq"
    _write_fastq(fq, [BIO])
    stats = qc.resolve_random_ends(_qc_args(tmp_path, random_ends="0", index_prefix="x"),
                                   [fq], tmp_path / "g.fa", tmp_path)
    assert stats["random_ends"] == 0 and fq.read_text().splitlines()[1] == BIO


def test_random_ends_auto_without_genome_or_index_is_zero(monkeypatch, tmp_path):
    _no_bowtie(monkeypatch)
    fq = tmp_path / "t.fastq"
    _write_fastq(fq, [BIO])
    stats = qc.resolve_random_ends(_qc_args(tmp_path), [fq], None, tmp_path)
    assert stats["random_ends"] == 0 and stats["random_ends_rate0"] is None


def test_random_ends_auto_chooses_4_and_trims(monkeypatch, tmp_path):
    fq = tmp_path / "t.fastq"
    _write_fastq(fq, ["ACGT" + BIO + "TGCA"])
    seen = {}

    def fake_rates(trimmed, index, threads):
        seen["index"] = index
        return 0.001, 0.84
    monkeypatch.setattr(qc, "measure_random_end_rates", fake_rates)
    stats = qc.resolve_random_ends(_qc_args(tmp_path, index_prefix="/idx/sp"),
                                   [fq], None, tmp_path)
    assert seen["index"] == "/idx/sp"
    assert stats["random_ends"] == 4 and stats["random_ends_rate4"] == 0.84
    assert (stats["random_ends_kept"], stats["random_ends_dropped"]) == (1, 0)
    assert fq.read_text().splitlines()[1] == BIO


def test_random_ends_auto_builds_index_where_mapping_expects_it(monkeypatch, tmp_path):
    built = []
    monkeypatch.setattr(mapping, "bowtie_build", lambda g, p, t=1: built.append((g, p)))
    monkeypatch.setattr(qc, "measure_random_end_rates", lambda *a: (0.9, 0.9))
    fq = tmp_path / "t.fastq"
    _write_fastq(fq, [BIO])
    genome = tmp_path / "lib.genome.std.fa"
    stats = qc.resolve_random_ends(_qc_args(tmp_path), [fq], genome, tmp_path)
    assert built == [(str(genome), str(tmp_path / "lib.bowtie_index"))]
    assert stats["random_ends"] == 0


def _six_index_files(prefix, content="x"):
    for n in ("1", "2", "3", "4", "rev.1", "rev.2"):
        Path(f"{prefix}.{n}.ebwt").write_text(content)


def _mapping_args(tmp_path, genome):
    return argparse.Namespace(
        outdir=str(tmp_path), prefix="lib", genome=str(genome), index_prefix=None,
        reads=str(tmp_path / "r.fq"), reads_paired=None, sam_out=None, fmt=None, mismatches=1,
        k=1, keep_unmapped=False, threads=1)


class _StopAfterIndex(Exception):
    pass


def _run_mapping_until_index(monkeypatch, tmp_path):
    genome = tmp_path / "g.fa"
    genome.write_text(">c\nACGT\n")
    calls = []
    monkeypatch.setattr(mapping, "validate_arguments", lambda a: None)
    monkeypatch.setattr(mapping, "check_disk_space", lambda *a: True)
    monkeypatch.setattr(mapping, "bowtie_build", lambda *a: calls.append(a))
    monkeypatch.setattr(mapping, "infer_format",
                        lambda *a: (_ for _ in ()).throw(_StopAfterIndex()))
    monkeypatch.setattr(mapping, "convert_uracil_to_thymine_optimized", lambda *a: False)
    monkeypatch.setattr(mapping, "get_read_stats", lambda *a: {})
    (tmp_path / "r.fq").write_text("@r\nACGT\n+\nIIII\n")
    with pytest.raises(_StopAfterIndex):
        mapping.run_mapping(_mapping_args(tmp_path, genome))
    return calls


def test_mapping_reuses_complete_index(monkeypatch, tmp_path):
    _six_index_files(tmp_path / "lib.bowtie_index")
    assert _run_mapping_until_index(monkeypatch, tmp_path) == []


def test_mapping_rebuilds_incomplete_index(monkeypatch, tmp_path):
    _six_index_files(tmp_path / "lib.bowtie_index")
    (tmp_path / "lib.bowtie_index.rev.2.ebwt").write_text("")
    assert len(_run_mapping_until_index(monkeypatch, tmp_path)) == 1


def test_run_passes_shared_index_to_qc(monkeypatch, tmp_path):
    from benchmark import run as run_mod
    from benchmark.manifest import Library
    argvs = []

    def fake_run(cmd, **kw):
        if cmd[0] == "git":
            raise subprocess.CalledProcessError(1, cmd)
        argvs.append(cmd)
        raise StopIteration  # stop after qc
    monkeypatch.setattr(run_mod.subprocess, "run", fake_run)
    for name, data in (("r.fq.gz", b"r"), ("g.fa", b">c\nA\n"), ("rules.yaml", b"a: 1\n")):
        (tmp_path / name).write_bytes(data)
    lib = Library("SRR9", "Bombyx mori", "test", "arthropod", "ovary", 1, "http://x", ["t"])
    shared = tmp_path / "idx" / "Bm"
    shared.parent.mkdir()
    for kw in ({}, {"index_prefix": shared}):
        argvs.clear()
        with pytest.raises(StopIteration):
            run_mod.run_library(lib, tmp_path / "r.fq.gz", tmp_path / "g.fa",
                                tmp_path / "w", tmp_path / "rules.yaml", **kw)
        assert "--index-prefix" not in argvs[0]      # shared index not built yet
    _six_index_files(shared)
    argvs.clear()
    with pytest.raises(StopIteration):
        run_mod.run_library(lib, tmp_path / "r.fq.gz", tmp_path / "g.fa",
                            tmp_path / "w", tmp_path / "rules.yaml", index_prefix=shared)
    qc_cmd = argvs[0]
    assert qc_cmd[2] == "qc" and qc_cmd[qc_cmd.index("--index-prefix") + 1] == str(shared)
    assert qc_cmd[qc_cmd.index("-g") + 1] == str(tmp_path / "g.fa")
