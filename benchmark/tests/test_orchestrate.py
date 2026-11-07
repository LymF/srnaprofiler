import json
from pathlib import Path

import pysam
import pytest
from bamutil import make_bam

from benchmark import orchestrate
from benchmark.manifest import Library
from benchmark.run import RunResult

CLUSTER_HEADER = "chrom\tstart\tend\treads\tcls\tconfidence\n"


def setup(tmp_path, monkeypatch, *, truth_chrom="chr1", index_rows=None, role="test",
          truth_sources=None, extra_files=None):
    repo = tmp_path / "repo"
    (repo / "benchmark" / "truth").mkdir(parents=True)
    (repo / "rules.yaml").write_text("a: 1\n")
    (repo / "rules.plant.yaml").write_text("a: 2\n")
    (repo / "pirna.bed").write_text(f"{truth_chrom}\t100\t200\tpic1\n")
    (repo / "mirna.bed").write_text(f"{truth_chrom}\t500\t560\tmir1\n")
    for name, text in (extra_files or {}).items():
        (repo / name).write_text(text)
    rows = ["path\tcls\tsource\tfeature_types\taliases"]
    if index_rows is None:
        index_rows = ["pirna.bed\tpiRNA\tSrcPi\t\t", "mirna.bed\tmiRNA\tSrcMi\t\t"]
    rows += index_rows
    index = repo / "benchmark" / "truth" / "truth_index.tsv"
    index.write_text("\n".join(rows) + "\n")
    config = repo / "config.yaml"
    config.write_text("tiny_bp: [50]\nlow_reads: [20]\nmin_overlap_frac: 0.5\n")
    monkeypatch.setattr(orchestrate, "REPO", repo)
    monkeypatch.setattr(orchestrate, "TRUTH_INDEX", index)
    monkeypatch.setattr(orchestrate, "CONFIG", config)

    sources = ["pirna.bed", "mirna.bed"] if truth_sources is None else truth_sources
    lib = Library("SRR9", "Bombyx mori", role, "arthropod", "ovary", 1, "http://x", sources)
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "accession\torganism\trole\ttaxon_group\ttissue\treplicate\tgenome_url\ttruth_sources\n"
        f"SRR9\tBombyx mori\t{role}\tarthropod\tovary\t1\thttp://x\t{','.join(sources)}\n")

    seen = {}

    def fake_genome(l, dest):
        seen["genome_dest"] = Path(dest)
        raw = tmp_path / "raw.fa"
        raw.write_text(">chr1 Homo sapiens chromosome 1\nACGT\n>chr2\nACGT\n")
        return raw

    def fake_reads(l, dest):
        seen["reads_dest"] = Path(dest)
        seen["order"] = seen.get("order", []) + ["reads"]
        return tmp_path / "r.fq.gz"

    def fake_run(l, reads, genome, workdir, rules, threads=8, index_prefix=None, classifier="rules", ml_model=None):
        seen["index_prefix"] = index_prefix
        seen["order"] = seen.get("order", []) + ["run"]
        seen["rules"] = Path(rules)
        seen["threads"] = threads
        workdir.mkdir(parents=True, exist_ok=True)
        (workdir / "SRR9.genome.std.fa").write_text(">chr1\n" + "A" * 1000 + "\n")
        # 12 alignments on the piRNA truth region (expressed), none on the miRNA region
        make_bam(workdir / "SRR9.bowtie.bam", [("chr1", 100 + i, False) for i in range(12)])
        clusters = workdir / "SRR9.clusters.tsv"
        clusters.write_text(CLUSTER_HEADER +
                            "chr1\t100\t200\t50\tpiRNA\t0.9\n"
                            "chr1\t500\t560\t30\tmiRNA\t0.8\n")
        meta = workdir / "benchmark_run.json"
        meta.write_text(json.dumps({"role": role, "rules_sha256": "x"}))
        return RunResult(clusters, workdir / "b.bam", meta)

    monkeypatch.setattr(orchestrate, "fetch_genome", fake_genome)
    monkeypatch.setattr(orchestrate, "fetch_library", fake_reads)
    monkeypatch.setattr(orchestrate, "run_library", fake_run)
    return manifest, seen


def go(tmp_path, manifest, acc="SRR9"):
    return orchestrate.run_benchmark_library(
        acc, manifest, tmp_path / "data", tmp_path / "results", 4)


def test_run_writes_metrics_json(tmp_path, monkeypatch):
    manifest, seen = setup(tmp_path, monkeypatch)
    path = go(tmp_path, manifest)
    assert path == tmp_path / "results" / "SRR9" / "metrics.json"
    m = json.loads(path.read_text())
    assert set(m) == {"accession", "organism", "role", "classification", "delimitation",
                      "by_source", "n_loci"}
    assert set(m["delimitation"]) == {"global", "by_source"}
    assert set(m["delimitation"]["global"]) == {"tiny_counts", "low_read_counts",
                                                   "n_fragments", "fragment_read_frac"}
    assert m["delimitation"]["global"]["tiny_counts"] == {"50": 0}
    assert m["n_loci"] == 2 and m["role"] == "test"
    assert set(m["by_source"]) == {"SrcPi", "SrcMi"}
    # delimitation only sees the piRNA cluster, not the miRNA precursor
    assert set(m["delimitation"]["by_source"]) == {"SrcPi"}
    assert m["delimitation"]["by_source"]["SrcPi"]["summary"]["n_clusters"] == 1
    assert "miRNA" in m["classification"]["per_class"]
    pi = m["classification"]["per_class"]["piRNA"]
    assert (pi["n_truth"], pi["n_truth_expressed"], pi["recall"], pi["recall_expressed"]) == (1, 1, 1.0, 1.0)
    mi = m["classification"]["per_class"]["miRNA"]
    assert (mi["n_truth"], mi["n_truth_expressed"], mi["recall"], mi["recall_expressed"]) == (1, 0, 1.0, None)
    assert m["by_source"]["SrcPi"]["piRNA"]["n_truth_expressed"] == 1
    assert (path.parent / "benchmark_run.json").exists()
    assert seen["rules"] == orchestrate.REPO / "rules.yaml"
    assert seen["threads"] == 4
    # one bowtie index per species, shared by its libraries
    assert seen["index_prefix"] == tmp_path / "results" / "_index" / "Bombyx_mori" / "Bombyx_mori"
    assert seen["genome_dest"] == tmp_path / "data" / "genomes" / "Bombyx_mori"
    assert seen["reads_dest"] == tmp_path / "data" / "reads"


def test_delimitation_is_per_source_and_independent_of_other_sources(tmp_path, monkeypatch):
    base_manifest, _ = setup(tmp_path / "a", monkeypatch)
    base = json.loads(go(tmp_path / "a", base_manifest).read_text())
    # a second source overlapping the piRNA cluster must not change SrcPi's numbers
    rows = ["pirna.bed\tpiRNA\tSrcPi\t\t", "mirna.bed\tmiRNA\tSrcMi\t\t",
            "other.bed\tsiRNA\tSrcOther\t\t"]
    manifest, _ = setup(tmp_path / "b", monkeypatch, index_rows=rows,
                        truth_sources=["pirna.bed", "mirna.bed", "other.bed"],
                        extra_files={"other.bed": "chr1\t90\t210\toth1\nchr1\t150\t300\toth2\n"})
    both = json.loads(go(tmp_path / "b", manifest).read_text())
    d = both["delimitation"]["by_source"]
    assert set(d) == {"SrcPi", "SrcOther"}
    assert d["SrcPi"] == base["delimitation"]["by_source"]["SrcPi"]
    assert d["SrcOther"]["summary"]["n_clusters"] == 2
    assert both["delimitation"]["global"] == base["delimitation"]["global"]


def test_run_aborts_on_contig_mismatch(tmp_path, monkeypatch):
    manifest, seen = setup(tmp_path, monkeypatch, truth_chrom="NC_000001.1")
    with pytest.raises(ValueError, match="NC_000001.1"):
        go(tmp_path, manifest)
    assert "run" not in seen.get("order", [])
    assert not (tmp_path / "results" / "SRR9" / "metrics.json").exists()


def test_contig_check_uses_raw_fasta_first_word_before_reads_and_pipeline(tmp_path, monkeypatch):
    # truth on "chr2" matches the first word of a header ">chr2 description"; "chr3" does not
    manifest, seen = setup(tmp_path, monkeypatch, truth_chrom="chr3")
    with pytest.raises(ValueError, match="chr3"):
        go(tmp_path, manifest)
    assert "order" not in seen  # neither reads download nor pipeline ran


def test_run_rejects_unknown_truth_source(tmp_path, monkeypatch):
    manifest, _ = setup(tmp_path, monkeypatch, truth_sources=["pirna.bed", "nope.bed"])
    with pytest.raises(ValueError, match="nope.bed"):
        go(tmp_path, manifest)


def test_run_rejects_unknown_accession(tmp_path, monkeypatch):
    manifest, _ = setup(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="UNKNOWN"):
        go(tmp_path, manifest, acc="UNKNOWN")


def test_metrics_only_rebuilds_from_existing_work_without_running(tmp_path, monkeypatch):
    manifest, _ = setup(tmp_path, monkeypatch)
    path = go(tmp_path, manifest)
    first = json.loads(path.read_text())
    path.unlink()
    (path.parent / "benchmark_run.json").unlink()

    def boom(*a, **k):
        raise AssertionError("must not run the pipeline or download")

    for name in ("run_library", "fetch_genome", "fetch_library"):
        monkeypatch.setattr(orchestrate, name, boom)
    out = orchestrate.rebuild_metrics("SRR9", manifest, tmp_path / "results")
    assert out == path
    assert json.loads(path.read_text()) == first
    assert (path.parent / "benchmark_run.json").exists()


def test_metrics_only_without_work_dir_fails_clearly(tmp_path, monkeypatch):
    manifest, _ = setup(tmp_path, monkeypatch)
    with pytest.raises(FileNotFoundError, match="SRR9"):
        orchestrate.rebuild_metrics("SRR9", manifest, tmp_path / "results")


def test_metrics_subcommand(tmp_path, monkeypatch, capsys):
    from benchmark.__main__ import main
    manifest, _ = setup(tmp_path, monkeypatch)
    path = go(tmp_path, manifest)
    path.unlink()
    rc = main(["metrics", "--library", "SRR9", "--manifest", str(manifest),
               "--results-dir", str(tmp_path / "results")])
    assert rc == 0 and path.exists()
