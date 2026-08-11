"""--prefix names files inside --outdir, and now also a tree that gets deleted."""
import pytest

from utils import validate_prefix


@pytest.mark.parametrize("prefix", ["at", "at.run", "at_2026-08-10", "at..x"])
def test_a_plain_prefix_is_accepted(prefix):
    assert validate_prefix(prefix) == prefix


@pytest.mark.parametrize("prefix", ["../foo", "a/b", "/abs", ".", "..", ""])
def test_a_prefix_that_escapes_the_outdir_is_refused(prefix):
    with pytest.raises(ValueError) as e:
        validate_prefix(prefix)
    assert "--prefix" in str(e.value)


def test_the_error_names_the_offending_value():
    """Refused, never rewritten: silently making it `foo` writes elsewhere."""
    with pytest.raises(ValueError) as e:
        validate_prefix("../foo")
    assert "../foo" in str(e.value)


def test_the_cli_refuses_before_any_subcommand_runs(monkeypatch, tmp_path):
    import srnaprofiler

    ran = []
    monkeypatch.setattr(srnaprofiler, "check_dependencies", lambda: True)
    monkeypatch.setattr(srnaprofiler, "print_header", lambda: None)
    monkeypatch.setattr(srnaprofiler, "run_mirna_analysis",
                        lambda args: ran.append(args))
    monkeypatch.setattr("sys.argv", [
        "srnaprofiler", "mirna", "-p", "../escape", "-o", str(tmp_path),
        "--reads", "r.fq", "--species", "ath", "--engine", "mirdeep-p2",
        "-g", "g.fa"])

    with pytest.raises(SystemExit) as e:
        srnaprofiler.main()
    assert e.value.code == 1
    assert ran == []


def test_collect_structures_refuses_a_dest_it_was_not_asked_to_fill(tmp_path):
    """Belt and braces at the delete itself."""
    from mirna_report import collect_structures

    with pytest.raises(ValueError):
        collect_structures(tmp_path / "mirscore", [], tmp_path / "..", tmp_path)


def test_the_module_entry_point_refuses_before_the_analysis_runs(
        monkeypatch, tmp_path, capsys):
    """`python mirna.py -p ../foo` never passes through srnaprofiler.main()."""
    import mirna

    ran = []
    monkeypatch.setattr(mirna, "run_mirna_analysis", lambda args: ran.append(args))
    monkeypatch.setattr("sys.argv", [
        "mirna.py", "-p", "../escape", "-o", str(tmp_path),
        "--reads", "r.fq", "--species", "ath", "-g", "g.fa"])

    with pytest.raises(SystemExit) as e:
        mirna.main()
    assert e.value.code == 1
    assert ran == []
    out = capsys.readouterr()
    assert "../escape" in out.out + out.err


def test_collect_structures_refuses_a_dest_outside_the_outdir(tmp_path):
    """Path() normalises `..` away, so the name check alone lets it through."""
    from mirna_report import collect_structures

    outdir = tmp_path / "out"
    outdir.mkdir()
    victim = tmp_path / "victim.mirna_structures"
    victim.mkdir()
    (victim / "keep.txt").write_text("mine")

    with pytest.raises(ValueError) as e:
        collect_structures(tmp_path / "mirscore", [],
                           outdir / ".." / "victim.mirna_structures", outdir)
    assert "victim.mirna_structures" in str(e.value)
    assert (victim / "keep.txt").exists()


def test_collect_structures_clears_a_stale_file_instead_of_aborting(tmp_path):
    """A file (not a directory) at the dest must not abort a written-out run."""
    from mirna_report import collect_structures

    outdir = tmp_path / "out"
    outdir.mkdir()
    dest = outdir / "at.mirna_structures"
    dest.write_text("not a directory")

    collect_structures(tmp_path / "mirscore", [], dest, outdir)
    assert not dest.exists()


def test_collect_structures_unlinks_a_symlinked_dest_without_touching_it(tmp_path):
    """The link goes; whatever it pointed at stays."""
    from mirna_report import collect_structures

    outdir = tmp_path / "out"
    outdir.mkdir()
    target = tmp_path / "elsewhere"
    target.mkdir()
    (target / "keep.txt").write_text("mine")
    dest = outdir / "at.mirna_structures"
    dest.symlink_to(target)

    collect_structures(tmp_path / "mirscore", [], dest, outdir)
    assert not dest.is_symlink()
    assert (target / "keep.txt").exists()


def test_dropping_deliverables_refuses_paths_outside_the_outdir(tmp_path):
    """The FASTA unlinks are aimed by the same prefix."""
    from mirna import deliverable_paths, drop_validated_deliverables

    outdir = tmp_path / "out"
    outdir.mkdir()
    victim = tmp_path / "victim.mirna.mature.fa"
    victim.write_text("mine")

    paths = deliverable_paths(outdir, "../victim")
    with pytest.raises(ValueError) as e:
        drop_validated_deliverables(paths, outdir)
    assert "victim" in str(e.value)
    assert victim.exists()


def test_dropping_deliverables_clears_a_stale_structures_file(tmp_path):
    from mirna import deliverable_paths, drop_validated_deliverables

    outdir = tmp_path / "out"
    outdir.mkdir()
    paths = deliverable_paths(outdir, "at")
    paths["structures"].write_text("not a directory")

    drop_validated_deliverables(paths, outdir)
    assert not paths["structures"].exists()
