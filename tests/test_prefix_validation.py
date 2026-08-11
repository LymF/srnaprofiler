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
        collect_structures(tmp_path / "mirscore", [], tmp_path / "..")
