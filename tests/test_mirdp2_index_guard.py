"""miRDeep-P2's fourth silent failure: no bowtie index.

The pipeline takes the index with `-x`, and `mirna.py` passed `args.index_prefix`
straight through. Omitted, it became the literal string "None": bowtie aligned
nothing, every .bst came out empty, and miRDeep-P2 ran to completion, wrote
`end...` and reported zero predictions with a success exit code.

Measured on a real Arabidopsis library: 5 minutes, 0 predictions, no error.
"""
import pytest

from mirna import require_bowtie_index


def make_index(tmp_path, prefix="idx", exts=None):
    exts = exts or ['.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt',
                    '.rev.1.ebwt', '.rev.2.ebwt']
    for ext in exts:
        (tmp_path / f"{prefix}{ext}").write_text("")
    return tmp_path / prefix


def test_a_complete_index_is_accepted(tmp_path):
    index = make_index(tmp_path)
    assert require_bowtie_index(index) == index


def test_no_index_at_all_is_refused(tmp_path):
    """The case that cost a five-minute run and produced a silent zero."""
    with pytest.raises(RuntimeError, match="--index-prefix"):
        require_bowtie_index(None)


def test_an_empty_string_is_refused_like_none(tmp_path):
    with pytest.raises(RuntimeError, match="--index-prefix"):
        require_bowtie_index("")


def test_a_path_with_no_index_beside_it_is_refused(tmp_path):
    with pytest.raises(RuntimeError, match="no bowtie index"):
        require_bowtie_index(tmp_path / "nothing_here")


def test_a_half_built_index_is_refused(tmp_path):
    """bowtie-build killed midway leaves some files; the run must not start."""
    index = make_index(tmp_path, exts=['.1.ebwt', '.2.ebwt'])
    with pytest.raises(RuntimeError, match="no bowtie index"):
        require_bowtie_index(index)


def test_the_large_index_flavour_is_accepted(tmp_path):
    """bowtie-build writes .ebwtl for a large genome; both are valid."""
    index = make_index(tmp_path, exts=['.1.ebwtl', '.2.ebwtl', '.3.ebwtl',
                                       '.4.ebwtl', '.rev.1.ebwtl', '.rev.2.ebwtl'])
    assert require_bowtie_index(index) == index


def test_the_refusal_names_the_path_it_looked_at(tmp_path):
    missing = tmp_path / "wrong_place"
    with pytest.raises(RuntimeError, match=str(missing)):
        require_bowtie_index(missing)
