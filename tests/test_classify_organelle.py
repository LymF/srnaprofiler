"""Organelle contigs are dropped by name, only when asked.

Chloroplast and mitochondrion give 336 loci and 7.9% of the read mass in the
Arabidopsis library (27 loci, 0.6% in Drosophila), almost all called hcsiRNA.
Plants have no organellar siRNA pathway: that is degradation being labelled.

Detecting them from the FASTA description was considered and dropped -- miRDeep2
refuses spaces in the identifier, so any genome prepared for it loses the word,
and an assembly from elsewhere carries no description at all. It would work on
our two genomes and fail silently on the first one outside them.
"""
import pytest

from classify import select_contigs, parse_contig_list


def test_named_contigs_are_removed():
    kept, dropped = select_contigs(["chr1", "chr2", "NC_000932.1"],
                                   ["NC_000932.1"])
    assert kept == ["chr1", "chr2"]
    assert dropped == ["NC_000932.1"]


def test_nothing_is_removed_without_the_flag():
    kept, dropped = select_contigs(["chr1", "NC_000932.1"], [])
    assert kept == ["chr1", "NC_000932.1"]
    assert dropped == []


def test_a_contig_that_is_not_there_is_an_error():
    """A typo would silently analyse the organelle and report a clean run."""
    with pytest.raises(ValueError, match="chrC"):
        select_contigs(["chr1", "chr2"], ["chrC"])


def test_no_guessing_from_the_name():
    """"chloroplast" in a contig name means nothing without the flag."""
    kept, dropped = select_contigs(["chr1", "chloroplast"], [])
    assert "chloroplast" in kept


def test_parses_a_comma_separated_list():
    assert parse_contig_list("NC_000932.1,NC_037304.1") == ["NC_000932.1",
                                                            "NC_037304.1"]


def test_tolerates_spaces_and_empty_fields():
    assert parse_contig_list(" NC_000932.1 , ,NC_037304.1,") == ["NC_000932.1",
                                                                 "NC_037304.1"]


def test_no_list_is_no_contigs():
    assert parse_contig_list(None) == []
