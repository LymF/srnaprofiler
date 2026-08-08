"""The 5' base null is the locus's own DNA composition, not 0.25.

A cut placed at random in the locus yields a 5'U as often as the locus has T on
the read's strand. Assuming 0.25 ignores that: a T-rich transposon would look
1U-biased with no selection at all. Measured, the genomic null is 0.272 in
Drosophila and 0.332 in Arabidopsis -- and reads' own positions 3-5, which look
like a cheaper proxy, disagree with the DNA in plants (0.255 against 0.332).

The null is strand-weighted because pct_1u pools reads from both strands: a read
on the minus strand starts with T where the plus strand carries A.
"""
import pytest

from samio import base_null


def test_all_T_on_the_plus_strand():
    # every plus-strand read starts with T; every minus-strand read with A
    assert base_null("TTTTTTTTTT", "T", frac_top=1.0) == 1.0
    assert base_null("TTTTTTTTTT", "T", frac_top=0.0) == 0.0


def test_balanced_strands_average_the_two():
    assert base_null("TTTTTTTTTT", "T", frac_top=0.5) == 0.5


def test_composition_is_counted_exactly():
    # ACGTACGTAC holds 2 T in 10 bases
    assert base_null("ACGTACGTAC", "T", frac_top=1.0) == pytest.approx(0.2, abs=0.001)


def test_a_T_rich_locus_raises_the_bar():
    """This is the point: 40% T means 40% of random cuts start with U."""
    seq = "T" * 40 + "A" * 20 + "C" * 20 + "G" * 20
    assert base_null(seq, "T", frac_top=1.0) == pytest.approx(0.40, abs=0.01)


def test_A_null_mirrors_T_on_the_other_strand():
    seq = "T" * 40 + "A" * 20 + "C" * 20 + "G" * 20
    assert base_null(seq, "A", frac_top=0.0) == pytest.approx(0.40, abs=0.01)


def test_lowercase_counts_and_N_is_excluded():
    """N cannot be cut into a 5' base, so it leaves the denominator."""
    assert base_null("ttttNNNN", "T", frac_top=1.0) == pytest.approx(1.0, abs=0.001)
    assert base_null("ttttcccc", "T", frac_top=1.0) == pytest.approx(0.5, abs=0.001)


def test_empty_sequence_falls_back_to_uniform():
    assert base_null("", "T", frac_top=1.0) == 0.25


def test_the_genome_actually_reaches_the_score(tmp_path):
    """Wiring test: -g must change s_1u end to end.

    It did not, once: the call in run_classify was left without genome_path and
    a whole validation run came back byte-identical before anyone noticed.
    """
    import types
    import pysam
    from classify import run_classify

    fa = tmp_path / "g.fa"
    # a T-rich locus: under a 0.25 null it looks 1U-biased with no selection
    fa.write_text(">c\n" + ("TTTTTTTTAC" * 60) + "\n")
    pysam.faidx(str(fa))

    sam = tmp_path / "in.sam"
    with open(sam, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c\tLN:600\n@PG\tID:x\tPN:other\n")
        for i in range(120):
            p = 10 + (i % 60) * 3
            fh.write(f"r{i}\t0\tc\t{p}\t255\t24M\t*\t0\t0\t{'T' + 'ACGT' * 5 + 'ACG'}\t"
                     f"{'I' * 24}\tXM:i:1\n")

    def run(genome):
        args = types.SimpleNamespace(
            sam=str(sam), genome=genome, prefix="g" if genome else "n",
            outdir=str(tmp_path), mincov=0.0, pad=50, mirna_loci=None,
            config="rules.yaml", norm=1_000_000, threads=1)
        run_classify(args)
        import csv
        rows = list(csv.DictReader(
            open(tmp_path / f"{'g' if genome else 'n'}.clusters.tsv"), delimiter="\t"))
        return max(float(r["s_1u"]) for r in rows)

    without = run(None)
    with_genome = run(str(fa))
    assert with_genome < without, (with_genome, without)
