"""The per-class aggregate has to be capped or a big class exhausts memory.

classplots holds every member's alignments at once to build the class-level
plots. On the Arabidopsis run the hcsiRNA class has 46431 loci covering most of
the genome, which is effectively the whole BAM as Python objects -- the run hung
without producing a single file. The cap keeps the most abundant loci, which is
where the read mass is.
"""
from classplots import top_n_loci, Locus


def _locus(name, reads):
    return Locus(name=name, chrom="c", start=0, end=100, cls="hcsiRNA",
                 reads=reads, dicer_call="24")


def test_cap_keeps_the_most_abundant():
    loci = [_locus(f"l{i}", reads=float(i)) for i in range(100)]
    kept = top_n_loci(loci, 10)
    assert len(kept) == 10
    assert {l.reads for l in kept} == {float(i) for i in range(90, 100)}


def test_cap_above_the_class_size_keeps_everything():
    loci = [_locus(f"l{i}", reads=float(i)) for i in range(5)]
    assert len(top_n_loci(loci, 1000)) == 5


def test_read_mass_covered_is_computable():
    # what the WARN line reports: the cap is only acceptable if we say how much
    # of the class the plot actually speaks for
    loci = [_locus(f"l{i}", reads=float(i)) for i in range(1, 11)]
    kept = top_n_loci(loci, 3)
    covered = sum(l.reads for l in kept) / sum(l.reads for l in loci)
    assert round(covered, 3) == round((10 + 9 + 8) / 55, 3)
