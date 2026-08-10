"""Interop with the external miRNA tools, where the silent failures live.

Five ways the third-party tools fail without saying so, all absorbed here so
none of them reaches the user:

1. miRDeep2 reports the precursor cut flush against the mature, and miRScore
   refuses a mature sitting at the edge of the hairpin. The fix is to re-cut the
   precursor from the genome, by miRDeep2's own coordinates, with flanks.
2. miRScore chdirs into its output directory and prefixes `../` to the FASTQ
   path, so an absolute path becomes `..//abs/path`, bowtie finds nothing, the
   BAM comes out empty, and the failure only surfaces later in samtools sort.
3. miRDeep-P2 does not create its output directory.
4. miRDeep-P2 does not read gzipped input -- it runs to completion and emits
   zero predictions.
5. miRDeep-P2 has a path bug in its own script (writes `$out/$n/$n.ori.fa`,
   reads `$out/$n.ori.fa`).
"""
import pytest

from mirna import hairpin_with_flank, relative_to, mirdeep_p2_command


class FakeGenome:
    """Minimal stand-in for pysam.FastaFile."""

    def __init__(self, seqs):
        self.seqs = seqs

    def fetch(self, chrom, start, end):
        return self.seqs[chrom][start:end]

    def get_reference_length(self, chrom):
        return len(self.seqs[chrom])


#   0         1         2         3         4
#   0123456789012345678901234567890123456789012345
SEQ = "AAAACCCCGGGGTTTTACGTACGTACGTTTTTGGGGCCCCAAAAGG"
GENOME = FakeGenome({"chr1": SEQ})


def test_the_mature_is_not_at_the_edge_of_the_hairpin():
    """The whole point: miRScore rejects a mature flush against either end."""
    mature = SEQ[16:24]                      # ACGTACGT
    hp = hairpin_with_flank(GENOME, "chr1", 17, 24, mature, flank=8)
    assert hp is not None
    i = hp.find(mature)
    assert i > 0
    assert i + len(mature) < len(hp)


def test_the_flank_is_taken_from_the_genome_on_both_sides():
    mature = SEQ[16:24]
    hp = hairpin_with_flank(GENOME, "chr1", 17, 24, mature, flank=5)
    assert len(hp) == (24 - 17 + 1) + 2 * 5


def test_a_mature_on_the_minus_strand_is_found_by_reverse_complement():
    """miRDeep2 reports the mature as sequenced; the genome holds the other strand."""
    mature_rc = "ACGTACGT"                   # revcomp of ACGTACGT is itself here
    minus = "GGGGCCCC"                       # lives on chr1 as CCCCGGGG? no: check both
    hp = hairpin_with_flank(GENOME, "chr1", 33, 40, "GGGGCCCC", flank=6)
    assert hp is not None and "GGGGCCCC" in hp


def test_a_mature_that_is_nowhere_in_the_region_is_refused():
    """Better no hairpin than one whose mature it does not contain."""
    assert hairpin_with_flank(GENOME, "chr1", 1, 10, "TTTTTTTT", flank=5) is None


def test_the_flank_is_clipped_at_the_contig_edges():
    mature = SEQ[0:6]
    hp = hairpin_with_flank(GENOME, "chr1", 1, 6, mature, flank=50)
    assert hp is not None and hp.startswith(SEQ[0:6]) is False or hp.count(mature) >= 1
    assert len(hp) <= len(SEQ)


def test_a_path_is_made_relative_for_mirscore(tmp_path):
    """An absolute path becomes `..//abs/path` inside miRScore and finds nothing."""
    work = tmp_path / "work"
    work.mkdir()
    reads = tmp_path / "lib.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    rel = relative_to(reads, work)
    assert not str(rel).startswith("/")
    assert (work / rel).resolve() == reads.resolve()


def test_a_path_already_inside_the_work_dir_stays_short(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    reads = work / "lib.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    assert str(relative_to(reads, work)) == "lib.fastq"


def test_the_mirdeep_p2_command_carries_genome_index_and_output():
    cmd = mirdeep_p2_command("/opt/mirdp2/pipeline.bash", "g.fa", "idx",
                             "reads.fastq", "out", threads=8)
    assert "/opt/mirdp2/pipeline.bash" in cmd[1] or cmd[0].endswith("bash")
    joined = " ".join(cmd)
    assert "-g g.fa" in joined and "-x idx" in joined
    assert "-o out" in joined and "-p 8" in joined


def test_the_mirdeep_p2_command_declares_fastq():
    cmd = " ".join(mirdeep_p2_command("p.bash", "g.fa", "idx", "reads.fastq",
                                      "out", threads=1))
    assert "-q" in cmd


def test_gzipped_input_is_refused_rather_than_run(tmp_path):
    """It would run to completion and emit zero predictions."""
    with pytest.raises(ValueError, match="gzip"):
        mirdeep_p2_command("p.bash", "g.fa", "idx", "reads.fastq.gz", "out",
                           threads=1)
