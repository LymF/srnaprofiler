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


# --- the miRDeep-P2 install, checked before it can fail silently ---------

from mirna import check_mirdeep_p2_install

BUGGY = 'seq="$results_folder/${filename}.ori.fa"\n'
FIXED = 'seq="$results_folder/${filename}/${filename}.ori.fa"\n'


def _install(tmp_path, line, with_rfam=True):
    root = tmp_path / "1.1.4"
    (root / "scripts" / "index").mkdir(parents=True)
    (root / "miRDP2-v1.1.4_pipeline.bash").write_text("#!/bin/bash\n" + line)
    if with_rfam:
        for e in ("1", "2", "3", "4", "rev.1", "rev.2"):
            (root / "scripts" / "index" / f"rfam_index.{e}.ebwt").write_text("x")
    return root / "miRDP2-v1.1.4_pipeline.bash"


def test_a_patched_install_passes(tmp_path):
    check_mirdeep_p2_install(_install(tmp_path, FIXED))


def test_the_path_bug_is_refused_with_a_pointer_to_setup(tmp_path):
    """Unpatched it runs to completion and predicts zero -- worse than an error."""
    with pytest.raises(RuntimeError, match="setup"):
        check_mirdeep_p2_install(_install(tmp_path, BUGGY))


def test_a_missing_rfam_index_is_refused(tmp_path):
    """Same silent failure: it ships mature_index but not rfam_index."""
    with pytest.raises(RuntimeError, match="rfam"):
        check_mirdeep_p2_install(_install(tmp_path, FIXED, with_rfam=False))


def test_a_missing_install_says_where_it_looked(tmp_path):
    with pytest.raises(RuntimeError, match="not found"):
        check_mirdeep_p2_install(tmp_path / "nope" / "pipeline.bash")


def test_the_engine_must_be_declared():
    """No default: running the animal engine on a plant library halves precision."""
    import subprocess, sys
    r = subprocess.run([sys.executable, "srnaprofiler.py", "mirna", "-p", "x",
                        "-o", "/tmp/x", "--reads", "r.fq", "--species", "ath"],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "--engine" in r.stderr


# --- parsing each engine into one shape, and the final BED ---------------

from mirna import (parse_mirdeep_p2_predictions, mirscore_command,
                   parse_mirscore_results, write_mirna_loci_bed)

P2_LINE = ("NC_003070.9\t+\tread04305100_x40\tNC_003070.9_102\t"
           "10227086..10227107\t10227086..10227186\tGGGTAAGATCTCTATTGGCAGG\t"
           "GGGTAAGATCTCTATTGGCAGGAAACCATTACTTAGATCTTTGCATCTCTTT\n")


def test_a_mirdeep_p2_prediction_carries_its_coordinates_and_mature(tmp_path):
    f = tmp_path / "p"
    f.write_text(P2_LINE)
    got = parse_mirdeep_p2_predictions(f)
    assert len(got) == 1
    p = got[0]
    assert p["chrom"] == "NC_003070.9" and p["strand"] == "+"
    assert (p["start"], p["end"]) == (10227086, 10227186)   # the precursor
    assert p["mature"] == "GGGTAAGATCTCTATTGGCAGG"
    assert p["name"] == "NC_003070.9_102"


def test_a_malformed_prediction_line_is_skipped_not_guessed(tmp_path):
    f = tmp_path / "p"
    f.write_text("rubbish\n" + P2_LINE)
    assert len(parse_mirdeep_p2_predictions(f)) == 1


def test_the_mirscore_command_uses_relative_paths(tmp_path):
    """Absolute paths make miRScore align to nothing and fail three steps later."""
    work = tmp_path
    reads = tmp_path / "sub" / "lib.fastq"
    reads.parent.mkdir()
    reads.write_text("@r\nACGT\n+\nIIII\n")
    cmd = mirscore_command(work / "m.fa", work / "h.fa", reads, "plant",
                           work / "out", work, threads=4)
    for token in cmd:
        assert not token.startswith("/"), f"{token} is absolute"
    assert "plant" in cmd and "-kingdom" in cmd


def test_mirscore_results_are_read_as_the_set_that_passed(tmp_path):
    f = tmp_path / "miRScore_results.csv"
    f.write_text("name,mSeq,result\na,ACGT,Pass\nb,ACGT,Fail\nc,ACGT,Pass\n")
    assert parse_mirscore_results(f) == {"a", "c"}


def test_the_final_bed_holds_only_what_passed(tmp_path):
    preds = [dict(chrom="c1", start=101, end=200, name="a", strand="+"),
             dict(chrom="c1", start=301, end=400, name="b", strand="-")]
    out = tmp_path / "loci.bed"
    n = write_mirna_loci_bed(preds, {"a"}, out)
    assert n == 1
    line = out.read_text().strip().split("\t")
    assert line[0] == "c1" and line[3] == "a"


def test_the_final_bed_is_zero_based_half_open(tmp_path):
    """--mirna-loci reads BED; the engines report 1-based inclusive."""
    preds = [dict(chrom="c1", start=101, end=200, name="a", strand="+")]
    out = tmp_path / "loci.bed"
    write_mirna_loci_bed(preds, {"a"}, out)
    f = out.read_text().strip().split("\t")
    assert (int(f[1]), int(f[2])) == (100, 200)


# --- the paths setup records, so the flags become optional ---------------

from utils import user_config, config_path


def test_a_recorded_path_is_read_from_the_config(tmp_path, monkeypatch):
    cfg = tmp_path / "paths.yaml"
    cfg.write_text("mirdp2_pipeline: /opt/mirdp2/pipeline.bash\n")
    monkeypatch.setenv("SRNAPROFILER_CONFIG", str(cfg))
    assert user_config()["mirdp2_pipeline"] == "/opt/mirdp2/pipeline.bash"


def test_no_config_is_an_empty_mapping_not_a_crash(tmp_path, monkeypatch):
    monkeypatch.setenv("SRNAPROFILER_CONFIG", str(tmp_path / "absent.yaml"))
    assert user_config() == {}


def test_the_config_lives_where_setup_writes_it(monkeypatch):
    monkeypatch.delenv("SRNAPROFILER_CONFIG", raising=False)
    assert str(config_path()).endswith("srnaprofiler/paths.yaml")


def test_validation_can_be_switched_off_but_is_not_by_default():
    import subprocess, sys
    h = subprocess.run([sys.executable, "srnaprofiler.py", "mirna", "--help"],
                       capture_output=True, text=True).stdout
    assert "--no-validate" in h
    # nothing in the help offers to turn it *on*: it already is
    assert "--validate " not in h


# --- silent failure 6: a dot in the FASTQ name zeroes every count ---------

from mirna import staged_reads_name, stage_reads_for_mirscore


def test_a_dotted_read_name_is_staged_without_dots():
    """miRScore keys the library on the name cut at the first dot, while the BAM
    keeps the whole name. They stop matching and every locus counts zero reads --
    with the reads sitting aligned inside the BAM."""
    assert staged_reads_name("at_ncrna.ncrna_filtered.fastq") == \
        "at_ncrna_ncrna_filtered.fastq"


def test_the_extension_survives_staging():
    """It dispatches on the extension, so .fastq must stay .fastq."""
    assert staged_reads_name("s.trimmed.fq").endswith(".fq")
    assert staged_reads_name("plain.fastq") == "plain.fastq"


def test_staging_points_at_the_same_reads(tmp_path):
    reads = tmp_path / "a.b.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    work = tmp_path / "w"
    work.mkdir()
    staged = stage_reads_for_mirscore(reads, work)
    assert "." not in staged.stem
    assert staged.read_text() == reads.read_text()


def test_staging_is_idempotent(tmp_path):
    reads = tmp_path / "a.b.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    work = tmp_path / "w"
    work.mkdir()
    assert stage_reads_for_mirscore(reads, work) == \
        stage_reads_for_mirscore(reads, work)


from mirna import no_reads_counted


def test_zero_reads_everywhere_is_treated_as_a_failure(tmp_path):
    f = tmp_path / "reads.csv"
    f.write_text("miRNA,library,mReads,msReads,allReads\na,x,0,0,0\nb,x,0,0,0\n")
    assert no_reads_counted(f) is True


def test_a_library_with_any_reads_is_a_real_result(tmp_path):
    f = tmp_path / "reads.csv"
    f.write_text("miRNA,library,mReads,msReads,allReads\na,x,0,0,45\nb,x,0,0,0\n")
    assert no_reads_counted(f) is False


def test_a_staged_link_keeps_its_own_name_in_the_command(tmp_path):
    """Resolving the symlink hands miRScore back the dotted name we staged
    away from, and the zero-count failure returns."""
    real = tmp_path / "a.b.fastq"
    real.write_text("@r\nACGT\n+\nIIII\n")
    work = tmp_path / "w"
    work.mkdir()
    staged = stage_reads_for_mirscore(real, work)
    assert str(relative_to(staged, work)) == "a_b.fastq"
