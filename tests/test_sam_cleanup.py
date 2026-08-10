"""The SAM is a subproduct: once sorted into a BAM, nothing reads it again.

One session left 24 GB of SAM behind. Deleting it is safe only under two rules,
and both are about the same thing -- never destroy the only copy of something:

1. never delete before the final artifact exists AND has been verified;
2. never delete what the user brought (only the SAM this run wrote).
"""
import subprocess
from pathlib import Path

import pytest

from mapping import bam_is_complete, drop_sam


def _bam(tmp_path, name="a.bam"):
    """A real, sorted, indexable BAM -- verification must be able to fail."""
    sam = tmp_path / "in.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:100\n"
        "r1\t0\tc1\t10\t255\t4M\t*\t0\t0\tACGT\tIIII\n")
    bam = tmp_path / name
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    return sam, bam


def test_a_sorted_bam_is_recognised_as_complete(tmp_path):
    _, bam = _bam(tmp_path)
    assert bam_is_complete(bam) is True


def test_a_truncated_bam_is_not_complete(tmp_path):
    """The failure this guards: a run killed mid-sort leaves a plausible BAM."""
    _, bam = _bam(tmp_path)
    data = bam.read_bytes()
    bam.write_bytes(data[:len(data) // 2])
    assert bam_is_complete(bam) is False


def test_an_absent_bam_is_not_complete(tmp_path):
    assert bam_is_complete(tmp_path / "nope.bam") is False


def test_an_empty_bam_is_not_complete(tmp_path):
    p = tmp_path / "empty.bam"
    p.write_bytes(b"")
    assert bam_is_complete(p) is False


def test_the_sam_goes_once_the_bam_is_verified(tmp_path):
    sam, bam = _bam(tmp_path)
    assert drop_sam(sam, bam) is True
    assert not sam.exists()
    assert bam.exists()


def test_the_sam_survives_an_unverifiable_bam(tmp_path):
    """If the BAM cannot be trusted, the SAM is the only copy of the alignment."""
    sam, bam = _bam(tmp_path)
    bam.write_bytes(b"not a bam")
    assert drop_sam(sam, bam) is False
    assert sam.exists()


def test_the_sam_survives_a_missing_bam(tmp_path):
    sam, bam = _bam(tmp_path)
    bam.unlink()
    assert drop_sam(sam, bam) is False
    assert sam.exists()


def test_keeping_the_sam_is_available_on_the_command_line():
    h = subprocess.run(["python3", "srnaprofiler.py", "mapping", "--help"],
                       capture_output=True, text=True).stdout
    assert "--keep-sam" in h


def test_the_bam_is_produced_by_default():
    """The BAM is what classify consumes, so making it opt-in was backwards:
    the default run left tens of GB of SAM and nothing the pipeline reads."""
    h = subprocess.run(["python3", "srnaprofiler.py", "mapping", "--help"],
                       capture_output=True, text=True).stdout
    assert "--no-bam" in h


def test_the_old_flag_is_still_accepted():
    """Existing command lines carry --create-bam; it is now the default."""
    h = subprocess.run(["python3", "srnaprofiler.py", "mapping", "--help"],
                       capture_output=True, text=True).stdout
    assert "--create-bam" in h
    # the option's own description, not the usage line at the top
    assert "default" in h.rsplit("--create-bam", 1)[1][:200]
