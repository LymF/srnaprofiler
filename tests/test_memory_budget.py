"""Memory budgets for the two steps that scale with --threads.

Both were found the hard way on an 8 GB machine:

  * `samtools sort` defaults to 768 MB PER THREAD, so --create-bam --threads 12
    asks for ~9 GB and dies with "couldn't allocate memory for bam_mem".
  * classify holds one whole chromosome's alignments per thread, so 4 threads
    on Arabidopsis got the process killed by the OOM killer.

Neither is a disk or a drvfs problem: they fail the same way on any filesystem.
"""
import mapping
import classify


def test_sort_share_shrinks_as_threads_grow():
    assert mapping.sort_memory_per_thread_mb(1, budget_mb=4096) == 4096
    assert mapping.sort_memory_per_thread_mb(12, budget_mb=4096) == 341
    # whatever the thread count, the total stays inside the budget
    for t in (1, 4, 12, 32):
        assert mapping.sort_memory_per_thread_mb(t, budget_mb=4096) * t <= 4096


def test_sort_share_has_a_floor():
    # a tiny share would make samtools spill constantly; better to overshoot
    assert mapping.sort_memory_per_thread_mb(64, budget_mb=1024) >= 64


def test_sort_budget_is_probed_when_not_given():
    assert mapping.sort_memory_per_thread_mb(4) > 0


def test_sort_command_carries_m_and_T(monkeypatch, tmp_path):
    seen = []
    monkeypatch.setattr(mapping, "run_command", lambda cmd, **kw: seen.append(cmd))
    monkeypatch.setattr(mapping, "which", lambda _: "/usr/bin/samtools")
    mapping.convert_sam_to_bam(tmp_path / "in.sam", tmp_path / "out.bam", threads=12)
    sort_cmd = seen[0]
    assert "-m" in sort_cmd and "-T" in sort_cmd


def test_classify_threads_capped_by_the_sum_of_the_biggest_refs():
    five_big = [10_000_000] * 5          # ~3.2 GB each
    # 2 GB cannot afford even one, and the floor is 1 (better slow than dead)
    assert classify.safe_thread_count(8, five_big, available_mb=2_000) == 1
    # 16 GB * 0.8 headroom = 12.8 GB, 3.2 GB each -> three of them
    assert classify.safe_thread_count(8, five_big, available_mb=16_000) == 3


def test_one_big_reference_does_not_throttle_the_small_ones():
    # sizing by the LARGEST alone would allow a single thread here; summing the
    # N largest lets the small references run in parallel behind the big one
    mixed = [10_000_000] + [100_000] * 7
    assert classify.safe_thread_count(8, mixed, available_mb=6_000) == 8


def test_never_more_threads_than_references():
    # one reference cannot use two workers, however much memory there is
    assert classify.safe_thread_count(2, [1000], available_mb=64_000) == 1
    assert classify.safe_thread_count(2, [1000, 1000], available_mb=64_000) == 2
