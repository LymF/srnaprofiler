"""Data contract for the classification engine: one record per locus."""
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class ClusterRecord:
    chrom: str
    start: int
    end: int
    name: str = ""
    # abundance
    reads: int = 0
    reads_rpm: float = 0.0
    distinct_seqs: int = 0
    major_rna: str = ""
    major_rna_reads: int = 0
    # size
    dicer_call: str = "N"
    frac_in_window: float = 0.0
    short: int = 0
    long: int = 0
    size_counts: dict = field(default_factory=dict)
    # strand
    frac_top: Optional[float] = None
    strand_call: str = "."
    # composition
    pct_1u: float = 0.0
    pct_10a: float = 0.0
    pct_5a: float = 0.0
    support: float = 0.0
    # evidence scores, log10(1/p) per signature
    s_1u: float = 0.0
    s_5a: float = 0.0
    s_fractop: float = 0.0
    s_pingpong: float = 0.0
    s_phase: float = 0.0
    s_density: float = 0.0
    score: float = 0.0
    gc: float = 0.0
    # biogenesis
    pingpong_z: float = 0.0
    pingpong_pairs10: float = 0.0
    pingpong_frac: float = 0.0
    phase_score: float = 0.0
    phase_pvalue: float = 1.0
    phase_number: int = 0
    phase_ratio: float = 0.0
    overlap19_frac: float = 0.0
    # repetition
    multimap_frac: float = 0.0
    # classification
    cls: str = ""
    confidence: float = 0.0
    review: bool = False
    class_scores: dict = field(default_factory=dict)
    evidence: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start
