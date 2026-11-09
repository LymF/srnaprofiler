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
    size_purity: float = 0.0
    dicer_call_2: str = ""
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
    # evidence scores, log10(1/p) per signature
    s_1u: float = 0.0
    s_5a: float = 0.0
    s_fractop: float = 0.0
    s_fractop_bal: float = 0.0
    s_fractop_bid: float = 0.0
    s_pingpong: float = 0.0
    s_phase: float = 0.0
    # Effect strengths, in [0,1] -- what `score` is built from. TSV columns too (the ML gate reads them).
    e_1u: float = 0.0
    e_5a: float = 0.0
    e_fractop_bid: float = 0.0
    e_pingpong: float = 0.0
    e_phase: float = 0.0
    score: Optional[float] = 0.0   # None when the class bypassed the scorer
    gc: float = 0.0
    # biogenesis
    pingpong_z: float = 0.0
    pingpong_pairs10: float = 0.0
    pingpong_frac: float = 0.0
    phase_score: float = 0.0
    phase_pvalue: float = 1.0
    phase_number: int = 0
    phase_ratio: float = 0.0
    overlap_dup_frac: float = 0.0
    # repetition
    multimap_frac: float = 0.0
    unique_frac: float = 0.0
    # bp from the first read start to the last read end (0 with no reads); None
    # when unknown (e.g. replayed .clusters.tsv without the column)
    read_span: Optional[int] = None
    # per-strand size histogram: fraction of the locus's 1/NH mass per length
    szp_18: float = 0.0
    szp_19: float = 0.0
    szp_20: float = 0.0
    szp_21: float = 0.0
    szp_22: float = 0.0
    szp_23: float = 0.0
    szp_24: float = 0.0
    szp_25: float = 0.0
    szp_26: float = 0.0
    szp_27: float = 0.0
    szp_28: float = 0.0
    szp_29: float = 0.0
    szp_30: float = 0.0
    szp_31: float = 0.0
    szp_32: float = 0.0
    szp_33: float = 0.0
    szp_34: float = 0.0
    szp_35: float = 0.0
    szm_18: float = 0.0
    szm_19: float = 0.0
    szm_20: float = 0.0
    szm_21: float = 0.0
    szm_22: float = 0.0
    szm_23: float = 0.0
    szm_24: float = 0.0
    szm_25: float = 0.0
    szm_26: float = 0.0
    szm_27: float = 0.0
    szm_28: float = 0.0
    szm_29: float = 0.0
    szm_30: float = 0.0
    szm_31: float = 0.0
    szm_32: float = 0.0
    szm_33: float = 0.0
    szm_34: float = 0.0
    szm_35: float = 0.0
    # read-start shape (None when unknown)
    start_entropy: Optional[float] = None
    top_start_frac: Optional[float] = None
    starts_per_kb: Optional[float] = None
    # proTRAC density significance vs the chromosome; None when unknown
    density_pvalue: Optional[float] = None
    # classification
    cls: str = ""
    confidence: float = 0.0
    review: bool = False
    class_scores: dict = field(default_factory=dict)
    evidence: str = ""
    p_real: Optional[float] = None   # ML mode only: probability the locus is a real small-RNA locus

    @property
    def length(self) -> int:
        return self.end - self.start
