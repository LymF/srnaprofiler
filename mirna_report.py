"""The mirna module's deliverables: one record per prediction, and the files.

The engine knows where a locus is and miRScore knows whether it is a miRNA;
neither knows both. This is where the two are joined.
"""
from dataclasses import dataclass

NA = "NA"

MIRNA_COLS = [
    "name", "chrom", "start", "end", "strand", "engine",
    "mature_seq", "mature_len", "star_seq", "star_len",
    "precursor_seq", "precursor_len",
    "mature_start", "mature_end", "star_start", "star_end",
    "result", "flags", "mismatches", "bulges",
    "mature_reads", "star_reads", "total_reads", "precision",
    "mirdeep_score", "mirdeep_prob", "rfam_alert", "randfold",
    "mirbase_hit", "mirbase_same_seed",
]

# miRScore answers in RNA and the predictions arrive in DNA; one alphabet wins,
# and it is the genome's, so the FASTAs and the BED describe the same thing.
_SEQUENCE_FIELDS = ("star_seq", "precursor_seq")

# miRScore column -> our column. Ours are snake_case and say what they hold.
_FROM_MIRSCORE = {
    "msSeq": "star_seq", "msLen": "star_len",
    "precSeq": "precursor_seq", "precLen": "precursor_len",
    "mStart": "mature_start", "mStop": "mature_end",
    "msStart": "star_start", "msStop": "star_end",
    "flags": "flags", "mismatches": "mismatches", "bulges": "bulges",
    "mReads": "mature_reads", "msReads": "star_reads",
    "totReads": "total_reads", "precision": "precision",
}

# miRDeep2 column -> our column. Passthrough: absent means NA, not an error.
_FROM_MIRDEEP = {
    "miRDeep2 score": "mirdeep_score",
    "estimated probability that the miRNA candidate is a true positive": "mirdeep_prob",
    "rfam alert": "rfam_alert",
    "significant randfold p-value": "randfold",
    "miRBase miRNA": "mirbase_hit",
    "example miRBase miRNA with the same seed": "mirbase_same_seed",
}


@dataclass
class MirnaRecord:
    """One prediction, with wherever it got to."""
    name: str
    chrom: str
    start: int                  # 0-based, half-open, like .clusters.tsv
    end: int
    strand: str
    engine: str
    mature_seq: str
    mature_len: int
    result: str                 # Pass | Fail | NotScored | NotValidated
    star_seq: str = NA
    star_len: str = NA
    precursor_seq: str = NA
    precursor_len: str = NA
    mature_start: str = NA
    mature_end: str = NA
    star_start: str = NA
    star_end: str = NA
    flags: str = NA
    mismatches: str = NA
    bulges: str = NA
    mature_reads: str = NA
    star_reads: str = NA
    total_reads: str = NA
    precision: str = NA
    mirdeep_score: str = NA
    mirdeep_prob: str = NA
    rfam_alert: str = NA
    randfold: str = NA
    mirbase_hit: str = NA
    mirbase_same_seed: str = NA


def _verdict(row):
    """miRScore's own result, normalised to Pass/Fail."""
    return "Pass" if str(row.get("result", "")).strip().lower().startswith("pass") else "Fail"


def build_records(predictions, engine, mirscore_rows=None, submitted=None,
                  mirdeep_rows=None):
    """One record per prediction, joined by name.

    `mirscore_rows=None` means miRScore did not run at all -- which is a
    different thing from a locus it rejected, so it gets its own verdict.
    """
    seen = set()
    for p in predictions:
        if p["name"] in seen:
            raise RuntimeError(
                f"prediction name {p['name']} appears more than once; the "
                f"table joins miRScore by name and would multiply rows")
        seen.add(p["name"])

    by_name = dict(mirdeep_rows or {})

    out = []
    for p in predictions:
        rec = MirnaRecord(
            name=p["name"], chrom=p["chrom"],
            start=p["start"] - 1, end=p["end"],       # 1-based inclusive in
            strand=p.get("strand", "."), engine=engine,
            mature_seq=p["mature"], mature_len=len(p["mature"]),
            result="NotValidated",
        )
        if mirscore_rows is not None:
            row = mirscore_rows.get(p["name"])
            if row is not None:
                rec.result = _verdict(row)
                for src, dst in _FROM_MIRSCORE.items():
                    val = row.get(src, NA) or NA
                    if dst in _SEQUENCE_FIELDS and val != NA:
                        val = val.upper().replace("U", "T")
                    setattr(rec, dst, val)
            elif submitted is not None and p["name"] in submitted:
                rec.result = "NotScored"
                rec.flags = "absent from miRScore results"
            else:
                rec.result = "NotScored"
                rec.flags = "precursor could not be re-cut from the genome"
        for src, dst in _FROM_MIRDEEP.items():
            val = by_name.get(p["name"], {}).get(src)
            if val not in (None, ""):
                setattr(rec, dst, val)
        out.append(rec)
    return out
