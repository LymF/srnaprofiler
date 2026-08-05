"""Engine outputs: TSV (feature matrix), GFF3, and JSON provenance."""
import json

_COLS = [
    "chrom", "start", "end", "name", "reads", "reads_rpm", "distinct_seqs",
    "major_rna", "dicer_call", "size_purity", "dicer_call_2", "frac_in_window", "frac_top", "strand_call",
    "pct_1u", "pct_10a", "pct_5a", "gc", "pingpong_z", "pingpong_pairs10", "pingpong_frac", "overlap_dup_frac",
    "phase_score", "phase_pvalue", "phase_number", "phase_ratio", "multimap_frac",
    "cls", "confidence", "score",
    "s_1u", "s_5a", "s_fractop", "s_fractop_bal", "s_pingpong", "s_phase", "s_density",
    "review", "evidence",
]


def _cell(value):
    """Serializes a field for the TSV, neutralizing TAB/newline that would break the row."""
    return str(value).replace("\t", " ").replace("\n", " ")


def write_tsv(records, path):
    with open(path, "w") as fh:
        fh.write("\t".join(_COLS) + "\n")
        for r in records:
            fh.write("\t".join(_cell(getattr(r, c)) for c in _COLS) + "\n")


def write_gff3(records, path):
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for r in records:
            attrs = f"ID={r.name};score={r.score};confidence={r.confidence};review={r.review}"
            fh.write("\t".join([
                r.chrom, "srnaprofiler", r.cls or "unclassified",
                str(r.start + 1), str(r.end), str(r.score),
                r.strand_call, ".", attrs,
            ]) + "\n")


def write_run_json(params, path):
    # default=str keeps this robust: any non-JSON-native value (e.g. a Path or
    # a stray callable) is stringified instead of crashing the whole run.
    with open(path, "w") as fh:
        json.dump(params, fh, indent=2, default=str)
