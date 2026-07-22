"""Saídas do motor: TSV (matriz de features), GFF3 e proveniência JSON."""
import json

_COLS = [
    "chrom", "start", "end", "name", "reads", "reads_rpm", "distinct_seqs",
    "major_rna", "dicer_call", "frac_in_window", "frac_top", "strand_call",
    "pct_1u", "pct_10a", "gc", "pingpong_z", "pingpong_pairs10", "overlap19_frac",
    "phase_score", "phase_pvalue", "phase_number", "phase_ratio", "multimap_frac",
    "cls", "confidence", "review", "evidence",
]


def _cell(value):
    """Serializa um campo para o TSV, neutralizando TAB/newline que quebrariam a linha."""
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
            attrs = f"ID={r.name};confidence={r.confidence};review={r.review}"
            fh.write("\t".join([
                r.chrom, "srnaprofiler", r.cls or "unclassified",
                str(r.start + 1), str(r.end), str(r.confidence),
                r.strand_call, ".", attrs,
            ]) + "\n")


def write_run_json(params, path):
    with open(path, "w") as fh:
        json.dump(params, fh, indent=2)
