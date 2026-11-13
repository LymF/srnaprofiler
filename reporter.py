"""Engine outputs: TSV (feature matrix), GFF3, and JSON provenance."""
import json

_COLS = [
    "chrom", "start", "end", "name", "reads", "reads_rpm", "distinct_seqs",
    "major_rna", "dicer_call", "size_purity", "dicer_call_2", "frac_in_window", "frac_top", "strand_call",
    "pct_1u", "pct_10a", "pct_5a", "gc", "pingpong_z", "pingpong_pairs10", "pingpong_frac", "overlap_dup_frac",
    "phase_score", "phase_pvalue", "phase_number", "phase_ratio", "multimap_frac", "unique_frac", "density_pvalue", "read_span",
    *[f"szp_{L}" for L in range(18, 36)],
    *[f"szm_{L}" for L in range(18, 36)],
    "start_entropy", "top_start_frac", "starts_per_kb",
    "cls", "confidence", "score", "p_real",
    "s_1u", "s_5a", "s_fractop", "s_fractop_bal", "s_fractop_bid", "s_pingpong", "s_phase",
    "e_1u", "e_5a", "e_fractop_bid", "e_pingpong", "e_phase",
    "review", "evidence",
]


def _cell(value):
    """Serializes a field for the TSV, neutralizing TAB/newline that would break the row."""
    if value is None:
        return "NA"
    return str(value).replace("\t", " ").replace("\n", " ")


def write_tsv(records, path, gate_classes=(), annotated=False):
    """`gate_classes` appends one gate_<class> column (0/1) per class, read off r.gates.

    `annotated` inserts annotation and annotation_support just before review;
    off, the file is byte-identical to a run without annotation layers.
    """
    gate_classes = list(gate_classes)
    cols = list(_COLS)
    if annotated:
        i = cols.index("review")
        cols[i:i] = ["annotation", "annotation_support"]
    head = cols + [f"gate_{c}" for c in gate_classes]
    with open(path, "w") as fh:
        fh.write("\t".join(head) + "\n")
        for r in records:
            cells = [_cell(getattr(r, c)) for c in cols]
            cells += [str(int(bool(r.gates.get(c)))) for c in gate_classes]
            fh.write("\t".join(cells) + "\n")


def write_fragments_tsv(fragments, path):
    """Regions below the locus floor: reported, never silently dropped."""
    with open(path, "w") as fh:
        fh.write("chrom\tstart\tend\treads\tread_span\treason\n")
        for f in fragments:
            fh.write(f"{f['chrom']}\t{f['start']}\t{f['end']}\t"
                     f"{f['reads']:.2f}\t{f['read_span']}\t{f['reason']}\n")


def _gff3_escape(text):
    # "%" first, so the codes below are not escaped a second time.
    for ch, code in (("%", "%25"), ("\t", "%09"), ("\n", "%0A"), ("\r", "%0D"),
                     (";", "%3B"), ("=", "%3D"), (",", "%2C"), ("&", "%26")):
        text = text.replace(ch, code)
    return text


def write_gff3(records, path, annotated=False):
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for r in records:
            # GFF3 column 6 takes a number or ".", never a string.
            score_col = "." if r.score is None else str(r.score)
            attrs = (f"ID={r.name};score={'NA' if r.score is None else r.score};"
                     f"confidence={r.confidence};review={r.review}")
            if annotated:
                attrs += (f";annotation_support={r.annotation_support}"
                          f";annotation={_gff3_escape(r.annotation)}")
            fh.write("\t".join([
                r.chrom, "srnaprofiler", r.cls or "unclassified",
                str(r.start + 1), str(r.end), score_col,
                r.strand_call, ".", attrs,
            ]) + "\n")


def write_run_json(params, path):
    # default=str keeps this robust: any non-JSON-native value (e.g. a Path or
    # a stray callable) is stringified instead of crashing the whole run.
    with open(path, "w") as fh:
        json.dump(params, fh, indent=2, default=str)
