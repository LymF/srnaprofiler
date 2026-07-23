"""Per-class plotting: read the classifier TSV + BAM, plot signals grouped by class."""
from collections import Counter, defaultdict
from dataclasses import dataclass

_REQUIRED = ("chrom", "start", "end", "name", "reads", "dicer_call", "cls")


@dataclass
class Locus:
    chrom: str
    start: int
    end: int
    name: str
    cls: str
    reads: float
    dicer_call: str


def load_clusters(tsv_path):
    """Parse the reporter TSV into Locus rows, skipping rows with a blank `cls`."""
    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        missing = [c for c in _REQUIRED if c not in idx]
        if missing:
            raise ValueError(f"clusters TSV missing columns: {missing}")
        loci = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cls = f[idx["cls"]]
            if not cls:
                continue
            loci.append(Locus(
                chrom=f[idx["chrom"]], start=int(f[idx["start"]]),
                end=int(f[idx["end"]]), name=f[idx["name"]],
                cls=cls, reads=float(f[idx["reads"]]),
                dicer_call=f[idx["dicer_call"]]))
    return loci


def group_by_class(loci, classes=None):
    """Group Locus rows by `cls`; restrict to `classes` if given."""
    grouped = defaultdict(list)
    for l in loci:
        if classes is None or l.cls in classes:
            grouped[l.cls].append(l)
    return dict(grouped)


def top_n_loci(loci, n):
    """Return the `n` loci with the largest `reads`, descending; `n <= 0` -> []."""
    if n <= 0:
        return []
    return sorted(loci, key=lambda l: l.reads, reverse=True)[:n]


def dominant_dicer(loci, default=21):
    """Return the most common integer `dicer_call` among loci, ignoring 'N'/''."""
    calls = Counter(int(l.dicer_call) for l in loci if l.dicer_call not in ("N", ""))
    return calls.most_common(1)[0][0] if calls else default


def parse_region(spec, name=None):
    """Parse 'chrom:start-end' (0-based half-open) into a Locus."""
    try:
        chrom, coords = spec.rsplit(":", 1)
        start_s, end_s = coords.split("-")
        start, end = int(start_s), int(end_s)
    except ValueError:
        raise ValueError(f"malformed region (expected chrom:start-end): {spec!r}")
    return Locus(chrom=chrom, start=start, end=end, name=name or spec,
                 cls="", reads=0.0, dicer_call="N")
