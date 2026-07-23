"""Per-class plotting: read the classifier TSV + BAM, plot signals grouped by class."""
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
