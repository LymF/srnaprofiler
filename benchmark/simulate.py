"""Spike-in simulation: synthetic small-RNA loci injected into a real background library."""
import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

PARAMS_PATH = Path(__file__).with_name("sim_params.yaml")
CLASSES = ("piRNA", "siRNA", "miRNA", "phasiRNA", "hcsiRNA")
COMP = str.maketrans("ACGTN", "TGCAN")
TRIES = 300
# a read is (start, end, strand, tag); coordinates are genomic, 0-based half-open
Read = tuple[int, int, str, str]


@dataclass
class SimSpec:
    classes: dict[str, int]
    seed: int


@dataclass
class SimResult:
    fastq: Path
    truth_bed: Path
    genome: Path


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def load_genome(path: Path) -> dict[str, str]:
    genome: dict[str, list[str]] = {}
    name = ""
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            genome[name] = []
        elif line.strip():
            genome[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in genome.items()}


def load_params() -> dict:
    raw = yaml.safe_load(PARAMS_PATH.read_text())
    return {cls: {k: v["value"] for k, v in block.items()} for cls, block in raw.items()}


def five_prime_base(seq: str, s: int, e: int, strand: str) -> str:
    return seq[s] if strand == "+" else seq[e - 1].translate(COMP)


def draw_length(rng, dist: dict[int, float]) -> int:
    return int(rng.choice(list(dist), p=np.array(list(dist.values())) / sum(dist.values())))


def draw_start(rng, seq, lo, hi, length, strand, base, fits=lambda s: True):
    """Random start in [lo, hi-length] with the wanted 5' base (None = any) and fits(start).

    Returns (start, satisfied); falls back to the last draw when nothing satisfies the constraints."""
    for _ in range(TRIES):
        s = int(rng.integers(lo, hi - length + 1))
        if (base is None or five_prime_base(seq, s, s + length, strand) == base) and fits(s):
            return s, True
    return s, False


def place_loci(rng, genome, n, length, taken):
    """n loci of `length`, inside contigs, N-free, disjoint from each other and from `taken`."""
    chroms = [c for c in genome if len(genome[c]) >= length]
    weights = np.array([len(genome[c]) for c in chroms], dtype=float)
    placed = []
    for _ in range(n):
        for _ in range(100000):
            chrom = chroms[int(rng.choice(len(chroms), p=weights / weights.sum()))]
            start = int(rng.integers(0, len(genome[chrom]) - length + 1))
            end = start + length
            clash = any(c == chrom and start < e and s < end for c, s, e, *_ in taken + placed)
            if "N" not in genome[chrom][start:end] and not clash:
                placed.append((chrom, start, end))
                break
        else:
            raise RuntimeError(f"cannot place {n} loci of {length} nt without overlap or N")
    return placed


def pirna_reads(rng, p, seq, start, end, strand, n) -> list[Read]:
    reads: list[Read] = []
    opp = "-" if strand == "+" else "+"
    for k in range(n):
        length = draw_length(rng, p["length_dist"])
        want_u = rng.random() < p["first_u_prob"]
        pingpong = rng.random() < p["pingpong_fraction"]
        plen = draw_length(rng, p["length_dist"])
        ov = p["pingpong_overlap"]

        def partner(s):   # opposite-strand read whose 5' end overlaps this read's 5' end by ov nt
            if strand == "+":
                return s + ov - plen, s + ov
            return s + length - ov, s + length - ov + plen

        def fits(s):
            ps, pe = partner(s)
            return not pingpong or (start <= ps and pe <= end)
        # the partner's position `ov` is the complement of the primary 5' base, so 1U gives A10
        s, ok = draw_start(rng, seq, start, end, length, strand, "T" if want_u or pingpong else None, fits)
        reads.append((s, s + length, strand, f"r{k}"))
        if pingpong and ok:
            ps, pe = partner(s)
            reads.append((ps, pe, opp, f"r{k}pp"))
    return reads


def sirna_reads(rng, p, seq, start, end, n) -> list[Read]:
    reads: list[Read] = []
    for k in range(n):
        length = draw_length(rng, p["length_dist"])
        strand = "+" if rng.random() < p["plus_strand_prob"] else "-"
        s, _ = draw_start(rng, seq, start, end, length, strand, None)
        reads.append((s, s + length, strand, f"r{k}"))
    return reads


def hcsirna_reads(rng, p, seq, start, end, n) -> list[Read]:
    reads: list[Read] = []
    for k in range(n):
        length = draw_length(rng, p["length_dist"])
        strand = "+" if rng.random() < p["plus_strand_prob"] else "-"
        base = "A" if rng.random() < p["first_a_prob"] else None
        s, _ = draw_start(rng, seq, start, end, length, strand, base)
        reads.append((s, s + length, strand, f"r{k}"))
    return reads


def mirna_reads(rng, p, start, end, strand, n) -> list[Read]:
    """5p mature arm plus a less abundant 3p star forming a duplex with 2-nt 3' overhangs."""
    mlen, a = p["mature_length"], p["five_prime_offset"]
    # arm positions pair as a + b + (mlen - 1 - overhang) == locus_length - 1 (hairpin symmetry)
    b = (end - start) - 1 - (mlen - 1 - p["overhang"]) - a
    arms = {"5p": (a, a + mlen), "3p": (b, b + mlen)}

    def genomic(ls, le):
        return (start + ls, start + le) if strand == "+" else (end - le, end - ls)
    reads: list[Read] = []
    for k in range(n):
        arm = "3p" if rng.random() < p["star_fraction"] else "5p"
        s, e = genomic(*arms[arm])
        reads.append((s, e, strand, f"r{k}{arm}"))
    return reads


def phasirna_reads(rng, p, start, end, n) -> list[Read]:
    """21-nt reads in a fixed register; minus-strand duplex partners are shifted by the 3' overhang."""
    period, ov = p["period"], p["overhang"]
    cycles = (end - start) // period
    reads: list[Read] = []
    for k in range(n):
        cycle = int(rng.integers(1, cycles))
        shift = 0 if rng.random() < p["in_phase_fraction"] else int(rng.integers(1, period))
        if rng.random() < p["minus_strand_prob"]:
            s = start + period * cycle - ov - shift
            s = max(s, start)
            reads.append((s, s + period, "-", f"r{k}"))
        else:
            s = min(start + period * cycle + shift, end - period)
            reads.append((s, s + period, "+", f"r{k}"))
    return reads


def hairpin(p, locus_seq: str, strand: str) -> str:
    """Genome stretch (as stored on the + strand) folding into a hairpin with the mirna_reads arms.

    The 5p arm and 5' flank stay as in the locus; the 3p arm side is their reverse complement
    (3p base b+j pairs with 5p base a+19-j, so 5p/3p carry 2-nt 3' overhangs); a short second stem
    closes a ~10-nt loop."""
    t = list(locus_seq if strand == "+" else revcomp(locus_seq))
    last = len(t) - 1
    pairs = list(range(0, p["five_prime_offset"] + p["mature_length"] - p["overhang"]))
    arm_end = p["five_prime_offset"] + p["mature_length"]
    pairs += list(range(arm_end, arm_end + 8))   # second stem, after a 2x2 internal loop at the arm ends
    for x in pairs:
        t[last - x] = t[x].translate(COMP)
    out = "".join(t)
    return out if strand == "+" else revcomp(out)


def make_reads(rng, cls, p, seq, start, end, strand, n) -> list[Read]:
    if cls == "piRNA":
        return pirna_reads(rng, p, seq, start, end, strand, n)
    if cls == "siRNA":
        return sirna_reads(rng, p, seq, start, end, n)
    if cls == "miRNA":
        return mirna_reads(rng, p, start, end, strand, n)
    if cls == "phasiRNA":
        return phasirna_reads(rng, p, start, end, n)
    return hcsirna_reads(rng, p, seq, start, end, n)


def open_binary(path: Path):
    return gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")


def simulate(spec: SimSpec, background_fastq: Path, genome: Path, out_dir: Path, *,
             reads_per_locus: int | None = None, adapter: str | None = None) -> SimResult:
    unknown = sorted(set(spec.classes) - set(CLASSES))
    if unknown:
        raise ValueError(f"unknown class(es) {unknown}; supported: {list(CLASSES)}")
    params = load_params()
    n_reads = reads_per_locus or params["defaults"]["reads_per_locus"]
    rng = np.random.default_rng(spec.seed)
    seqs = load_genome(genome)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    loci: list[tuple] = []        # (chrom, start, end, name, cls, strand)
    records: list[str] = []
    for cls in spec.classes:
        p = params[cls]
        for i, (chrom, start, end) in enumerate(place_loci(rng, seqs, spec.classes[cls], p["locus_length"], loci)):
            strand = ("+" if rng.random() < 0.5 else "-") if cls in ("piRNA", "miRNA") else "."
            loci.append((chrom, start, end, f"{cls}_{i}", cls, strand))
            if cls == "miRNA":
                hp = hairpin(p, seqs[chrom][start:end], strand)
                seqs[chrom] = seqs[chrom][:start] + hp + seqs[chrom][end:]
            for s, e, st, tag in make_reads(rng, cls, p, seqs[chrom], start, end, strand, n_reads):
                seq = seqs[chrom][s:e] if st == "+" else revcomp(seqs[chrom][s:e])
                seq += adapter or ""
                records.append(f"@sim_{cls}_{i}_{tag}\n{seq}\n+\n{'I' * len(seq)}\n")

    fastq = out_dir / "sim.fastq.gz"
    with open(fastq, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        last = b"\n"
        with open_binary(Path(background_fastq)) as bg:
            while chunk := bg.read(1 << 20):
                gz.write(chunk)
                last = chunk[-1:]
        if last != b"\n":
            gz.write(b"\n")
        gz.write("".join(records).encode())

    truth_bed = out_dir / "sim_truth.bed"
    truth_bed.write_text("".join(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\n" for r in loci))
    for cls in spec.classes:
        lines = [f"{c}\t{s}\t{e}\t{n}\t0\t{st}\n" for c, s, e, n, k, st in loci if k == cls]
        (out_dir / f"sim_truth.{cls}.bed").write_text("".join(lines))
    genome_out = out_dir / "sim_genome.fa"
    with open(genome_out, "w") as fh:
        for chrom, seq in seqs.items():
            fh.write(f">{chrom}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    return SimResult(fastq, truth_bed, genome_out)
