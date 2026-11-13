"""EXPERIMENTAL: predicted ORFs for the genomic-context figure (report --predict-orfs).

A removable spike: this module, the report flags, the gm "o" key and one JS track. ORFs are
predicted from the genome sequence alone, never enter an annotation verdict, a class, a
confidence or the review flag.

find_orfs: ATG-to-stop on the six frames of a sequence (standard stops TAA/TAG/TGA, start only
ATG). Per stop codon the LONGEST ORF is kept (the most upstream ATG after the previous stop in
that frame), so no nested ORFs. The reported span includes the stop codon. ORFs that touch the
sequence edge are not reported: an ORF needs a stop codon of its own and a previous stop in
the same frame inside the sequence (otherwise its true start may lie outside). A codon with a
base other than A, C, G, T (N, IUPAC codes) breaks the ORF and forgets the previous stop.
"""
import os

import pysam

import report_genome
from utils import log

_STOPS = frozenset(("TAA", "TAG", "TGA"))
_ACGT = frozenset("ACGT")
_RC = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp(seq):
    # every other character (N, IUPAC codes) becomes N: it breaks an ORF either way
    return "".join(c if c in _ACGT else "N" for c in seq.upper().translate(_RC))[::-1]


def _scan(seq, min_aa):
    """[(start, end, frame)] of the ORFs of one strand, in the coordinates of `seq`."""
    out = []
    n = len(seq)
    for frame in range(3):
        anchored, first = False, None
        for i in range(frame, n - 2, 3):
            codon = seq[i:i + 3]
            if codon in _STOPS:
                if anchored and first is not None and (i - first) // 3 >= min_aa:
                    out.append((first, i + 3, frame))
                anchored, first = True, None
            elif codon == "ATG":
                if anchored and first is None:
                    first = i
            elif not _ACGT.issuperset(codon):
                anchored, first = False, None
    return out


def find_orfs(seq, offset, min_aa=100):
    """[(start, end, strand, frame)], genomic 0-based half-open (offset is where seq starts).

    The minus strand is translated back to forward coordinates. frame is the codon phase
    (0, 1 or 2) of the ORF start counted from the 5' end of `seq` on the strand it is on
    (for the minus strand, from the 3' end of seq). min_aa counts the amino acids, the stop
    not included.
    """
    seq = seq.upper()
    n = len(seq)
    out = [(offset + a, offset + b, "+", f) for a, b, f in _scan(seq, min_aa)]
    out += [(offset + n - b, offset + n - a, "-", f) for a, b, f in _scan(_revcomp(seq), min_aa)]
    return sorted(out)


def _cds(layers, chrom, ws, we):
    """CDS features of the gene layers in the window, straight from the layer index."""
    out = []
    for la in layers or ():
        if la.kind != "gene":
            continue
        idx = la.index.get(chrom)
        if idx is not None:
            out.extend((f.start, f.end) for f in idx.overlaps(ws, we) if f.ftype.lower() == "cds")
    return out


def window_orfs(fasta, chrom, ws, we, layers, min_aa):
    """([[start, end, strand, aa], ...] by start, n dropped by the cap) for one window.

    ORFs overlapping an annotated CDS are left out; at most FEATURE_CAP are kept, longest first."""
    we = min(we, fasta.get_reference_length(chrom))
    if we <= ws:
        return [], 0
    cds = _cds(layers, chrom, ws, we)
    rows = []
    for s, e, st, _ in find_orfs(fasta.fetch(chrom, ws, we), ws, min_aa):
        if not any(s < ce and e > cs for cs, ce in cds):
            rows.append([s, e, st, (e - s) // 3 - 1])
    dropped = max(0, len(rows) - report_genome.FEATURE_CAP)
    if dropped:
        rows = sorted(rows, key=lambda r: (-r[3], r[0]))[:report_genome.FEATURE_CAP]
    rows.sort(key=lambda r: (r[0], r[1]))
    return rows, dropped


def check_genome(path):
    """None when the FASTA is usable, else the error message."""
    if not path:
        return "--predict-orfs needs --genome FASTA"
    if not os.path.isfile(path):
        return f"genome FASTA not found: {path}"
    if not os.path.isfile(path + ".fai"):
        return (f"genome FASTA {path} has no .fai index; create it with "
                f"`samtools faidx {path}`")
    return None


def add_orfs(gm, clusters, layers, genome_path, min_aa):
    """Adds "o" (and mo["o"] when capped) to gm[cluster name] for each cluster, in place."""
    warned = set()
    with pysam.FastaFile(genome_path) as fa:
        refs = set(fa.references)
        for c in clusters:
            g = gm.get(c.name)
            if g is None:
                continue
            if c.chrom not in refs:
                if c.chrom not in warned:
                    warned.add(c.chrom)
                    log(f"report: contig {c.chrom} is not in the genome FASTA; no predicted "
                        f"ORFs there", "WARN")
                continue
            rows, dropped = window_orfs(fa, c.chrom, g["w"][0], g["w"][1], layers, min_aa)
            if rows:
                g["o"] = rows
            if dropped:
                g.setdefault("mo", {})["o"] = dropped
