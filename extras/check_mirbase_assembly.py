#!/usr/bin/env python3
"""Do miRBase's coordinates belong to this genome?

miRBase publishes one GFF3 per species against whatever assembly it used --
`osa` is on MSU7, `ath` on TAIR10, `dme` on Release 6 -- and a run that mixes
assemblies fails silently: --mirna-loci carves loci at the wrong places, asserts
miRNA on them, and every validation still returns a number.

The check is the only one that cannot be fooled: pull the sequence the GFF
points at out of the genome and compare it with the precursor miRBase says
lives there. Matching names prove nothing; matching sequence proves everything.

Usage:
    python3 extras/check_mirbase_assembly.py genome.fa species.gff3 hairpin.fa osa

The species prefix selects the precursors (osa, ath, dme, ...). Chromosome
naming often differs between the two files ("1" vs "Chr1" vs "NC_029256.1"),
so names are matched loosely before giving up.
"""
import gzip
import re
import sys


def _open(path):
    """Open by content, not by name: servers hand out plain text named .gz."""
    with open(path, "rb") as fh:
        gzipped = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gzipped else open(path)


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks)
    return seqs


def read_precursors(path, prefix):
    """Precursor sequences by miRBase accession (MI......), as DNA."""
    out = {}
    name = None
    chunks = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(chunks).upper().replace("U", "T")
                fields = line[1:].split()
                name = fields[1] if (len(fields) > 1
                                     and fields[0].startswith(prefix + "-")) else None
                chunks = []
            elif name:
                chunks.append(line.strip())
    if name:
        out[name] = "".join(chunks).upper().replace("U", "T")
    return out


def read_gff(path):
    rows = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "miRNA_primary_transcript":
                continue
            attrs = dict(re.findall(r"(\w+)=([^;]+)", f[8]))
            rows.append((f[0], int(f[3]) - 1, int(f[4]), f[6],
                         attrs.get("ID", ""), attrs.get("Name", "")))
    return rows


def match_contigs(rows, precursors, genome, probes=40):
    """Find each GFF chromosome in the genome by content, not by name.

    Names are useless across sources -- the same chromosome is "chr1" in
    miRBase, "Chr1" in sRNAanno and "NC_003070.9" in RefSeq -- and a name table
    would need one entry per assembly anyone ever uses. Instead, for a few
    precursors per chromosome, ask which contig actually carries that sequence
    at those coordinates. A chromosome that no contig matches is reported as
    unmatched, which is itself the answer: wrong assembly.
    """
    by_chrom = {}
    for chrom, start, end, strand, acc, name in rows:
        by_chrom.setdefault(chrom, []).append((start, end, strand, acc))

    mapping = {}
    for chrom, entries in by_chrom.items():
        votes = {}
        for start, end, strand, acc in entries[:probes]:
            seq = precursors.get(acc)
            if not seq:
                continue
            for contig, dna in genome.items():
                if len(dna) < end:
                    continue
                got = dna[start:end].upper()
                if strand == "-":
                    got = revcomp(got)
                if got == seq:
                    votes[contig] = votes.get(contig, 0) + 1
                    break
        if votes:
            mapping[chrom] = max(votes, key=votes.get)
    return mapping


_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def main(genome_path, gff_path, hairpin_path, prefix):
    genome = read_fasta(genome_path)
    precursors = read_precursors(hairpin_path, prefix)
    rows = read_gff(gff_path)
    if not rows:
        sys.exit(f"no miRNA_primary_transcript in {gff_path}")
    if not precursors:
        sys.exit(f"no {prefix}- precursors in {hairpin_path}")

    mapping = match_contigs(rows, precursors, genome)
    unmapped = sorted({r[0] for r in rows} - set(mapping))
    print(f"{len(rows)} precursors in the GFF, {len(precursors)} sequences, "
          f"{len(genome)} contigs in the genome")
    if unmapped:
        print(f"chromosomes with no counterpart in the genome: "
              f"{', '.join(unmapped[:8])}")

    checked = exact = 0
    misses = []
    for chrom, start, end, strand, acc, name in rows:
        seq = precursors.get(acc)
        contig = mapping.get(chrom)
        if not seq or not contig:
            continue
        got = genome[contig][start:end].upper()
        if strand == "-":
            got = revcomp(got)
        checked += 1
        if got == seq:
            exact += 1
        elif len(misses) < 5:
            misses.append((name, chrom, start, len(seq), len(got)))

    if not checked:
        sys.exit("nothing could be compared: no shared accession or contig")
    pct = 100.0 * exact / checked
    print(f"\nprecursors whose genome sequence matches miRBase exactly: "
          f"{exact}/{checked} ({pct:.1f}%)")
    for name, chrom, start, want, got in misses:
        print(f"   mismatch: {name} at {chrom}:{start} "
              f"(precursor {want} nt, genome slice {got} nt)")

    print()
    if pct >= 95:
        print("SAME ASSEMBLY -- the coordinates can be used as they are.")
    elif pct >= 50:
        print("PARTIAL: most match, some do not. Look at the mismatches before "
              "trusting the annotation; a liftover may be needed for part of it.")
    else:
        print("DIFFERENT ASSEMBLY -- do NOT pass this GFF to --mirna-loci. It "
              "would carve loci at the wrong places and assert miRNA on them.")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    main(*sys.argv[1:5])
