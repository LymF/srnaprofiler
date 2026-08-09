#!/usr/bin/env python3
"""Rebuild a reference annotation's coordinates from its own sequence.

A published annotation carries coordinates on whatever assembly its authors
used, and that assembly is often unobtainable, renamed, or silently different.
The failure is invisible: every tool runs, every validator returns a number,
and the number is wrong. Measured on the Tribolium piRNA clusters of
piRNAclusterDB -- validating against the published GTF gave recall 31/61 with
27 clusters holding zero reads and no enrichment at all, which reads as "the
tool does not generalise". The sequences were in the genome the whole time,
somewhere else.

When the reference ships sequences as well as coordinates, the coordinates are
redundant: cut the sequences into short chunks, map them, and let the genome say
where each feature is. What comes out belongs to the assembly in hand by
construction.

Usage:
    python3 extras/remap_reference_by_sequence.py features.fasta bowtie_index out.bed

Each FASTA record becomes one BED line. A feature whose chunks disagree about
where it is -- rearranged between assemblies, or too repetitive to place -- is
dropped and counted, never guessed.
"""
import collections
import os
import statistics
import subprocess
import sys
import tempfile

CHUNK = 60          # long enough to be unique, short enough for bowtie
CHUNKS_PER_FEATURE = 25
MIN_ANCHORS = 3     # unique chunks needed to place a feature
MIN_AGREEMENT = 0.6  # share of anchors that must land on the same contig
MAX_SPREAD = 2000   # bp: median disagreement about the start, above which we refuse


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
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


def write_chunks(seqs, path):
    n = 0
    with open(path, "w") as out:
        for name, seq in seqs.items():
            if len(seq) < CHUNK * 2:
                continue
            step = max(CHUNK, (len(seq) - CHUNK) // CHUNKS_PER_FEATURE)
            for off in range(0, len(seq) - CHUNK, step):
                piece = seq[off:off + CHUNK].upper()
                if piece.count("N") > CHUNK // 10:
                    continue
                out.write(f">{name}|{off}\n{piece}\n")
                n += 1
    return n


def run_bowtie(index, chunks_fa, hits_path, threads=8):
    cmd = ["bowtie", "-f", "-v", "1", "-m", "1", "--best", "-p", str(threads),
           "-x", index, chunks_fa, hits_path]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"bowtie failed:\n{res.stderr}")
    return res.stderr


def place(hits, lengths):
    """One BED row per feature, or a reason it was dropped."""
    rows, dropped = [], []
    for name, anchors in hits.items():
        by_contig = collections.Counter(a[0] for a in anchors)
        contig, count = by_contig.most_common(1)[0]
        if count < MIN_ANCHORS or count / len(anchors) < MIN_AGREEMENT:
            dropped.append((name, "anchors disagree on the contig"))
            continue
        same = [a for a in anchors if a[0] == contig]
        plus = [pos - off for _, pos, off, strand in same if strand == "+"]
        minus = [pos + off for _, pos, off, strand in same if strand == "-"]
        starts = plus if len(plus) >= len(minus) else minus
        if len(starts) < MIN_ANCHORS:
            dropped.append((name, "too few anchors on one strand"))
            continue
        start = int(statistics.median(starts))
        spread = statistics.median([abs(s - start) for s in starts])
        if spread > MAX_SPREAD:
            dropped.append((name, f"anchors imply starts {spread:.0f} bp apart"))
            continue
        rows.append((contig, max(0, start), max(0, start) + lengths[name],
                     name, int(spread)))
    return sorted(rows), dropped


def main(fasta_path, index, out_path, threads=8):
    seqs = read_fasta(fasta_path)
    if not seqs:
        sys.exit(f"no sequence in {fasta_path}")
    lengths = {k: len(v) for k, v in seqs.items()}

    tmp = tempfile.mkdtemp()
    chunks_fa = os.path.join(tmp, "chunks.fa")
    hits_path = os.path.join(tmp, "hits.txt")
    n_chunks = write_chunks(seqs, chunks_fa)
    print(f"{len(seqs)} features, {n_chunks} chunks of {CHUNK} nt")

    stderr = run_bowtie(index, chunks_fa, hits_path, threads)
    for line in stderr.strip().splitlines():
        if "reads" in line.lower():
            print("  " + line.strip())

    hits = collections.defaultdict(list)
    with open(hits_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            name, off = f[0].rsplit("|", 1)
            hits[name].append((f[2], int(f[3]), int(off), f[1]))

    rows, dropped = place(hits, lengths)
    with open(out_path, "w") as out:
        for contig, start, end, name, _ in rows:
            out.write(f"{contig}\t{start}\t{end}\t{name}\n")

    total = sum(e - s for _, s, e, _, _ in rows)
    print(f"\n{len(rows)}/{len(seqs)} features placed -> {out_path}")
    print(f"  they span {total/1e6:.2f} Mb; the sequences total "
          f"{sum(lengths.values())/1e6:.2f} Mb")
    if rows:
        print(f"  median disagreement among a feature's own chunks: "
              f"{statistics.median([r[4] for r in rows])} bp")
    for name, why in dropped[:5]:
        print(f"  dropped {name}: {why}")
    if len(dropped) > 5:
        print(f"  ... and {len(dropped) - 5} more dropped")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(*sys.argv[1:4])
