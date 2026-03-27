#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
from collections import defaultdict, Counter
from typing import List, Tuple

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
from Bio import SeqIO


BASES = ["A", "C", "G", "U"]

COLORS = {
    "A": "#4daf4a",
    "C": "#377eb8",
    "G": "#ff7f00",
    "U": "#e41a1c"
}

FONT = FontProperties(family="DejaVu Sans", weight="bold")


# ------------------------------------------------------------
# IO
# ------------------------------------------------------------

def read_sam_sequences_with_strand(path: Path) -> List[Tuple[str, str]]:

    out = []

    with open(path) as f:

        for line in f:

            if line.startswith("@"):
                continue

            parts = line.rstrip().split("\t")

            if len(parts) < 10:
                continue

            try:
                flag = int(parts[1])
            except:
                continue

            seq = parts[9].upper()

            if seq == "*" or not seq:
                continue

            seq = seq.replace("T", "U")

            strand = "antisense" if (flag & 16) else "sense"

            out.append((seq, strand))

    return out


def read_fasta_sequences(path: Path) -> List[str]:

    return [
        str(r.seq).upper().replace("T", "U")
        for r in SeqIO.parse(str(path), "fasta")
    ]


# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------

def uniq_list(lst):

    seen = set()
    out = []

    for s in lst:

        if s not in seen:
            seen.add(s)
            out.append(s)

    return out


def trim_first_n(lst, n):

    if not n or n <= 0:
        return lst

    return [s[:n] for s in lst]


# ------------------------------------------------------------
# logo math
# ------------------------------------------------------------

def seqs_to_matrix(seqs):

    L = len(seqs[0])
    matrix = []

    for i in range(L):

        col = Counter(s[i] for s in seqs if len(s) == L)
        total = sum(col.values())

        freqs = {b: col.get(b, 0)/total for b in BASES}

        H = -sum(p * math.log2(p) for p in freqs.values() if p > 0)

        R = 2 - H  # RNA alphabet

        heights = {b: freqs[b] * R for b in BASES}

        matrix.append(heights)

    return matrix


# ------------------------------------------------------------
# drawing engine
# ------------------------------------------------------------

def draw_letter(ax, letter, x, y, height):

    ax.text(
        x,
        y,
        letter,
        fontsize=height * 120,   # escala visual
        ha="center",
        va="bottom",
        color=COLORS[letter],
        fontweight="bold"
    )


def draw_logo(matrix, outfile, title):

    fig, ax = plt.subplots(figsize=(12, 4))

    for pos, heights in enumerate(matrix):

        letters = sorted(BASES, key=lambda b: heights[b])

        y = 0

        for base in letters:

            h = heights[base]

            if h > 0:
                draw_letter(ax, base, pos, y, h)
                y += h

    ax.set_xlim(-0.5, len(matrix)-0.5)
    ax.set_ylim(0, 2)
    ax.set_xlabel("Position")
    ax.set_ylabel("Information (bits)")
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(outfile)
    plt.savefig(outfile, format="svg")
    
    plt.close()


def generate_logo(seqs, outfile, title):

    if len(seqs) < 2:
        return

    matrix = seqs_to_matrix(seqs)

    draw_logo(matrix, outfile, title)

    print(f"[OK] {outfile}")


# ------------------------------------------------------------
# main pipeline
# ------------------------------------------------------------

def main():

    ap = argparse.ArgumentParser()

    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--outprefix", required=True)

    ap.add_argument("--min", type=int, default=20)
    ap.add_argument("--max", type=int, default=30)

    ap.add_argument("--uniq", action="store_true")
    ap.add_argument("--limit-positions", type=int, default=11)

    args = ap.parse_args()

    in_path = Path(args.input)
    outprefix = Path(args.outprefix)

    outprefix.parent.mkdir(parents=True, exist_ok=True)

    total = 0

    # ---------- SAM ----------

    if in_path.suffix.lower() == ".sam":

        seqs = read_sam_sequences_with_strand(in_path)

        buckets = defaultdict(lambda: {"sense": [], "antisense": []})

        for seq, strand in seqs:

            L = len(seq)

            if args.min <= L <= args.max:
                buckets[L][strand].append(seq)

        for L in range(args.min, args.max + 1):

            group = buckets.get(L)

            if not group:
                continue

            sense = group["sense"]
            anti = group["antisense"]

            if args.uniq:
                sense = uniq_list(sense)
                anti = uniq_list(anti)

            sense = trim_first_n(sense, args.limit_positions)
            anti = trim_first_n(anti, args.limit_positions)

            combined = sense + anti

            if len(combined) >= 2:

                generate_logo(
                    combined,
                    outprefix.parent / f"{outprefix.name}.{L}.combined.pdf",
                    f"{L} nt combined"
                )

                total += 1

            if len(sense) >= 2:

                generate_logo(
                    sense,
                    outprefix.parent / f"{outprefix.name}.{L}.sense.pdf",
                    f"{L} nt sense"
                )

                total += 1

            if len(anti) >= 2:

                generate_logo(
                    anti,
                    outprefix.parent / f"{outprefix.name}.{L}.antisense.pdf",
                    f"{L} nt antisense"
                )

                total += 1

    # ---------- FASTA ----------

    elif in_path.suffix.lower() in [".fa", ".fasta"]:

        seqs = read_fasta_sequences(in_path)

        if args.uniq:
            seqs = uniq_list(seqs)

        buckets = defaultdict(list)

        for s in seqs:

            L = len(s)

            if args.min <= L <= args.max:
                buckets[L].append(s)

        for L in range(args.min, args.max + 1):

            reads = buckets.get(L, [])

            if len(reads) < 2:
                continue

            reads = trim_first_n(reads, args.limit_positions)

            generate_logo(
                reads,
                outprefix.parent / f"{outprefix.name}.{L}.pdf",
                f"{L} nt"
            )

            total += 1

    if total == 0:
        print("No logos generated.")


if __name__ == "__main__":
    main()