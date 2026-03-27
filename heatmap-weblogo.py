#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import Counter


def read_sam(sam, min_len, max_len):

    sense = []
    antisense = []

    with open(sam) as f:
        for line in f:
            if line.startswith("@"):
                continue

            cols = line.strip().split("\t")
            flag = int(cols[1])
            seq = cols[9]

            L = len(seq)
            if L < min_len or L > max_len:
                continue

            seq = seq.upper().replace("T", "U")[:11]

            if len(seq) < 11:
                continue

            if flag & 16:
                antisense.append(seq)
            else:
                sense.append(seq)

    return sense, antisense


def freq_matrix(seqs):

    bases = ["A","C","G","U"]
    mat = np.zeros((4,11))

    for i in range(11):
        col = [s[i] for s in seqs if len(s) > i]
        count = Counter(col)
        total = sum(count.values()) or 1

        for j,b in enumerate(bases):
            mat[j,i] = count[b] / total

    return mat


def plot_pair(sense_mat, anti_mat, out):

    fig, axes = plt.subplots(
        2,1,
        figsize=(6,2.5),
        sharex=True
    )

    bases = ["A","C","G","U"]

    for ax, mat, title in zip(
        axes,
        [sense_mat, anti_mat],
        ["Sense strand","Antisense strand"]
    ):

        im = ax.imshow(mat, aspect="auto", cmap="viridis")

        ax.set_yticks(range(4))
        ax.set_yticklabels(bases)
        ax.set_title(title, fontsize=9)

        ax.axvline(0, linestyle="--", linewidth=0.25)
        ax.axvline(9, linestyle="--", linewidth=0.25)

    axes[1].set_xticks(range(11))
    axes[1].set_xticklabels(range(1,12))
    axes[1].set_xlabel("Position")

    cbar = fig.colorbar(im, ax=axes, shrink=0.8)
    cbar.set_label("Frequency")

    plt.tight_layout()

    with PdfPages(out) as pdf:
        pdf.savefig(fig)

    plt.close()


def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("-i", required=True)
    ap.add_argument("-o", required=True)
    ap.add_argument("--min", type=int, default=21)
    ap.add_argument("--max", type=int, default=22)

    args = ap.parse_args()

    sense, anti = read_sam(args.i, args.min, args.max)

    if not sense or not anti:
        print("Warning: one strand empty")

    sense_mat = freq_matrix(sense)
    anti_mat = freq_matrix(anti)

    plot_pair(sense_mat, anti_mat, args.o)

    print("PDF saved:", args.o)


if __name__ == "__main__":
    main()