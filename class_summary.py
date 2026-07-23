"""Matplotlib QC summary panels from the classifier .clusters.tsv (no R needed)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def write_summary(tsv_path, out_dir, prefix):
    os.makedirs(out_dir, exist_ok=True)
    df = pd.read_csv(tsv_path, sep="\t")
    df = df[df["cls"].notna() & (df["cls"].astype(str) != "")]
    written = []

    def save(fig, tag):
        path = os.path.join(out_dir, f"{prefix}.summary.{tag}.pdf")
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        written.append(path)

    # counts per class
    fig, ax = plt.subplots(figsize=(5, 4))
    df["cls"].value_counts().plot.bar(ax=ax, color="grey")
    ax.set_ylabel("loci"); ax.set_title("Loci per class")
    save(fig, "counts")

    # dicer_call distribution per class
    if "dicer_call" in df:
        fig, ax = plt.subplots(figsize=(6, 4))
        (df.groupby("cls")["dicer_call"].value_counts()
           .unstack(fill_value=0).T.plot.bar(ax=ax))
        ax.set_xlabel("DicerCall"); ax.set_ylabel("loci")
        save(fig, "dicer")

    # scalar signals per class (box)
    signals = [c for c in ("pingpong_frac", "phase_score", "pct_1u", "gc", "multimap_frac")
               if c in df]
    if signals:
        fig, axes = plt.subplots(1, len(signals), figsize=(3 * len(signals), 4))
        if len(signals) == 1:
            axes = [axes]
        for ax, col in zip(axes, signals):
            groups = [g[col].dropna().values for _, g in df.groupby("cls")]
            ax.boxplot(groups, labels=list(df.groupby("cls").groups))
            ax.set_title(col)
        save(fig, "signals")

    # confidence + review rate per class
    if {"confidence", "review"} <= set(df.columns):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 4))
        for cls, g in df.groupby("cls"):
            a1.hist(g["confidence"].dropna(), bins=10, alpha=0.5, label=cls)
        a1.set_title("Confidence"); a1.legend()
        rate = df.groupby("cls")["review"].apply(
            lambda s: (s.astype(str) == "True").mean())
        rate.plot.bar(ax=a2, color="salmon"); a2.set_title("Review rate")
        save(fig, "review")

    return written
