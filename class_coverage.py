"""IGV-style mirror coverage track (matplotlib, native PDF, no Ghostscript)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_coverage_track(sense, antisense, chrom, start, out_path):
    """Sense up (#996631), antisense down (#0099CC), over genomic position."""
    x = [start + i for i in range(len(sense))]
    neg = [-a for a in antisense]
    m = max(max(sense, default=0), max(antisense, default=0), 1) * 1.1
    fig, ax = plt.subplots(figsize=(15, 7))
    ax.fill_between(x, sense, color="#996631", label="Sense (+)")
    ax.fill_between(x, neg, color="#0099CC", label="Antisense (-)")
    ax.set_ylim(-m, m)
    ax.set_xlabel(f"Position on {chrom}", fontsize=14)
    ax.set_ylabel("sRNA density (coverage)", fontsize=14)
    ax.set_title(f"Coverage - {chrom}:{start}-{start + len(sense)}",
                 fontsize=16, fontweight="bold")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.6)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
