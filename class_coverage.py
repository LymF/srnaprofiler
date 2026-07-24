"""IGV-style mirror coverage track (matplotlib, native PDF, no Ghostscript).

The look replicates `extras/rna_density_plot.R` (ggplot theme_classic): white
background, no grid, only the left/bottom spines, large fonts, and the mirror
fill — sense up in brown (#996631), antisense down in blue (#0099CC).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_SENSE = "#996631"      # brown, drawn upward
_ANTISENSE = "#0099CC"  # blue, drawn downward


def plot_coverage_track(sense, antisense, chrom, start, out_path):
    """Mirror coverage over genomic position, in the rna_density_plot.R style."""
    x = [start + i for i in range(len(sense))]
    neg = [-a for a in antisense]
    m = max(max(sense, default=0), max(antisense, default=0), 1) * 1.1

    fig, ax = plt.subplots(figsize=(12, 5))
    # filled areas with a thin outline, like ggplot geom_density(stat="identity")
    ax.fill_between(x, sense, color=_SENSE, linewidth=0)
    ax.fill_between(x, neg, color=_ANTISENSE, linewidth=0)
    ax.plot(x, sense, color=_SENSE, linewidth=0.4)
    ax.plot(x, neg, color=_ANTISENSE, linewidth=0.4)

    ax.set_ylim(-m, m)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("RNA density", fontsize=18)
    ax.set_xlabel(chrom, fontsize=18)
    ax.set_title(f"{chrom}:{start}-{start + len(sense)}", fontsize=18)
    ax.tick_params(labelsize=13)

    # theme_classic: keep only the left and bottom axis lines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")

    # minimal frameless legend (up = sense, down = antisense)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=_SENSE, label="Sense (+)"),
                       Patch(color=_ANTISENSE, label="Antisense (-)")],
              frameon=False, fontsize=13, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
