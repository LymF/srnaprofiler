"""Sequence logos via logomaker (pure matplotlib; no external CLI, no Ghostscript)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils import log

# Nucleotide colours (A green, C blue, G gold, U red) — the classic look.
_COLORS = {"A": "#00A000", "C": "#0000C0", "G": "#F2A900", "U": "#C00000"}


def plot_logo(sequences, title, out_path):
    """Render an information-content (bits) logo of equal-length RNA `sequences`.

    Returns False (with a warning) if logomaker is unavailable or there is nothing
    to plot. The output format is inferred from `out_path`'s extension by matplotlib
    (pdf/svg/png/eps) and needs no Ghostscript.
    """
    if not sequences:
        log(f"no sequences for logo ({title}); skipping", "WARN")
        return False
    try:
        import logomaker
    except ImportError:
        log("logomaker not installed; skipping logo (pip install logomaker)", "WARN")
        return False

    matrix = logomaker.alignment_to_matrix(sequences, to_type="information")
    fig, ax = plt.subplots(figsize=(max(4, 0.35 * len(matrix)), 2.4))
    logomaker.Logo(matrix, ax=ax, color_scheme=_COLORS)
    ax.set_title(title)
    ax.set_ylabel("bits")
    ax.set_xlabel("position (nt)")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return True
