"""Sinais de biogênese: ping-pong z-score e phasing."""
import statistics


def pingpong_zscore(counts, target=10, lo=1, hi=20):
    """Z-score do overlap 5'-5' no `target` (default 10 nt) contra o fundo.

    counts: {overlap_length: n_pares}. Fundo = offsets em [lo,hi] menos o target.
    """
    background = [counts.get(o, 0) for o in range(lo, hi + 1) if o != target]
    peak = counts.get(target, 0)
    if len(background) < 2:
        return 0.0, peak
    mean = statistics.fmean(background)
    if mean == 0:
        # Fundo sem nenhuma leitura: dado insuficiente para um z-score.
        return 0.0, peak
    sd = statistics.pstdev(background)
    if sd == 0:
        # Fundo constante (variância nula): aproximação de Poisson (sd ~ sqrt(mean))
        # evita divisão por zero e ainda detecta picos que destoam do fundo.
        sd = mean ** 0.5
    return round((peak - mean) / sd, 3), peak
