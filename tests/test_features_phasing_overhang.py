"""FIX 3 (TDD): correção do overhang de 2 nt no phasing.

Produtos fásicos de Dicer têm overhang 3' de 2 nt, então o 5' da fita - fica
deslocado (+2) em relação ao registro da fita +. Sem a correção, um locus
fásico bifita racha a abundância em dois registros (ratio ~0.5). Com a correção,
sense e antisense caem no MESMO registro dominante (ratio == 1.0).
"""
from samio import read_alignments
from features import build_record

CONFIG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def test_overhang_pools_both_strands_in_one_register():
    alns = list(read_alignments("tests/data/phasing_overhang.sam"))
    # region_start=1000, pl detectado = 21 (todos os reads têm 21 nt)
    rec = build_record(("chr1", 1000, 1200), alns, lib_size=1_000_000, config=CONFIG)
    assert rec.dicer_call == "21"
    # 3 ciclos fásicos, cobertos por + e - juntos após a correção de 2 nt
    assert rec.phase_number == 3
    # sem a correção, +/- cairiam em registros diferentes -> ratio == 0.5;
    # com a correção os dois pooling no mesmo registro -> ratio == 1.0
    assert rec.phase_ratio == 1.0
