"""FIX 5+6 (TDD): atribuição read->região por SOBREPOSIÇÃO DE SPAN.

A cobertura em clustering.cluster_regions é construída pelo span inteiro do read.
A atribuição em features.extract deve ser consistente: um read que começa ANTES
de `start` mas entra na ilha precisa ser contado (antes era descartado por usar só
a posição leftmost).
"""
from features import extract

CONFIG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def test_read_starting_before_region_start_is_counted():
    # edge_read: pos0=90, 20 nt -> span [90,110), começa antes de start=100.
    # inside_read: pos0=120, 21 nt -> claramente dentro.
    region = ("chr1", 100, 160)
    recs = extract("tests/data/edge_reads.sam", [region], lib_size=1_000_000, config=CONFIG)
    r = recs[0]
    # ambos os reads sobrepõem a região por span; o edge_read NÃO pode ser descartado
    assert r.reads == 2
    assert r.distinct_seqs == 2
