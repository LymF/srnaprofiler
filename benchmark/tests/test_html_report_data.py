"""report_data: TSV reading, depth floor and library summary."""
import report_data as rd

OLD = ("chrom start end name reads reads_rpm distinct_seqs major_rna dicer_call size_purity "
       "dicer_call_2 frac_in_window frac_top strand_call pct_1u pct_10a pct_5a gc pingpong_z "
       "pingpong_pairs10 pingpong_frac overlap_dup_frac phase_score phase_pvalue phase_number "
       "phase_ratio multimap_frac unique_frac cls confidence score s_1u s_5a s_fractop "
       "s_fractop_bal s_fractop_bid s_pingpong s_phase review evidence").split()


def _write(path, rows, cols=OLD):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")
    return str(path)


def _row(name="c1", reads=500, conf=0.5, **kw):
    base = dict(chrom="chr1", start=0, end=22, name=name, reads=reads, reads_rpm=10.0,
                major_rna="ACGU", dicer_call=22, size_purity=0.5, frac_top=0.4, pct_1u=0.7,
                pct_10a=0.3, pingpong_z=1.5, phase_pvalue=0.2, multimap_frac=0.1,
                unique_frac=0.9, cls="miRNA", confidence=conf, review="True", evidence="x")
    base.update(kw)
    return base


def _cl(n, reads, conf):
    return rd.Cluster(name=n, chrom="c", start=0, end=1, reads=reads, reads_rpm=None, cls="x",
                      confidence=conf, dicer_call=None, size_purity=None, frac_top=None,
                      pct_1u=None, pct_10a=None, pingpong_z=None, phase_pvalue=None,
                      multimap_frac=None, unique_frac=None, review=False, major_rna="",
                      evidence="", e={}, gates={}, has_evidence=False)


def test_load_tolerates_old_columns(tmp_path):
    cl = rd.load_clusters(_write(tmp_path / "a.tsv", [_row()]))
    assert len(cl) == 1
    assert cl[0].has_evidence is False
    assert set(cl[0].e) == {"1u", "5a", "fractop_bid", "pingpong", "phase"}
    assert all(v is None for v in cl[0].e.values())
    assert cl[0].gates == {}
    assert cl[0].review is True and cl[0].reads == 500.0 and cl[0].dicer_call == 22


def test_na_and_N_become_none(tmp_path):
    cl = rd.load_clusters(_write(tmp_path / "a.tsv", [_row(dicer_call="N", score="NA", frac_top="")]))
    assert cl[0].dicer_call is None and cl[0].frac_top is None


def test_evidence_and_gates(tmp_path):
    cols = OLD + ["e_1u", "e_5a", "e_fractop_bid", "e_pingpong", "e_phase", "gate_miRNA"]
    row = _row(e_1u=0.1, e_5a=0.2, e_fractop_bid=0.3, e_pingpong=0.4, e_phase="NA", gate_miRNA=1)
    cl = rd.load_clusters(_write(tmp_path / "a.tsv", [row], cols))[0]
    assert cl.has_evidence and cl.e["5a"] == 0.2 and cl.e["phase"] is None
    assert cl.gates == {"miRNA": True}


def test_select_floor():
    cls = [_cl("a", 299, .5), _cl("b", 300, .5), _cl("c", 301, .5)]
    s = rd.select(cls, 300)
    assert [c.name for c in s.shown] == ["c", "b"]
    assert s.n_hidden == 1 and s.n_total == 3 and s.floor == 300 and not s.capped


def test_select_sorted_by_confidence_then_reads():
    cls = [_cl("a", 400, .5), _cl("b", 900, .5), _cl("c", 350, .9), _cl("d", 500, None)]
    assert [c.name for c in rd.select(cls, 0).shown] == ["c", "b", "a", "d"]


def test_select_cap():
    s = rd.select([_cl(str(i), 1000, .5) for i in range(3001)], 0, max_detail=3000)
    assert s.capped and len(s.detail) == 3000 and len(s.shown) == 3001


def test_floor_above_all():
    s = rd.select([_cl("a", 10, .5), _cl("b", 20, .5)], 1000)
    assert s.shown == [] and s.n_hidden == s.n_total == 2


def test_library_summary_without_fragments():
    a, b = _cl("a", 100, .5), _cl("b", 300, .5)
    a.reads_rpm, b.reads_rpm = 100.0, 300.0
    a.unique_frac, b.unique_frac = 1.0, 0.0
    s = rd.library_summary([a, b], None)
    assert s["fragments"] is None and s["fragment_reads"] is None
    assert s["library_reads"] == 1_000_000 and s["cluster_reads"] == 400
    assert abs(s["unique_share"] - 0.25) < 1e-9


def test_library_summary_with_fragments(tmp_path):
    p = tmp_path / "f.tsv"
    p.write_text("chrom\tstart\tend\treads\tread_span\treason\nc\t0\t9\t5.5\t9\tx\nc\t0\t9\t4.5\t9\tx\n")
    s = rd.library_summary([], str(p))
    assert s["fragments"] == 2 and s["fragment_reads"] == 10.0 and s["library_reads"] is None


# ---- final review fixes ----

def test_start_end_are_ints_and_none_safe(tmp_path):
    cl = rd.load_clusters(_write(tmp_path / "a.tsv", [_row(start="100.0", end="250.0"),
                                                      _row(name="c2", start="NA", end="")]))
    assert cl[0].start == 100 and isinstance(cl[0].start, int)
    assert cl[0].end == 250 and isinstance(cl[0].end, int)
    assert cl[1].start is None and cl[1].end is None


def test_detect_kingdom_from_gates_or_classes():
    a = _cl("a", 500, .9)
    assert rd.detect_kingdom([a]) == "animal"
    b = _cl("b", 500, .9)
    b.gates = {"hcsiRNA": False, "siRNA": True}
    assert rd.detect_kingdom([a, b]) == "plant"
    c = _cl("c", 500, .9)
    c.cls = "phasiRNA"
    assert rd.detect_kingdom([c]) == "plant"
    assert rd.select([a, b], 0).kingdom == "plant"
    assert rd.select([a], 0).kingdom == "animal"


def test_detect_mode():
    a = _cl("a", 500, .9)
    assert rd.detect_mode([a]) == "rules"
    a.evidence = "ml:hgb"
    assert rd.detect_mode([a]) == "ml"


def test_select_puts_ambiguous_after_resolved():
    amb = _cl("amb", 900, .99)
    amb.cls = "ambiguous"
    ok = _cl("ok", 400, .6)
    s = rd.select([amb, ok], 0, max_detail=1)
    assert [c.name for c in s.shown] == ["ok", "amb"]
    assert [c.name for c in s.detail] == ["ok"]


def test_class_specs_follow_the_rules_files():
    a = rd.class_specs("animal")
    assert a["piRNA"]["defining"] == {"1u", "pingpong", "size"}
    assert a["siRNA"]["defining"] == {"strand", "size"}
    assert a["piRNA"]["size"] == (24, 31) and a["piRNA"]["min_length"] == 1000
    p = rd.class_specs("plant")
    assert p["hcsiRNA"]["defining"] == {"5a", "size"}
    assert p["phasiRNA"]["defining"] == {"phasing", "size"}
    assert p["phasiRNA24"]["phase_max_p"] == 0.01
    assert "piRNA" not in p
