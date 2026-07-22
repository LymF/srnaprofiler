from qc import prepare_genome

def test_prepare_genome_does_not_tile(tmp_path):
    g = tmp_path / "g.fa"
    g.write_text(">chr1\n" + "ACGT" * 1000 + "\n")   # 4000 nt, 1 sequência
    out = prepare_genome(str(g), str(tmp_path), "t")
    headers = [l for l in open(out) if l.startswith(">")]
    assert len(headers) == 1                      # continua 1 sequência (não tileou em janelas)
    assert headers[0].strip() == ">chr1"
