from srnaprofiler import build_parser


def test_classplots_subcommand_parses():
    parser = build_parser()
    args = parser.parse_args([
        "classplots", "--clusters", "c.tsv", "--bam", "x.bam", "-o", "out",
        "--prefix", "run", "--top-n", "3", "--region", "chr1:100-400",
        "--region", "chr2:5-9", "--no-summary"])
    assert args.clusters == "c.tsv"
    assert args.bam == "x.bam"
    assert args.outdir == "out"
    assert args.top_n == 3
    assert args.region == ["chr1:100-400", "chr2:5-9"]
    assert args.no_summary is True
    assert args.func.__name__ == "run_classplots"
