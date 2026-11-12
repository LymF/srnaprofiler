import pysam


def make_bam(path, reads):
    """reads: list of (chrom, start, secondary: bool). 20-bp reads, coordinate sorted + indexed."""
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": n, "LN": 5000} for n in ("c1", "c2", "chr1")]}
    unsorted = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for i, (chrom, start, secondary) in enumerate(reads):
            a = pysam.AlignedSegment(out.header)
            a.query_name = f"r{i}"
            a.query_sequence = "A" * 20
            a.query_qualities = pysam.qualitystring_to_array("I" * 20)
            a.reference_id = out.get_tid(chrom)
            a.reference_start = start
            a.cigar = [(0, 20)]
            a.flag = 256 if secondary else 0
            out.write(a)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))


def make_seq_bam(path, reads, chroms=(("c1", 5000),)):
    """reads: list of (chrom, start, seq, strand, nh[, flag]). `seq` in reference
    orientation (as stored in BAM); nh=None writes no NH tag. Sorted + indexed."""
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": n, "LN": ln} for n, ln in chroms]}
    unsorted = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for i, r in enumerate(reads):
            chrom, start, seq, strand, nh = r[:5]
            a = pysam.AlignedSegment(out.header)
            a.query_name = f"r{i}"
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.reference_id = out.get_tid(chrom)
            a.reference_start = start
            a.cigar = [(0, len(seq))]
            a.flag = (16 if strand == "-" else 0) | (r[5] if len(r) > 5 else 0)
            if nh is not None:
                a.set_tag("NH", nh)
            out.write(a)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))
