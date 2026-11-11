"""Plain-language explanations for the HTML report.

DISPLAY_PRESENT and DISPLAY_WEAK only decide how an evidence bar is worded in the
report (present / weak / absent). They are display-only and are NOT the classifier's
thresholds.
"""
DISPLAY_PRESENT = 0.5
DISPLAY_WEAK = 0.25
# Display-only cutoffs for pct_10a (raw fraction of reads with A at position 10).
DISPLAY_10A_PRESENT = 0.5
DISPLAY_10A_WEAK = 0.25
DISPLAY_SIZE_PRESENT = 0.5
DISPLAY_SIZE_WEAK = 0.3
_MIN_CLUSTERS = 20
_LOW_UNIQUE = 0.3

STRINGS = {
    "verdict.piRNA": "This region looks like a piRNA cluster: small RNAs of the length typical of animal germline defence, usually starting with U and amplified in pairs.",
    "verdict.siRNA": "This region looks like siRNA-producing: small RNAs of one tight size cut by the Dicer enzyme from double-stranded RNA.",
    "verdict.miRNA": "This region looks like a microRNA locus: a precise size and a one-sided pile of reads, as expected from a hairpin precursor.",
    "verdict.phasiRNA": "This region looks like phasiRNA-producing: small RNAs released in regular steps along one transcript.",
    "verdict.phasiRNA24": "This region looks like 24-nt phasiRNA-producing: 24-nt small RNAs released in regular steps along one transcript.",
    "verdict.hcsiRNA": "This region looks like heterochromatic siRNA-producing: 24-nt small RNAs typical of plant DNA silencing.",
    "verdict.ambiguous": "The evidence does not point to a single class, so this region was left unresolved.",
    "verdict.other": "This region has small RNAs, but they do not match any class the tool knows.",
    "verdict.generic": "This region was assigned to the class named above; no class-specific explanation is available.",
    "status.present": "present",
    "status.weak": "weak",
    "status.absent": "absent",
    "status.unavailable": "not available",
    "check.size.title": "Size consistency",
    "check.size.detail": "Share of reads with the dominant length: {v}. Dicer-made RNAs share one length.",
    "check.1u.title": "Starts with U",
    "check.1u.detail": "Evidence score {v}. Many piRNAs and miRNAs start with the base U.",
    "check.10a.title": "A at position 10",
    "check.10a.detail": "Share of reads with A at position 10: {v}. Partner piRNAs in ping-pong amplification often carry an A at position 10.",
    "check.strand.title": "Strand balance",
    "check.strand.detail": "Evidence score {v}. Hairpin and piRNA loci are often one-sided; Dicer products of long duplexes use both strands.",
    "check.pingpong.title": "Ping-pong pairing",
    "check.pingpong.detail": "Evidence score {v}. Reads on opposite strands overlap by 10 bases when piRNAs amplify each other.",
    "check.phasing.title": "Regular phasing",
    "check.phasing.detail": "Evidence score {v}. Phased RNAs are cut at regular steps along one transcript.",
    "check.multimap.title": "Mapping uniqueness",
    "check.multimap.detail": "Share of reads mapping to several places: {v}. High values make the call less certain.",
    "check.na": "This table has no such value.",
    "whynot.default": "The checks for {cls} were not met.",
    "whynot.unknown": "Per-class checks are not available in this table.",
    "caution.noevidence": "This table predates the evidence columns, so some checks cannot be shown.",
    "caution.review": "The tool flagged this region for manual review.",
    "caution.multimap": "Many reads map to several places, so counts here may be shared with other regions.",
    "caution.lowconf": "Confidence is low; treat the class as a suggestion.",
    "missing.1u": "no clear U at the first base",
    "missing.10a": "no A at position 10",
    "missing.pingpong": "no ping-pong pairing",
    "missing.phasing": "no regular phasing",
    "summary.empty": "No clusters passed the depth filter (at least {floor} reads).",
    "summary.main": "{n} clusters passed the depth filter (at least {floor} reads); the most common class is {top} ({k}).",
    "warn.few": "Fewer than {min} clusters passed the depth filter: only {n} (at least {floor} reads); {hidden} were hidden. Lower the minimum depth to see more.",
    "warn.cap": "Only the top {n} clusters are shown in detail; the rest appear in the table only.",
    "warn.unique": "Only {pct}% of reads map to a single place; many small RNAs here may come from repeats, so classes are less certain.",
    "ui.brand": "srnaprofiler",
    "ui.title.default": "srnaprofiler report",
    "ui.na": "n/a",
    "ui.unclassified": "unclassified",
    "ui.howto": "How to read: ",
    "ui.tab.overview": "Overview",
    "ui.tab.clusters": "Clusters",
    "ui.tab.glossary": "Glossary",
    "ui.tab.references": "References",
    "ui.dot.size": "Size",
    "ui.dot.1u": "Starts with U",
    "ui.dot.10a": "A at position 10",
    "ui.dot.pingpong": "Ping-pong",
    "ui.dot.phasing": "Phasing",
    "ui.donut.empty": "No clusters to draw.",
    "ui.sizes.empty": "No dominant read size is available.",
    "ui.ov.title": "What is in this library?",
    "ui.ov.sub": "A plain-language summary before any technical detail.",
    "ui.ov.main": "Main numbers",
    "ui.ov.kinds": "Which kinds of small RNA?",
    "ui.ov.library": "The library",
    "ui.stat.analysed": "Clusters analysed",
    "ui.stat.analysed.sub": "regions with many reads",
    "ui.stat.hidden": "Hidden by the depth filter",
    "ui.stat.hidden.sub": "of {total} (below {floor} reads)",
    "ui.stat.amb": "Ambiguous",
    "ui.stat.amb.sub": "not enough evidence for one class",
    "ui.stat.reads": "Reads in analysed clusters",
    "ui.card.classes": "Classes found",
    "ui.help.classes": "each slice is the share of clusters in one class. Ambiguous is not an error: the evidence did not point to a single class.",
    "ui.card.sizes": "Dominant read size (nt)",
    "ui.help.sizes": "each cluster is counted once, at its most common read length, coloured by class. About 21 to 22 nt suggests siRNA or miRNA; about 24 to 30 nt suggests piRNA.",
    "ui.lib.size": "Library size",
    "ui.lib.size.sub": "reads, estimated from RPM",
    "ui.lib.unique": "Reads mapping to one place",
    "ui.lib.unique.sub": "weighted by cluster reads",
    "ui.lib.reads": "Reads in all clusters",
    "ui.lib.reads.sub": "before the depth filter",
    "ui.lib.frag": "Fragments below the pipeline floor",
    "ui.lib.frag.sub": "{n} reads",
    "ui.list.title": "Which clusters deserve attention?",
    "ui.list.sub": "Ordered by confidence by default. Click a row to see why.",
    "ui.all": "All",
    "ui.sort.conf": "Sort: confidence",
    "ui.sort.reads": "Sort: reads",
    "ui.sort.name": "Sort: name",
    "ui.floor.label": "Minimum reads ",
    "ui.floor.note": "This box can only narrow the table to the {n} clusters already in this file (at least {floor} reads). To see clusters with fewer reads, rerun the report with --min-reads set lower.",
    "ui.shown": "{n} shown (of {m} embedded)",
    "ui.list.none.filter": "No clusters match these filters.",
    "ui.list.none.rows": "No clusters passed the depth filter, so there is nothing to list.",
    "ui.row.missing": "Missing: {list}",
    "ui.row.missing.nodetail": "Missing: {list} (no detail page, cap)",
    "ui.row.nodetail": "no detail page (cap)",
    "ui.th.cluster": "Cluster",
    "ui.th.class": "Class",
    "ui.th.conf": "Confidence",
    "ui.th.evidence": "Evidence",
    "ui.th.reads": "Reads",
    "ui.page.prev": "Previous",
    "ui.page.next": "Next",
    "ui.page.info": "Page {p} of {pages} ({n} clusters)",
    "ui.dots.legend": "Evidence dots, in order: size, starts with U, A at position 10, ping-pong, phasing.",
    "ui.strand.plus": "+ strand",
    "ui.strand.minus": "- strand",
    "ui.pairs": "pairs",
    "ui.fig.na": "Figures are not available for this cluster.",
    "ui.fig.sizes": "Read sizes by strand",
    "ui.fig.sizes.help": "each bar is the number of reads of one length. Tight clusters have one tall bar.",
    "ui.fig.nt": "Nucleotide at each position",
    "ui.fig.nt.help": "reads of the dominant size, base by base. A tall U at position 1 and A at position 10 is the 1U/10A signature.",
    "ui.fig.cov": "Coverage along the cluster",
    "ui.fig.cov.help": "reads per bin along the region. Above the line is the + strand, below is the - strand.",
    "ui.fig.pp": "Overlap between strands (ping-pong)",
    "ui.fig.pp.help": "how often reads on opposite strands overlap by 1 to {n} bases. A peak at 10 means ping-pong.",
    "ui.detail.back": "Back to the list",
    "ui.detail.sub": "{where} | {reads} reads",
    "ui.detail.size": " | dominant size {z} nt",
    "ui.detail.conf": "confidence, from 0 to 1; it is not a probability of being right",
    "ui.detail.why": "Why this class?",
    "ui.check.line": "{title}: {status}",
    "ui.detail.figs": "The signals, in figures",
    "ui.detail.whynot": "Why not another class?",
    "ui.detail.cautions": "Cautions",
    "ui.detail.nocautions": "No specific cautions for this cluster.",
    "ui.detail.seq": "Dominant sequence",
    "ui.detail.noseq": "No dominant sequence is recorded for this cluster.",
    "ui.gl.title": "Glossary",
    "ui.gl.sub": "For readers new to small RNAs. Numbers in brackets point to the References tab.",
    "ui.ref.title": "References",
    "ui.footer": "The coloured evidence dots use display-only cutoffs, not the classifier's thresholds. 1U, ping-pong and phasing scores: present at {present} or more, weak at {weak} or more. A at position 10: present at {p10} or more, weak at {w10} or more. Size consistency: present at {psize} or more, weak at {wsize} or more. They only decide how a bar is worded in this report.",
}

REFERENCES = {
    "brennecke2007": "Brennecke, J., Aravin, A. A., Stark, A., Dus, M., Kellis, M., Sachidanandam, R., & Hannon, G. J. (2007). Discrete Small RNA-Generating Loci as Master Regulators of Transposon Activity in Drosophila. *Cell*, 128(6), 1089-1103. https://doi.org/10.1016/j.cell.2007.01.043",
    "czech2016": "Czech, B., & Hannon, G. J. (2016). One Loop to Rule Them All: The Ping-Pong Cycle and piRNA-Guided Silencing. *Trends in Biochemical Sciences*, 41(4), 324-337. https://doi.org/10.1016/j.tibs.2015.12.008",
    "bernstein2001": "Bernstein, E., Caudy, A. A., Hammond, S. M., & Hannon, G. J. (2001). Role for a bidentate ribonuclease in the initiation step of RNA interference. *Nature*, 409(6818), 363-366. https://doi.org/10.1038/35053110",
    "hamilton1999": "Hamilton, A. J., & Baulcombe, D. C. (1999). A Species of Small Antisense RNA in Posttranscriptional Gene Silencing in Plants. *Science*, 286(5441), 950-952. https://doi.org/10.1126/science.286.5441.950",
    "bartel2004": "Bartel, D. P. (2004). MicroRNAs. *Cell*, 116(2), 281-297. https://doi.org/10.1016/s0092-8674(04)00045-5",
    "allen2005": "Allen, E., Xie, Z., Gustafson, A. M., & Carrington, J. C. (2005). microRNA-Directed Phasing during Trans-Acting siRNA Biogenesis in Plants. *Cell*, 121(2), 207-221. https://doi.org/10.1016/j.cell.2005.04.004",
    "axtell2013": "Axtell, M. J. (2013). Classification and Comparison of Small RNAs from Plants. *Annual Review of Plant Biology*, 64(1), 137-159. https://doi.org/10.1146/annurev-arplant-050312-120043",
}

GLOSSARY = [
    {"term": "small RNA", "text": "A short RNA molecule, roughly 20 to 30 bases long, that helps switch genes or mobile DNA on or off.", "cite": ["axtell2013"]},
    {"term": "cluster", "text": "A stretch of the genome where many small RNA reads pile up. Each cluster is classified as a unit.", "cite": []},
    {"term": "piRNA", "text": "A small RNA, found mostly in animal germlines, that silences mobile DNA and is produced from discrete genomic loci.", "cite": ["brennecke2007"]},
    {"term": "siRNA", "text": "A small RNA of one tight size cut by the Dicer enzyme from double-stranded RNA.", "cite": ["bernstein2001", "hamilton1999"]},
    {"term": "miRNA", "text": "A small RNA of about 21 to 22 bases cut from a hairpin-shaped precursor, which tunes gene activity.", "cite": ["bartel2004"]},
    {"term": "phasiRNA", "text": "Plant small RNAs released in regular steps from one transcript, usually after a microRNA cuts it.", "cite": ["allen2005", "axtell2013"]},
    {"term": "1U/10A", "text": "1U means the RNA starts with U; 10A means it has an A at position 10. Both are typical of piRNAs.", "cite": ["brennecke2007"]},
    {"term": "ping-pong", "text": "A loop in which two piRNAs on opposite strands make each other, leaving pairs that overlap by 10 bases.", "cite": ["brennecke2007", "czech2016"]},
    {"term": "phasing", "text": "Cuts at regular steps along a transcript, so reads fall into a repeating pattern.", "cite": ["allen2005"]},
    {"term": "ambiguous", "text": "The evidence fits more than one class, or none, so the tool does not choose.", "cite": []},
    {"term": "confidence", "text": "How strongly the evidence supports the assigned class, from 0 to 1. It is not a probability of being right.", "cite": []},
    {"term": "multimapping", "text": "A read that fits several places in the genome, so its true origin is uncertain.", "cite": []},
]

_KNOWN = {"piRNA", "siRNA", "miRNA", "phasiRNA", "phasiRNA24", "hcsiRNA", "ambiguous", "other"}


def _level(v, present=DISPLAY_PRESENT, weak=DISPLAY_WEAK):
    if v is None:
        return "unavailable"
    return "present" if v >= present else "weak" if v >= weak else "absent"


def evidence_status(c):
    """Display-only status of each piece of evidence for a cluster."""
    e = c.e or {}
    return {
        "size": _level(c.size_purity, DISPLAY_SIZE_PRESENT, DISPLAY_SIZE_WEAK),
        "1u": _level(e.get("1u")),
        "10a": _level(c.pct_10a, DISPLAY_10A_PRESENT, DISPLAY_10A_WEAK),
        "pingpong": _level(e.get("pingpong")),
        "phasing": _level(e.get("phase")),
    }


def _fmt(v):
    return "n/a" if v is None else f"{v:.2f}"


def explain(c):
    st = evidence_status(c)
    e = c.e or {}
    S = STRINGS
    strand_v = e.get("fractop_bid")
    mm = c.multimap_frac
    rows = [
        ("size", st["size"], c.size_purity),
        ("1u", st["1u"], e.get("1u")),
        ("10a", st["10a"], c.pct_10a),
        ("strand", _level(strand_v), strand_v),
        ("pingpong", st["pingpong"], e.get("pingpong")),
        ("phasing", st["phasing"], e.get("phase")),
        ("multimap", "unavailable" if mm is None else
         ("absent" if mm < 0.3 else "weak" if mm < 0.6 else "present"), mm),
    ]
    checks = [{"key": k, "status": s, "title": S[f"check.{k}.title"],
               "detail": S[f"check.{k}.detail"].format(v=_fmt(v)) if v is not None else S["check.na"]}
              for k, s, v in rows]
    verdict = S.get(f"verdict.{c.cls}", S["verdict.generic"])
    if c.gates:
        why_not = [{"cls": k, "reason": S["whynot.default"].format(cls=k)}
                   for k, ok in c.gates.items() if not ok and k != c.cls]
    else:
        why_not = [{"cls": "", "reason": S["whynot.unknown"]}]
    cautions = []
    if not c.has_evidence:
        cautions.append(S["caution.noevidence"])
    if c.review:
        cautions.append(S["caution.review"])
    if mm is not None and mm >= 0.6:
        cautions.append(S["caution.multimap"])
    if c.confidence is not None and c.confidence < 0.3:
        cautions.append(S["caution.lowconf"])
    missing = []
    if c.cls == "ambiguous":
        for key, name in (("1u", "1u"), ("10a", "10a"), ("pingpong", "pingpong"), ("phasing", "phasing")):
            if st[key] in ("absent", "unavailable"):
                missing.append(S[f"missing.{name}"])
    return {"verdict": verdict, "checks": checks, "why_not": why_not,
            "cautions": cautions, "missing": missing}


def summary_sentence(counts, floor):
    S = STRINGS
    f = f"{floor:g}"
    total = sum(counts.values())
    if not total:
        return S["summary.empty"].format(floor=f)
    top, k = max(counts.items(), key=lambda kv: kv[1])
    return S["summary.main"].format(n=total, floor=f, top=top, k=k)


def warnings(sel, summary):
    S = STRINGS
    out = []
    if len(sel.shown) < _MIN_CLUSTERS:
        out.append(S["warn.few"].format(min=_MIN_CLUSTERS, n=len(sel.shown), floor=f"{sel.floor:g}", hidden=sel.n_hidden))
    if sel.capped:
        out.append(S["warn.cap"].format(n=len(sel.detail)))
    u = (summary or {}).get("unique_share")
    if u is not None and u < _LOW_UNIQUE:
        out.append(S["warn.unique"].format(pct=round(u * 100)))
    return out
