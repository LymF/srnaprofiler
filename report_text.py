"""Plain-language explanations for the HTML report.

DISPLAY_* values only decide how an evidence bar is worded in the report
(present / weak / absent, low / moderate / high). They are display-only and are
NOT the classifier's thresholds.
"""
import re

import report_data

DISPLAY_PRESENT = 0.5
DISPLAY_WEAK = 0.25
# Display-only cutoffs for pct_10a (raw fraction of reads with A at position 10).
DISPLAY_10A_PRESENT = 0.5
DISPLAY_10A_WEAK = 0.25
DISPLAY_SIZE_PRESENT = 0.5
DISPLAY_SIZE_WEAK = 0.3
# Display-only cutoffs for multimap_frac (share of reads mapping to several places).
DISPLAY_MULTIMAP_MODERATE = 0.3
DISPLAY_MULTIMAP_HIGH = 0.6
_MIN_CLUSTERS = 20
_LOW_UNIQUE = 0.3

STRINGS = {
    "verdict.piRNA": "This region looks like a piRNA cluster: many reads of 24 to 31 nt, most starting with U or paired across the strands by ping-pong. Animals make piRNAs this way in the germline to silence transposons.",
    "verdict.siRNA": "This region looks like siRNA-producing: reads of about 21 nt from both strands, as Dicer leaves when it cuts long double-stranded RNA.",
    "verdict.miRNA": "This region looks like a microRNA locus: a precise size and a one-sided pile of reads, as expected from a hairpin precursor.",
    "verdict.phasiRNA": "This region looks like phasiRNA-producing: small RNAs released in regular steps along one transcript.",
    "verdict.phasiRNA24": "This region looks like 24-nt phasiRNA-producing: 24-nt small RNAs released in regular steps along one transcript.",
    "verdict.hcsiRNA": "This region looks like heterochromatic siRNA-producing: reads of 23 to 24 nt, many starting with A, typical of the plant small RNAs that guide DNA silencing at repeats.",
    "verdict.ambiguous": "This region was left unresolved: {reason}",
    "verdict.other": "This region has small RNAs, but they do not match any class the tool knows.",
    "verdict.generic": "This region was assigned to the class named above; no class-specific explanation is available.",
    "status.present": "present",
    "status.weak": "weak",
    "status.absent": "absent",
    "status.unavailable": "not available",
    "status.low": "low",
    "status.moderate": "moderate",
    "status.high": "high",
    "check.size.title": "Size consistency",
    "check.size.detail": "Share of reads within 1 nt of the dominant length: {v}. Small RNAs of one class share a narrow range of lengths.",
    "check.1u.title": "Starts with U",
    "check.1u.detail": "Evidence score {v}. Most piRNAs and miRNAs start with the base U.",
    "check.1u.detail.plant": "Evidence score {v}. Plant miRNAs and many 21-nt siRNAs start with U, which sorts them into the AGO1 protein.",
    "check.10a.title": "A at position 10",
    "check.5a.title": "Starts with A",
    "check.5a.detail": "Evidence score {v}. Most plant 24-nt heterochromatic siRNAs start with A, which sorts them into the AGO4 protein.",
    "check.10a.detail": "Share of reads with A at position 10: {v}. Partner piRNAs in ping-pong amplification often carry an A at position 10.",
    "check.strand.title": "Strand balance",
    "check.strand.detail": "Evidence score {v}. Hairpin and many piRNA loci are one-sided; Dicer products of long double-stranded RNA come from both strands.",
    "check.strand.detail.plant": "Evidence score {v}. Hairpin (miRNA) loci are one-sided; Dicer products of long double-stranded RNA come from both strands.",
    "check.pingpong.title": "Ping-pong pairing",
    "check.pingpong.detail": "Evidence score {v}. Reads on opposite strands overlap by 10 bases when piRNAs amplify each other.",
    "check.phasing.title": "Regular phasing",
    "check.phasing.detail": "Evidence score {v}. Phased RNAs are cut at regular steps along one transcript.",
    "check.multimap.title": "Multimapping",
    "check.multimap.detail": "Share of reads mapping to several places: {v}. Low is better; high values make the call less certain.",
    "check.context": "not needed for this class",
    "check.note.none": "None of these checks decides this class on its own; they are shown as context.",
    "check.note.miRNA": "miRNA calls come from the hairpin search (miRDeep2), not from these checks; they are shown as context.",
    "check.na": "This table has no such value.",
    "whynot.default": "The basic checks for {cls} were not met.",
    "whynot.nosize": "No dominant read size could be called, so {cls} could not be tested.",
    "whynot.size": "The dominant size ({z} nt) is outside the {cls} window ({lo} to {hi} nt).",
    "whynot.length": "The region spans {span} bp, less than the {min} bp a {cls} cluster needs.",
    "whynot.phasep": "The phasing test is not significant (p = {p}; {cls} needs {lim} or less).",
    "whynot.phasing": "No regular phasing, which {cls} requires.",
    "whynot.none": "No other class failed its basic checks (size window, cluster length, phasing test); the call above came from comparing the strength of the evidence.",
    "whynot.unknown": "Per-class checks are not available in this table.",
    "caution.noevidence": "This table predates the evidence columns, so some checks cannot be shown.",
    "caution.review": "The tool flagged this region for review: the evidence for this class was not clearly stronger than for another class, or the reads were few. Treat the class as a suggestion.",
    "caution.multimap": "Many reads map to several places, so counts here may be shared with other regions.",
    "amb.size": "the dominant read size ({z} nt) fits no class's size window; the closest class is {cls}.",
    "amb.nosize": "no dominant read size could be called, so no class could be tested.",
    "amb.criterion": "it looks most like {cls}, but it fails the {cls} cluster criterion: {list}.",
    "amb.crit.length": "the reads span less than the minimum cluster length",
    "amb.crit.density": "the reads are not denser than expected by chance",
    "amb.ml": "the trained model did not pick a class: its top answer was 'other', no class it favoured passed that class's basic checks, or its score was below the cutoff.",
    "amb.generic": "it did not meet the requirements of any class.",
    "summary.empty": "No clusters passed the depth filter (at least {floor} reads).",
    "summary.main": "{n} clusters passed the depth filter (at least {floor} reads); the most common class is {top} ({k}).",
    "warn.few": "Fewer than {min} clusters passed the depth filter: only {n} (at least {floor} reads); {hidden} were hidden. To see more, rerun the report with a lower --min-reads value.",
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
    "ui.dot.5a": "Starts with A",
    "ui.dot.strand": "Both strands",
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
    "ui.stat.amb.sub": "no class's requirements met",
    "ui.stat.reads": "Reads in analysed clusters",
    "ui.card.classes": "Classes found",
    "ui.help.classes": "each slice is the share of clusters in one class. Ambiguous is not an error: the region did not meet the requirements of any one class.",
    "ui.card.sizes": "Dominant read size (nt)",
    "ui.help.sizes.animal": "each cluster is counted once, at its most common read length, coloured by class. About 20 to 22 nt suggests siRNA or miRNA; about 24 to 31 nt suggests piRNA.",
    "ui.help.sizes.plant": "each cluster is counted once, at its most common read length, coloured by class. About 21 to 22 nt suggests miRNA, phasiRNA or siRNA; 23 to 24 nt suggests heterochromatic siRNA (or, in monocots such as rice and maize, 24-nt phasiRNA).",
    "ui.lib.size": "Library size",
    "ui.lib.size.sub": "reads, estimated from RPM",
    "ui.lib.unique": "Reads mapping to one place",
    "ui.lib.unique.sub": "weighted by cluster reads",
    "ui.lib.reads": "Reads in all clusters",
    "ui.lib.reads.sub": "before the depth filter",
    "ui.lib.frag": "Fragments below the pipeline floor",
    "ui.lib.frag.sub": "{n} reads",
    "ui.list.title": "Which clusters deserve attention?",
    "ui.list.sub": "Resolved clusters first, by confidence; ambiguous ones after them. Click a row to see why.",
    "ui.all": "All",
    "ui.sort.conf": "Sort: confidence",
    "ui.sort.reads": "Sort: reads",
    "ui.sort.name": "Sort: name",
    "ui.floor.label": "Minimum reads ",
    "ui.floor.note": "This box can only narrow the table to the {n} clusters already in this file (at least {floor} reads). To see clusters with fewer reads, rerun the report with --min-reads set lower.",
    "ui.shown": "{n} shown (of {m} embedded)",
    "ui.list.none.filter": "No clusters match these filters.",
    "ui.list.none.rows": "No clusters passed the depth filter, so there is nothing to list.",
    "ui.row.unresolved": "Unresolved: {reason}",
    "ui.row.unresolved.nodetail": "Unresolved: {reason} (no detail page, cap)",
    "ui.row.nodetail": "no detail page (cap)",
    "ui.th.cluster": "Cluster",
    "ui.th.class": "Class",
    "ui.th.conf": "Confidence",
    "ui.th.evidence": "Evidence",
    "ui.th.reads": "Reads",
    "ui.conf.amb": "n/a (nearest class: {cls})",
    "ui.conf.amb.none": "n/a",
    "ui.where": "{chrom}:{start}-{end}",
    "ui.page.prev": "Previous",
    "ui.page.next": "Next",
    "ui.page.info": "Page {p} of {pages} ({n} clusters)",
    "ui.dots.legend": "Evidence dots, in order: {list}.",
    "ui.strand.plus": "+ strand",
    "ui.strand.minus": "- strand",
    "ui.pairs": "pairs",
    "ui.fig.na": "Figures are not available for this cluster.",
    "ui.fig.noreads": "No reads were found in this region of the BAM, so there is nothing to draw.",
    "ui.fig.err": "This region could not be read from the BAM (see the warning in the log), so there is nothing to draw.",
    "ui.fig.sizes": "Read sizes by strand",
    "ui.fig.sizes.help": "each bar is the number of reads of one length. Tight clusters have one tall bar.",
    "ui.fig.nt": "Nucleotide at each position",
    "ui.fig.nt.help": "reads of the dominant size, base by base. A tall U at position 1 and A at position 10 is the 1U/10A signature of piRNAs.",
    "ui.fig.nt.help.plant": "reads of the dominant size, base by base. The first base shows which Argonaute protein the RNAs are sorted into: U for AGO1 (miRNAs), A for AGO4 (24-nt siRNAs).",
    "ui.fig.cov": "Coverage along the cluster",
    "ui.fig.cov.help": "mean number of reads covering each base, per bin along the region, not weighted by multimapping. Above the line is the + strand, below is the - strand.",
    "ui.fig.pp": "Overlap between strands (ping-pong)",
    "ui.fig.pp.help": "how often reads on opposite strands overlap by 1 to {n} bases. A peak at 10 means ping-pong.",
    "ui.fig.pp.plant": "Overlap between strands",
    "ui.fig.pp.help.plant": "how often reads on opposite strands overlap by 1 to {n} bases.",
    "ui.detail.back": "Back to the list",
    "ui.detail.sub": "{where} | {reads} reads",
    "ui.detail.size": " | dominant size {z} nt",
    "ui.detail.conf": "confidence: how strongly the evidence favours this class over the others ({k} classes compared). The lowest possible value is {min}, 1 divided by {k}, which means no preference; it is not a probability of being right.",
    "ui.detail.conf.ml": "confidence: the trained model's score for its top class, from 0 to 1; it is not a probability of being right.",
    "ui.detail.conf.amb": "No confidence is shown for an unresolved region: the number the tool stores belongs to the nearest class, which was not assigned.",
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
    "ui.footer": "The coloured evidence dots use display-only cutoffs, not the classifier's thresholds. 1U, 5'A, strand, ping-pong and phasing scores: present at {present} or more, weak at {weak} or more. A at position 10: present at {p10} or more, weak at {w10} or more. Size consistency: present at {psize} or more, weak at {wsize} or more. Multimapping: moderate at {mmod} or more, high at {mhigh} or more. They only decide how a bar is worded in this report.",
}

REFERENCES = {
    "brennecke2007": "Brennecke, J., Aravin, A. A., Stark, A., Dus, M., Kellis, M., Sachidanandam, R., & Hannon, G. J. (2007). Discrete Small RNA-Generating Loci as Master Regulators of Transposon Activity in Drosophila. *Cell*, 128(6), 1089-1103. https://doi.org/10.1016/j.cell.2007.01.043",
    "czech2016": "Czech, B., & Hannon, G. J. (2016). One Loop to Rule Them All: The Ping-Pong Cycle and piRNA-Guided Silencing. *Trends in Biochemical Sciences*, 41(4), 324-337. https://doi.org/10.1016/j.tibs.2015.12.008",
    "bernstein2001": "Bernstein, E., Caudy, A. A., Hammond, S. M., & Hannon, G. J. (2001). Role for a bidentate ribonuclease in the initiation step of RNA interference. *Nature*, 409(6818), 363-366. https://doi.org/10.1038/35053110",
    "hamilton1999": "Hamilton, A. J., & Baulcombe, D. C. (1999). A Species of Small Antisense RNA in Posttranscriptional Gene Silencing in Plants. *Science*, 286(5441), 950-952. https://doi.org/10.1126/science.286.5441.950",
    "bartel2004": "Bartel, D. P. (2004). MicroRNAs. *Cell*, 116(2), 281-297. https://doi.org/10.1016/s0092-8674(04)00045-5",
    "allen2005": "Allen, E., Xie, Z., Gustafson, A. M., & Carrington, J. C. (2005). microRNA-Directed Phasing during Trans-Acting siRNA Biogenesis in Plants. *Cell*, 121(2), 207-221. https://doi.org/10.1016/j.cell.2005.04.004",
    "mi2008": "Mi, S., Cai, T., Hu, Y., Chen, Y., Hodges, E., Ni, F., Wu, L., Li, S., Zhou, H., Long, C., Chen, S., Hannon, G. J., & Qi, Y. (2008). Sorting of Small RNAs into Arabidopsis Argonaute Complexes Is Directed by the 5' Terminal Nucleotide. *Cell*, 133(1), 116-127. https://doi.org/10.1016/j.cell.2008.02.034",
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
    {"term": "hcsiRNA", "text": "Heterochromatic siRNA: plant small RNAs of about 24 nt, made by Dicer from double-stranded RNA at repeats and transposons, that guide DNA silencing. Most start with A, which sorts them into the AGO4 protein. They are not phased.", "cite": ["axtell2013", "mi2008"]},
    {"term": "ambiguous", "text": "The tool did not assign a class. Either the dominant read size fits no class's size range (or no dominant size could be called), or the region looks most like one class but fails that class's cluster criterion (a minimum length and, for piRNA, a read-density test), or, with the trained model, the model abstained. The region page gives the reason. This differs from the review flag, which marks an assigned class whose evidence was close to another class's.", "cite": []},
    {"term": "confidence", "text": "How strongly the evidence favours the assigned class over the other classes the tool compared. The lowest possible value is 1 divided by the number of classes (0.5 with two classes), which means no preference; values near 1 mean one class clearly wins. With the trained model it is the model's score for its top class. It is not a probability of being right, and it is not shown for ambiguous regions.", "cite": []},
    {"term": "multimapping", "text": "A read that fits several places in the genome, so its true origin is uncertain.", "cite": []},
]

# Report check keys shown per kingdom; the first five are the list's evidence dots.
DOTS = {"animal": ("size", "1u", "10a", "pingpong", "phasing"),
        "plant": ("size", "1u", "5a", "strand", "phasing")}
_CHECKS = {"animal": ("size", "1u", "10a", "strand", "pingpong", "phasing", "multimap"),
           "plant": ("size", "1u", "5a", "strand", "phasing", "multimap")}
# kingdom-specific text: role -> STRINGS key
_KTEXT = {"animal": {"sizes": "ui.help.sizes.animal", "nt": "ui.fig.nt.help",
                     "pp_title": "ui.fig.pp", "pp_help": "ui.fig.pp.help"},
          "plant": {"sizes": "ui.help.sizes.plant", "nt": "ui.fig.nt.help.plant",
                    "pp_title": "ui.fig.pp.plant", "pp_help": "ui.fig.pp.help.plant"}}
_NEAREST = re.compile(r"(?:^|;\s*)nearest=([^;\s]+)")
_CRITERION = re.compile(r"below_cluster_criterion\(([^)]*)\)")


def _level(v, present=DISPLAY_PRESENT, weak=DISPLAY_WEAK):
    if v is None:
        return "unavailable"
    return "present" if v >= present else "weak" if v >= weak else "absent"


def _mm_level(v):
    if v is None:
        return "unavailable"
    return ("high" if v >= DISPLAY_MULTIMAP_HIGH else
            "moderate" if v >= DISPLAY_MULTIMAP_MODERATE else "low")


def _values(c):
    e = c.e or {}
    return {"size": (c.size_purity, DISPLAY_SIZE_PRESENT, DISPLAY_SIZE_WEAK),
            "1u": (e.get("1u"), DISPLAY_PRESENT, DISPLAY_WEAK),
            "5a": (e.get("5a"), DISPLAY_PRESENT, DISPLAY_WEAK),
            "10a": (c.pct_10a, DISPLAY_10A_PRESENT, DISPLAY_10A_WEAK),
            "strand": (e.get("fractop_bid"), DISPLAY_PRESENT, DISPLAY_WEAK),
            "pingpong": (e.get("pingpong"), DISPLAY_PRESENT, DISPLAY_WEAK),
            "phasing": (e.get("phase"), DISPLAY_PRESENT, DISPLAY_WEAK)}


def evidence_status(c, kingdom="animal"):
    """Display-only status of the list's evidence dots for a cluster, in dot order."""
    vals = _values(c)
    return {k: _level(*vals[k]) for k in DOTS[kingdom]}


def dot_legend(kingdom="animal"):
    names = ", ".join(STRINGS[f"ui.dot.{k}"] for k in DOTS[kingdom])
    return STRINGS["ui.dots.legend"].format(list=names)


def size_help(kingdom="animal"):
    return STRINGS[_KTEXT[kingdom]["sizes"]]


def kingdom_text(kingdom="animal"):
    """Kingdom-specific help strings the front end shows (role -> text)."""
    return {role: STRINGS[k] for role, k in _KTEXT[kingdom].items()}


def confidence_floor(kingdom="animal"):
    """Lowest softmax confidence in rules mode: 1 / number of enabled classes."""
    return 1.0 / len(report_data.class_specs(kingdom))


def _fmt(v):
    return "n/a" if v is None else f"{v:.2f}"


def nearest_class(c):
    """Nearest class named in the evidence trail; None without a size call or a trail."""
    if c.dicer_call is None:
        return None
    m = _NEAREST.search(c.evidence or "")
    return m.group(1) if m else None


def ambiguous_reason(c):
    """Why the classifier left the region unresolved, from its evidence trail."""
    S = STRINGS
    ev = c.evidence or ""
    if ev.startswith("ml:"):
        return S["amb.ml"]
    if not _NEAREST.search(ev):
        return S["amb.generic"]
    if c.dicer_call is None:
        return S["amb.nosize"]
    near = nearest_class(c)
    m = _CRITERION.search(ev)
    if m:
        names = [n.strip() for n in m.group(1).split(",") if n.strip()]
        items = [S.get(f"amb.crit.{n}", n) for n in names]
        return S["amb.criterion"].format(cls=near, list="; ".join(items))
    return S["amb.size"].format(z=c.dicer_call, cls=near)


def _why_not(c, k, spec):
    """The first basic check of class k this region fails, where the TSV shows it."""
    S = STRINGS
    if spec is None:
        return S["whynot.default"].format(cls=k)
    if c.dicer_call is None:
        return S["whynot.nosize"].format(cls=k)
    if spec["size"] and not spec["size"][0] <= c.dicer_call <= spec["size"][1]:
        lo, hi = spec["size"]
        return S["whynot.size"].format(z=c.dicer_call, cls=k, lo=lo, hi=hi)
    span = c.read_span
    if span is None and c.start is not None and c.end is not None:
        span = c.end - c.start
    if spec["min_length"] is not None and span is not None and span < spec["min_length"]:
        return S["whynot.length"].format(span=span, min=spec["min_length"], cls=k)
    lim = spec["phase_max_p"]
    if lim is not None and (c.phase_pvalue is None or c.phase_pvalue > lim):
        return S["whynot.phasep"].format(p=_fmt(c.phase_pvalue), cls=k, lim=f"{lim:g}")
    if spec["by_phase"] and _level((c.e or {}).get("phase")) != "present":
        return S["whynot.phasing"].format(cls=k)
    return S["whynot.default"].format(cls=k)


def explain(c, kingdom="animal"):
    S = STRINGS
    specs = report_data.class_specs(kingdom)
    vals = _values(c)
    defining = specs.get(c.cls, {}).get("defining", set())
    mm = c.multimap_frac
    checks = []
    for k in _CHECKS[kingdom]:
        if k == "multimap":
            v, status, coloured = mm, _mm_level(mm), True
        else:
            v = vals[k][0]
            status, coloured = _level(*vals[k]), k in defining
        detail = S.get(f"check.{k}.detail.{kingdom}", S[f"check.{k}.detail"])
        checks.append({"key": k, "status": status, "title": S[f"check.{k}.title"],
                       "detail": detail.format(v=_fmt(v)) if v is not None else S["check.na"],
                       "coloured": coloured})
    note = "" if defining else S.get(f"check.note.{c.cls}", S["check.note.none"])
    unresolved = ambiguous_reason(c) if c.cls == "ambiguous" else ""
    if unresolved:
        verdict = S["verdict.ambiguous"].format(reason=unresolved)
    else:
        verdict = S.get(f"verdict.{c.cls}", S["verdict.generic"])
    if c.gates:
        why_not = [{"cls": k, "reason": _why_not(c, k, specs.get(k))}
                   for k, ok in c.gates.items() if not ok and k != c.cls]
        if not why_not:
            why_not = [{"cls": "", "reason": S["whynot.none"]}]
    else:
        why_not = [{"cls": "", "reason": S["whynot.unknown"]}]
    cautions = []
    if not c.has_evidence:
        cautions.append(S["caution.noevidence"])
    # ambiguous always carries review; its reason is given above instead
    if c.review and c.cls != "ambiguous":
        cautions.append(S["caution.review"])
    if mm is not None and mm >= DISPLAY_MULTIMAP_HIGH:
        cautions.append(S["caution.multimap"])
    return {"verdict": verdict, "checks": checks, "why_not": why_not, "note": note,
            "cautions": cautions, "unresolved": unresolved}


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
