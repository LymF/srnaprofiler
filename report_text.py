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
    "check.1u.title": "Starts with U",
    "check.10a.title": "A at position 10",
    "check.5a.title": "Starts with A",
    "check.strand.title": "Strand balance",
    "check.pingpong.title": "Ping-pong pairing",
    "check.phasing.title": "Regular phasing",
    "check.multimap.title": "Multimapping",
    "check.context": "not needed for this class",
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
    "warn.cap": "Only {n} clusters are shown in detail (each class keeps a share, then the highest confidence); the rest appear in the table only.",
    "warn.unique": "Only {pct}% of reads map to a single place; many small RNAs here may come from repeats, so classes are less certain.",
    "ui.about.title": "About",
    "ui.about.name": "sRNAProfiler",
    "ui.about.tagline": "Small RNA loci discovery and classification from sequencing data",
    "ui.about.version": "Version {v}",
    "ui.about.github": "https://github.com/LymF/srnaprofiler",
    "ui.about.gh.title": "GitHub",
    "ui.about.gh.sub": "Source code & documentation",
    "ui.about.install": "Install",
    "ui.about.install.1": "conda env create -f environment.yml",
    "ui.about.install.2": "conda activate srnaprofiler",
    "ui.fig.sizes": "Read size and 5' base",
    "ui.fig.sizes.tip": "Reads 18 to 35 nt, weighted by 1/NH, stacked by first base.",
    "ui.fig.phasing": "Phasing",
    "ui.fig.phasing.tip": "Pairs of 5' ends on the same strand, by distance.",
    "ui.fig.cov": "Coverage along the cluster",
    "ui.fig.cov.help": "Mean reads per base in each bin, not weighted by multimapping.",
    "ui.fig.pp": "Ping-pong overlap",
    "ui.fig.pp.plant": "Overlap between strands",
    "ui.fig.pp.tip": "Pairs of reads on opposite strands, by 5'-5' overlap.",
    "ui.fig.ends": "5' ends along the cluster",
    "ui.fig.ends.tip": "Read 5' ends in each bin, weighted by 1/NH.",
    "ui.fig.logo": "Base logo, {z}\u00a0nt reads",
    "ui.fig.logo.tip": "Share of each base at each position, reads of the dominant size.",
    "ui.fig.seqs": "Most abundant sequences",
    "ui.fig.seqs.tip": "Distinct sequences, reads weighted by 1/NH.",
    "ui.fig.ctx": "Library context",
    "ui.fig.ctx.tip": "Among the clusters in this report.",
    "ui.fig.n": "n\u00a0=\u00a0{n} reads",
    "ui.fig.nseq": "top {n}",
    "ui.fig.nclusters": "n\u00a0=\u00a0{n} clusters",
    "ui.fig.table": "Table",
    "ui.fig.chart": "Figure",
    "ui.fig.svg": "Download SVG",
    "ui.fig.precision": "Values reproduce the drawn ones within 1% of the panel maximum.",
    "ui.fig.exact": "Exact counts.",
    "ui.fig.mirror": "Sense (+) up, antisense (\u2212) down",
    "ui.fig.period": "Period {p}\u00a0nt",
    "ui.fig.bin": "Bins of {b}\u00a0nt",
    "ui.fig.pp10": "10\u00a0nt overlap",
    "ui.fig.ppother": "Other overlaps",
    "ui.fig.others": "Other clusters",
    "ui.fig.nostrand": "No reads on this strand",
    "ui.fig.nodata": "Nothing to draw on this strand",
    "ui.fig.nopairs": "No pairs between the strands",
    "ui.fig.noseqs": "No sequences were recorded for this cluster.",
    "ui.strand.sense": "Sense (+)",
    "ui.strand.antisense": "Antisense (\u2212)",
    "ui.axis.reads": "Reads",
    "ui.axis.length": "Small RNA length (nt)",
    "ui.axis.distance": "5'-5' distance (nt)",
    "ui.axis.overlap": "5'-5' overlap (nt)",
    "ui.axis.pairs": "Weighted pairs",
    "ui.axis.position": "Position on {chrom}",
    "ui.axis.depth": "Reads per base",
    "ui.axis.ends": "5' ends",
    "ui.axis.logo": "Position in the read (nt)",
    "ui.axis.share": "Share of reads",
    "ui.axis.logreads": "Reads per cluster (log scale)",
    "ui.axis.conf": "Confidence within the class",
    "ui.axis.clusters": "Clusters",
    "ui.unit.nt": "{v}\u00a0nt",
    "ui.unit.pct": "{v}%",
    "ui.unit.plus": "{v}% on +",
    "ui.unit.z": "z\u00a0=\u00a0{v}",
    "ui.unit.p": "p\u00a0=\u00a0{v}",
    "ui.unit.k": "{v}K",
    "ui.unit.m": "{v}M",
    "ui.tip.dist": "{x}\u00a0nt apart",
    "ui.tip.overlap": "{x}\u00a0nt overlap",
    "ui.tip.pos": "{start}\u2013{end}",
    "ui.tip.logo": "Position {x}",
    "ui.tip.range": "{lo}\u2013{hi}",
    "ui.tip.this": "This cluster",
    "ui.tip.ctx": "Top {pct}% of clusters by reads",
    "ui.tip.conf": "Higher than {pct}% of {cls} clusters",
    "ui.ctx.reads": "Reads rank",
    "ui.ctx.top": "top {pct}%",
    "ui.ctx.sub": "{lb}: {v}",
    "ui.ctx.conf": "Confidence rank",
    "ui.ctx.pct": "percentile {pct}",
    "ui.ctx.conf.na": "No rank within the class",
    "ui.an.legend": ("Annotation is evidence: it can raise or lower the confidence and flag a "
                     "locus for review, but it never changes the class. Annotation is "
                     "incomplete, so no annotation means nothing."),
    "ui.an.title": "Annotation evidence per class",
    "ui.an.note": ("Loci listed in this report, by whether the annotation supports the class "
                   "call, conflicts with it, or says nothing."),
    "ui.an.th.class": "Class",
    "ui.an.th.supports": "Supports",
    "ui.an.th.conflicts": "Conflicts",
    "ui.an.th.neutral": "Neutral",
    "ui.an.filter": "Annotation",
    "ui.an.all": "All",
    "ui.gm.title": "Genomic context",
    "ui.gm.track.cv": "Coverage",
    "ui.gm.track.ends": "5' ends",
    "ui.gm.track.g": "Genes",
    "ui.gm.track.t": "Transposons",
    "ui.gm.track.r": "Structural RNA",
    "ui.gm.track.p": "PHAS loci",
    "ui.gm.track.orf": "Predicted ORFs",
    "ui.gm.orf.note": "experimental",
    "ui.gm.orf.empty": "No ORFs of at least {n} aa outside annotated CDS",
    "ui.gm.orf.len": "{n} aa",
    "ui.gm.orf.tip": "Predicted, not annotated",
    "ui.gm.orf.legend": "Predicted ORF (experimental)",
    "ui.gm.th.orf": "Predicted ORFs",
    "ui.gm.empty": "Nothing annotated in this window",
    "ui.gm.nolayer": "No {kind} layer given",
    "ui.gm.more": "{n} more not shown",
    "ui.gm.tip.type": "Type",
    "ui.gm.tip.strand": "Strand",
    "ui.gm.tip.pos": "Position",
    "ui.gm.tip.overlap": "Overlaps the locus",
    "ui.gm.tip.nostrand": "Strand unknown",
    "ui.gm.zoom.in": "Zoom in",
    "ui.gm.zoom.out": "Zoom out",
    "ui.gm.zoom.reset": "Reset",
    "ui.gm.pan.left": "Pan left",
    "ui.gm.pan.right": "Pan right",
    "ui.gm.locus": "Zoom to locus",
    "ui.gm.band": "Locus",
    "ui.gm.ends.unit": "Starts/base",
    "ui.gm.legend.overlap": "Overlaps the locus",
    "ui.gm.legend.band": "Locus",
    "ui.gm.legend.cont": "Continues past the view",
    "ui.gm.legend.marks": "Marks",
    "ui.gm.range": "{chrom}:{start}-{end}, {n}\u00a0bp",
    "ui.gm.legend.te": "Transposon class",
    "ui.gm.legend.size": "Read size class",
    "ui.gm.kind.g": "Gene",
    "ui.gm.te.other": "Other",
    "ui.gm.th.track": "Track",
    "ui.gm.th.name": "Name",
    "ui.gm.th.start": "Start",
    "ui.gm.th.end": "End",
    "ui.review": "Flagged for review",
    "ui.review.tip": "The classifier marked this locus for a manual look.",
    "ui.brand": "srnaprofiler",
    "ui.title.default": "srnaprofiler report",
    "ui.na": "n/a",
    "ui.unclassified": "unclassified",
    "ui.tab.overview": "Overview",
    "ui.tab.clusters": "Clusters",
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
    "ui.help.classes": "each bar is the share of clusters in one class. Ambiguous is not an error: the region did not meet the requirements of any one class.",
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
    "ui.sort.evidence": "Sort: evidence",
    "ui.sort.reads": "Sort: reads",
    "ui.sort.name": "Sort: name",
    "ui.floor.label": "Minimum reads ",
    "ui.sort.aria": "Sort clusters by",
    "ui.floor.aria": "Minimum reads",
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
    "ui.th.length": "Length (nt)",
    "ui.th.strand": "Strand",
    "ui.th.total": "Total",
    "ui.th.seq": "Sequence",
    "ui.th.base": "5' base",
    "ui.th.bin": "Range",
    "ui.conf.amb": "n/a (nearest class: {cls})",
    "ui.conf.amb.none": "n/a",
    "ui.where": "{chrom}:{start}-{end}",
    "ui.page.prev": "Previous",
    "ui.page.next": "Next",
    "ui.page.info": "Page {p} of {pages} ({n} clusters)",
    "ui.dots.legend": "Evidence dots, in order: {list}.",
    "ui.fig.na": "Figures are not available for this cluster.",
    "ui.fig.noreads": "No reads were found in this region of the BAM, so there is nothing to draw.",
    "ui.fig.err": "This region could not be read from the BAM (see the warning in the log), so there is nothing to draw.",
    "ui.detail.back": "Back to the list",
    "ui.detail.sub": "{where} \u00b7 {reads} reads",
    "ui.detail.size": " \u00b7 dominant size {z}\u00a0nt",
    "ui.detail.conflabel": "confidence {v}",
    "ui.detail.conf": "Evidence for this class against the others ({k} classes compared); {min} means no preference. Not a probability.",
    "ui.detail.conf.plant": "Evidence for this class against the others. Up to {k} classes compete; a class that fails its phasing test drops out, so the no-preference value is 1 divided by the classes left (as low as {min} when all {k} compete). Not a probability.",
    "ui.detail.conf.ml": "The trained model's score for its top class. Not a probability.",
    "ui.detail.conf.mirna": "Set to 1 for miRNAs: the class comes from the miRDeep2 hairpin search.",
    "ui.detail.conf.amb": "Not shown: the stored value belongs to the nearest class, which was not assigned.",
    "ui.footer.animal": "The coloured evidence dots use display-only cutoffs, not the classifier's thresholds. 1U, strand, ping-pong and phasing scores: present at {present} or more, weak at {weak} or more. A at position 10: present at {p10} or more, weak at {w10} or more. Size consistency: present at {psize} or more, weak at {wsize} or more. Multimapping: moderate at {mmod} or more, high at {mhigh} or more. They only decide how a bar is worded in this report.",
    "ui.footer.plant": "The coloured evidence dots use display-only cutoffs, not the classifier's thresholds. 1U, 5'A, strand and phasing scores: present at {present} or more, weak at {weak} or more. Size consistency: present at {psize} or more, weak at {wsize} or more. Multimapping: moderate at {mmod} or more, high at {mhigh} or more. They only decide how a bar is worded in this report.",
    "check.1u_pingpong.title": "Starts with U or ping-pong",
    "check.length.title": "Cluster length",
    "check.density.title": "Read density",
    "ck.title": "Why {cls}",
    "ck.title.amb": "Why this is unresolved: checks for the nearest class, {cls}",
    "ck.mirna": "The class comes from the miRDeep2 hairpin search, so no read-signature checks apply.",
    "ck.ml": "The class came from the trained model. These are the class's basic checks; the evidence rows are for reference.",
    "ck.note": "Evidence checks have no cutoff of their own: the classifier weighs them in a score, and a row is met here when its effect is above chance. The size window and cluster criteria are hard limits.",
    "ck.th.check": "Check",
    "ck.th.observed": "Observed",
    "ck.th.rule": "Rule",
    "ck.th.status": "Result",
    "ck.st.1": "Met",
    "ck.st.0": "Not met",
    "ck.st.2": "Not available",
    "ck.st.3": "Not tested",
    "ck.thr.size": "Dominant size of {lo} to {hi} nt",
    "ck.thr.size1": "Dominant size of {lo} nt",
    "ck.thr.1u": "Share of 5' U above chance (95% confidence)",
    "ck.thr.5a": "Share of 5' A above chance (95% confidence)",
    "ck.thr.strand": "Minority strand above 10% of the reads (95% confidence)",
    "ck.thr.pingpong": "10 nt overlap above the other overlaps (95% confidence)",
    "ck.thr.1u_pingpong": "Either 5' U above chance or the 10 nt overlap above the other overlaps (95% confidence)",
    "ck.thr.phasing": "Cycles in register above chance (95% confidence)",
    "ck.thr.phasing_p": "Phasing test p of {p} or less, and cycles in register above chance (95% confidence)",
    "ck.thr.length": "Reads span at least {n} bp",
    "ck.thr.density": "Reads denser than chance (p of {p} or less)",
    "ck.obs.size": "{z} nt",
    "ck.obs.nosize": "no size call",
    "ck.obs.1u": "{pct}% start with U",
    "ck.obs.5a": "{pct}% start with A",
    "ck.obs.strand": "{pct}% on the minority strand",
    "ck.obs.pingpong": "z = {z}",
    "ck.obs.phasing": "p = {p}, {n} cycles",
    "ck.obs.strength": ", strength {e}",
    "ck.obs.span": "{n} bp",
    "ck.obs.density": "p = {p}",
    "ui.ov.lq": "Library quality",
    "ui.ov.lq.note": "Read lengths of the mapped reads in the BAM, each read counted once (weighted by 1/NH). Shares are of mapped reads.",
    "ui.fig.lq": "Read length, whole library",
    "ui.fig.lq.tip": "Mapped reads of 16 to 40 nt, weighted by 1/NH, with one bin each for shorter and longer reads.",
    "lq.band": "20 to 25 nt",
    "lq.band.label": "20 to 25 nt: {pct}",
    "lq.share.label": "Reads of 20 to 25 nt",
    "lq.share.1": "OK: at least 25% of the reads",
    "lq.share.0": "Check: under 25% of the reads",
    "lq.dominant.label": "Dominant read length",
    "lq.dominant.1": "OK: within 18 to 35 nt, the sizes the classifier keeps",
    "lq.dominant.0": "Check: outside 18 to 35 nt, the sizes the classifier keeps",
    "lq.short.label": "Reads shorter than 18 nt",
    "lq.short.1": "OK: 25% of the reads or fewer",
    "lq.short.0": "Check: more than 25% of the reads",
    "lq.mapped": "Mapped reads counted",
    "lq.mapped.sub": "weighted by 1/NH",
    "lq.unmapped": "Unmapped records",
    "lq.unmapped.sub": "in the BAM, not in the figure",
    "lq.tr.before": "Reads before trimming",
    "lq.tr.after": "Reads after trimming",
    "lq.tr.after.sub": "{pct} kept",
    "lq.tr.adapter": "Reads with an adapter trimmed",
    "lq.tr.adapter.sub": "{pct} of the reads before trimming",
    "mi.title": "miRNA hairpin",
    "mi.tip": "The precursor with its mature, loop and star segments, the predicted fold and the reads inside it.",
    "mi.strand": "Shown 5' to 3' on the {st} strand of {chrom}.",
    "mi.seg.mature": "mature",
    "mi.seg.loop": "loop",
    "mi.seg.star": "star",
    "mi.kind.other": "other",
    "mi.lg.mismatch": "Mismatch to the precursor (upper case)",
    "mi.axis.cov": "Reads per base",
    "mi.struct": "Fold by RNAfold, minimum free energy {mfe} kcal/mol.",
    "mi.nostruct.m": "Structure not available: RNAfold not installed.",
    "mi.nostruct.f": "Structure not available: RNAfold could not fold this sequence.",
    "mi.noseq": "The miRNA table has no precursor sequence for this locus (miRScore did not score it), so the hairpin cannot be drawn. The numbers are below.",
    "mi.nowin": "The precursor could not be placed on the genome from its coordinates, so its reads are not shown.",
    "mi.nostrand": "The strand of this precursor is not known (the miRDeep2 table gave '.'), so its reads are not shown.",
    "mi.noreads": "No reads of the precursor strand fall fully inside the precursor in the BAM.",
    "mi.more": "{n} more sequences ({reads} reads) are not drawn; the table view lists the drawn ones.",
    "mi.anti": "Not shown: {anti} reads on the opposite strand and {ovh} reads that run past the precursor ends.",
    "mi.reads.note": "Reads of the precursor strand that lie fully inside the precursor, collapsed by position and sequence and weighted by 1/NH. Lower case matches the precursor, upper case is a mismatch. A read is mature or star when at least half of it overlaps that segment. The sequence includes the flanks given to miRScore, and the depth counts every collapsed read.",
    "mi.th.pos": "Position",
    "mi.th.base": "Base",
    "mi.th.fold": "Fold",
    "mi.th.seg": "Segment",
    "mi.th.pair": "Paired with",
    "mi.d2.mfe": "MFE {mfe} kcal/mol",
    "mi.d2.base": "Base {pos}: {base}",
    "mi.d2.pair": "paired with {pos}",
    "mi.d2.note": "Secondary structure drawn from the RNAfold dot-bracket with the ViennaRNA layout. Every 10th base is numbered; letters are drawn only when there is room, and each base shows its position and letter on hover.",
    "mi.th.depth": "Depth",
    "mi.th.start": "Start",
    "mi.th.end": "End",
    "mi.th.seq": "Sequence",
    "mi.th.reads": "Reads",
    "mi.th.type": "Type",
    "mi.num.title": "miRNA numbers",
    "mi.num.note": "From the miRNA table of the same run.",
    "mi.num.result": "miRScore result",
    "mi.num.mature_reads": "Mature reads",
    "mi.num.star_reads": "Star reads",
    "mi.num.precision": "Precision",
    "mi.num.mismatches": "Mismatches in the duplex",
    "mi.num.bulges": "Bulges",
    "mi.num.precursor_len": "Precursor length (with flanks)",
    "mi.num.locus_len": "Hairpin locus length",
    "mi.num.mirdeep_score": "miRDeep2 score",
    "mi.num.mirdeep_prob": "miRDeep2 probability",
    "mi.num.randfold": "Randfold p-value",
    "mi.num.rfam_alert": "Rfam alert",
    "mi.num.mirbase_hit": "miRBase hit",
    "mi.num.mirbase_same_seed": "miRBase hit with the same seed",
    "mi.crit.title": "Plant criteria, shown for reference",
    "mi.crit.note": "Three criteria of Axtell and Meyers 2018 that the table can answer, shown for every kingdom. They are not used to call the class.",
    "mi.crit.len": "Precursor length",
    "mi.crit.mm": "Mismatches in the duplex",
    "mi.crit.prec": "Precision",
    "mi.rule.len": "Hairpin locus of {n} nt or less",
    "mi.rule.mm": "At most {n} mismatches",
    "mi.rule.prec": "At least {p}% of the reads at the miRNA and its star",
    "mi.obs.nt": "{n} nt",
    "mi.obs.pct": "{p}%",
}

# Report check keys shown per kingdom; the first five are the list's evidence dots.
DOTS = {"animal": ("size", "1u", "10a", "pingpong", "phasing"),
        "plant": ("size", "1u", "5a", "strand", "phasing")}
_CHECKS = {"animal": ("size", "1u", "10a", "strand", "pingpong", "phasing", "multimap"),
           "plant": ("size", "1u", "5a", "strand", "phasing", "multimap")}
# kingdom-specific text: role -> STRINGS key
_KTEXT = {"animal": {"sizes": "ui.help.sizes.animal", "pp_title": "ui.fig.pp"},
          "plant": {"sizes": "ui.help.sizes.plant", "pp_title": "ui.fig.pp.plant"}}
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


def trail_nearest(c):
    """Nearest class named in the evidence trail, with or without a size call."""
    m = _NEAREST.search(c.evidence or "")
    return m.group(1) if m else None


def nearest_class(c):
    """Nearest class named in the evidence trail; None without a size call or a trail."""
    return None if c.dicer_call is None else trail_nearest(c)


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


def verdict_key(c):
    """STRINGS key of the one-line verdict; the front end formats it."""
    key = f"verdict.{c.cls}"
    return key if key in STRINGS else "verdict.generic"


def strip(c, kingdom="animal"):
    """Measured values for the strip: [key, level, value], value None when absent.

    The level is the display-only status (see DISPLAY_*), the value the raw number.
    """
    vals = _values(c)
    raw = {"size": c.dicer_call, "1u": c.pct_1u, "10a": c.pct_10a, "5a": c.pct_5a,
           "strand": c.frac_top, "pingpong": c.pingpong_z, "phasing": c.phase_pvalue}
    out = []
    for k in _CHECKS[kingdom]:
        if k == "multimap":
            lv, v = _mm_level(c.multimap_frac), c.multimap_frac
        else:
            lv, v = _level(*vals[k]), raw[k]
        out.append([k, lv, v])
    return out


def explain(c, kingdom="animal"):
    """Verdict key, strip and (ambiguous only) the reason the class was withheld."""
    return {"verdict": verdict_key(c), "strip": strip(c, kingdom),
            "unresolved": ambiguous_reason(c) if c.cls == "ambiguous" else ""}


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
