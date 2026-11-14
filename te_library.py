"""Build a transposon library for a genome: RepeatModeler2, then drop contaminants.

`annotate-loci --te-library` needs a library. RepeatModeler2 families include
multicopy gene families and structural RNAs; this module removes them by
aligning each family (minimap2 -x asm20) to the exons of the non-TE genes and
to any user FASTA. Optional and separate: classify/report never call it.
"""
import csv
import gzip
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

from utils import log, which

MIN_IDENTITY = 0.80
MIN_COVERAGE = 0.50
RNA_CLASSES = ("rRNA", "tRNA", "snRNA", "scRNA", "srpRNA")


def _die(msg):
    log(msg, "ERROR")
    sys.exit(1)


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path):
    """Stream (header, sequence) pairs; header is the line without '>'."""
    head, parts = None, []
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if head is not None:
                    yield head, "".join(parts)
                head, parts = line[1:], []
            elif head is not None:
                parts.append(line.strip())
    if head is not None:
        yield head, "".join(parts)


def family_name(header):
    return header.split(None, 1)[0] if header.split() else header


def family_class(name):
    """RepeatModeler naming `name#LTR/Gypsy` -> 'LTR'; no '#' -> ''."""
    if "#" not in name:
        return ""
    return name.split("#", 1)[1].split("/", 1)[0]


def write_nuclear_fasta(genome, out_fa, exclude):
    """Copy the genome without the excluded contigs, streaming."""
    exclude = set(exclude)
    seen, keep = set(), True
    with _open_text(genome) as fh, open(out_fa, "w") as out:
        for line in fh:
            if line.startswith(">"):
                name = family_name(line[1:].rstrip("\r\n"))
                keep = name not in exclude
                if name in exclude:
                    seen.add(name)
            if keep:
                out.write(line.rstrip("\r\n") + "\n")
    return sorted(exclude - seen)


def _contig_names(genome):
    with _open_text(genome) as fh:
        return {family_name(l[1:].rstrip("\r\n")) for l in fh if l.startswith(">")}


def write_gene_exons(genes, genome, out_fa):
    """FASTA of the exons of the non-TE genes of an annotation. Returns the count.

    The genome is streamed contig by contig, so no .fai is needed.
    """
    from annotation import _drop_te_trees, _gene_class, parse_annotation
    contigs = _contig_names(genome)
    unknown = {}
    feats = parse_annotation(genes, known_contigs=contigs, unknown=unknown)
    n = 0
    with open(out_fa, "w") as out:
        for head, seq in read_fasta(genome):
            contig = family_name(head)
            for f in _drop_te_trees(feats.get(contig, [])):
                if _gene_class(f.ftype) != "exon" or f.end <= f.start:
                    continue
                n += 1
                out.write(f">exon{n}|{contig}:{f.start}-{f.end}\n{seq[f.start:f.end]}\n")
    if n == 0:
        skipped = sum(unknown.values())
        eg = ", ".join(list(unknown)[:3])
        _die(f"--genes {genes}: no exon of a non-TE gene was extracted from the genome. "
             + (f"{skipped} features lie on {len(unknown)} contigs absent from the "
                f"genome (e.g. {eg}); check contig naming, e.g. Chr1 vs 1"
                if unknown else "the file has no gene exon rows"))
    return n


def concat_fastas(paths, out_fa):
    with open(out_fa, "w") as out:
        for p in paths:
            with _open_text(p) as fh:
                for line in fh:
                    out.write(line.rstrip("\r\n") + "\n")


def _union_len(spans):
    total, end = 0, None
    start = None
    for s, e in sorted(spans):
        if end is None or s > end:
            if end is not None:
                total += end - start
            start, end = s, e
        else:
            end = max(end, e)
    if end is not None:
        total += end - start
    return total


def parse_paf(paf_path, lengths):
    """{family: (coverage, matched target)} using hits with identity >= MIN_IDENTITY.

    identity = matches / alignment block length; coverage = union of the query
    spans of the family's passing hits / family length.
    """
    spans, best = {}, {}
    with open(paf_path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 11:
                continue
            try:
                qs, qe, nm, blk = int(c[2]), int(c[3]), int(c[9]), int(c[10])
            except ValueError:
                continue
            if blk <= 0 or nm / blk < MIN_IDENTITY:
                continue
            q = c[0]
            spans.setdefault(q, []).append((qs, qe))
            if qe - qs > best.get(q, (0, ""))[0]:
                best[q] = (qe - qs, c[5])
    out = {}
    for q, sp in spans.items():
        ln = lengths.get(q, 0)
        if ln:
            out[q] = (_union_len(sp) / ln, best[q][1])
    return out


def _run(cmd, cwd=None, stdout_path=None):
    """The one subprocess seam; output goes to a file, never into memory."""
    with open(stdout_path, "w") if stdout_path else open(os.devnull, "w") as fh:
        subprocess.run(cmd, check=True, cwd=cwd, stdout=fh,
                       stderr=subprocess.STDOUT if stdout_path else None)


def _run_dirs(outdir):
    """RepeatModeler run directories, newest first.

    A run leaves many RM_* directories (RepeatClassifier makes temporary ones too); only a real
    run directory has a `rmod.log` or a `round-1` folder."""
    dirs = [d for d in outdir.glob("RM_*")
            if (d / "rmod.log").exists() or (d / "round-1").exists()]
    return sorted(dirs, key=lambda d: d.stat().st_mtime, reverse=True)


def _run_dirs(outdir):
    """RepeatModeler run directories, newest first.

    A run leaves many RM_* directories (RepeatClassifier makes temporary ones too); only a real
    run directory has a `rmod.log` or a `round-1` folder."""
    dirs = [d for d in outdir.glob("RM_*")
            if (d / "rmod.log").exists() or (d / "round-1").exists()]
    return sorted(dirs, key=lambda d: d.stat().st_mtime, reverse=True)


def _repeatmodeler_library(outdir, prefix):
    fam = outdir / f"{prefix}-families.fa"
    if fam.exists():
        return fam
    for d in _run_dirs(outdir):
        c = d / "consensi.fa.classified"
        if c.exists():
            return c
    return None


def _pick_run_dir(outdir, prefix):
    """The RM_* dir to resume: the one whose rmod.log mentions the prefix, else the newest."""
    dirs = _run_dirs(outdir)
    for d in dirs:
        lg = d / "rmod.log"
        try:
            if lg.exists() and prefix in lg.read_text(errors="replace"):
                return d
        except OSError:
            pass
    return dirs[0] if dirs else None


def _validate(args, filtering, genes, excl, families_in):
    for label, path in [("--genes", genes), ("--families", families_in),
                        ("--genome", args.genome)] + [("--exclude-fasta", p) for p in excl]:
        if path and not os.path.isfile(path):
            _die(f"{label} {path}: no such file")
    for p in excl:
        if not any(seq for _, seq in read_fasta(p)):
            _die(f"--exclude-fasta {p} holds no sequence")


def run_te_library(args):
    outdir = Path(args.outdir)
    prefix = args.prefix
    final = outdir / f"{prefix}-families.fa"
    raw = outdir / f"{prefix}-families.raw.fa"
    report = outdir / f"{prefix}.te-library.report.tsv"
    genes = getattr(args, "genes", None)
    excl = list(getattr(args, "exclude_fasta", None) or [])
    families_in = getattr(args, "families", None)
    dfam = getattr(args, "dfam_species", None)
    filtering = bool(genes or excl)

    if report.exists():
        _die(f"{report} exists; refusing to overwrite a finished library. "
             f"Remove it or pick another --prefix/--outdir")
    if families_in and final.exists() and not os.path.samefile(families_in, final):
        _die(f"{final} exists; refusing to overwrite it. Remove it or pick "
             f"another --prefix/--outdir")
    _validate(args, filtering, genes, excl, families_in)
    reuse = (not families_in) and (raw.exists() or final.exists())
    need = [] if (families_in or reuse) else [
        ("BuildDatabase", "RepeatModeler (`conda install -c bioconda repeatmodeler`)"),
        ("RepeatModeler", "RepeatModeler (`conda install -c bioconda repeatmodeler`)")]
    if filtering:
        need.append(("minimap2", "`conda install -c bioconda minimap2`"))
    for tool, hint in need:
        if not which(tool):
            _die(f"{tool} is not on PATH. Install it with {hint}; it is "
                 f"optional for the rest of srnaprofiler")
    organelle = [x.strip() for x in (getattr(args, "organelle_contigs", None) or "").split(",")
                 if x.strip()]
    outdir.mkdir(parents=True, exist_ok=True)

    # 1-2. library to filter
    if families_in:
        src = Path(families_in)
        log(f"using the library {src}; RepeatModeler is skipped")
    elif reuse:
        if not raw.exists():
            os.replace(final, raw)
        src = raw
        log(f"{src.name} already exists and no report: RepeatModeler finished earlier "
            f"but filtering did not; RepeatModeler is skipped and this library is filtered")
    else:
        genome = args.genome
        if organelle:
            nuc = outdir / f"{prefix}.nuclear.fa"
            missing = write_nuclear_fasta(args.genome, str(nuc), organelle)
            if missing:
                nuc.unlink()
                _die(f"--organelle-contigs: {', '.join(missing)} not in {args.genome}")
            genome = str(nuc)
            log(f"--organelle-contigs: dropped {', '.join(organelle)} -> {nuc}")
        db_log = outdir / f"{prefix}.builddatabase.log"
        rm_log = outdir / f"{prefix}.repeatmodeler.log"
        log(f"BuildDatabase on {genome} (log {db_log})")
        try:
            _run(["BuildDatabase", "-name", prefix, os.path.abspath(genome)],
                 cwd=str(outdir), stdout_path=str(db_log))
        except (subprocess.CalledProcessError, OSError) as exc:
            _die(f"BuildDatabase failed: {exc} (see {db_log})")
        cmd = ["RepeatModeler", "-database", prefix, "-threads", str(args.threads),
               "-LTRStruct"]
        if getattr(args, "resume", False):
            rd = _pick_run_dir(outdir, prefix)
            if rd is not None:
                cmd += ["-recoverDir", str(rd)]
                log(f"resuming RepeatModeler from {rd}")
            else:
                log("--resume: no RM_* run directory to resume; starting fresh")
        log(f"RepeatModeler running (hours); log: {rm_log}. If interrupted, rerun "
            f"the same command with --resume")
        try:
            _run(cmd, cwd=str(outdir), stdout_path=str(rm_log))
        except (subprocess.CalledProcessError, OSError) as exc:
            _die(f"RepeatModeler failed: {exc} (see {rm_log}; rerun with --resume)")
        src = _repeatmodeler_library(outdir, prefix)
        if src is None:
            _die(f"RepeatModeler produced no {final.name} or RM_*/consensi.fa.classified "
                 f"in {outdir}")
    if src.exists() and final.exists() and os.path.samefile(src, final):
        os.replace(final, raw)
        src = raw
        log(f"the input library is kept as {raw}")

    # families metadata (streamed; sequences are not kept)
    lengths, classes, order, verdict = {}, {}, [], {}
    dup = 0
    for head, seq in read_fasta(src):
        nm = family_name(head)
        if nm in lengths:
            dup += 1
            continue
        lengths[nm] = len(seq)
        classes[nm] = family_class(nm)
        order.append(nm)
        if not seq:
            verdict[nm] = ("empty", "")
    if not order:
        _die(f"the input library {src} has 0 families")
    if dup:
        log(f"{dup} duplicate family names in {src}; the first of each is kept", "WARN")
    if any(r == "empty" for r, _ in verdict.values()):
        log(f"{sum(1 for r, _ in verdict.values() if r == 'empty')} zero-length "
            f"families dropped", "WARN")
    log(f"{len(order)} input families from {src}")

    # 3. filter
    for nm in order:
        if nm not in verdict and classes[nm] in RNA_CLASSES:
            verdict[nm] = (classes[nm], "")
    note = ""
    if filtering:
        excl_fa = outdir / f"{prefix}.exclude.fa"
        parts = []
        if genes:
            ex = outdir / f"{prefix}.gene-exons.fa"
            n = write_gene_exons(genes, args.genome, str(ex))
            log(f"{n} exons of non-TE genes extracted from {genes}")
            parts.append(ex)
        parts += [Path(p) for p in excl]
        concat_fastas(parts, str(excl_fa))
        paf = outdir / f"{prefix}.te-library.paf"
        try:
            _run(["minimap2", "-x", "asm20", "-c", "-s", "60", "-N", "5",
                  "-t", str(args.threads), str(excl_fa), str(src)],
                 stdout_path=str(paf))
        except (subprocess.CalledProcessError, OSError) as exc:
            _die(f"minimap2 failed: {exc}")
        for nm, (cov, tname) in parse_paf(paf, lengths).items():
            if cov >= MIN_COVERAGE and nm not in verdict:
                verdict[nm] = ("gene", tname)
        for p in (excl_fa, paf):
            p.unlink(missing_ok=True)
    else:
        note = ("no --genes or --exclude-fasta: no sequence filtering was done "
                "(only RNA-class families were dropped)")
        log(note, "WARN")

    # write library: kept families streamed from the source
    tmp = outdir / f"{prefix}-families.fa.tmp"
    kept = Counter()
    names = set()
    n_dfam = 0
    with open(tmp, "w") as out:
        for head, seq in read_fasta(src):
            nm = family_name(head)
            if nm in names:
                continue
            names.add(nm)
            if nm in verdict:
                continue
            out.write(f">{head}\n{seq}\n")
            kept[classes[nm] or "unclassified"] += 1
        if dfam:
            n_dfam = _append_dfam(dfam, outdir, prefix, out, names)
            if n_dfam:
                kept["dfam"] += n_dfam
    os.replace(tmp, final)

    with open(report, "w", newline="") as fh:
        fh.write(f"# source library: {src}\n")
        if note:
            fh.write(f"# {note}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["family", "length", "status", "reason", "matched_sequence"])
        for nm in order:
            reason, matched = verdict.get(nm, ("", ""))
            w.writerow([nm, lengths[nm], "dropped" if nm in verdict else "kept",
                        reason, matched])
    dropped = Counter(r for r, _ in verdict.values())
    log(f"kept {sum(kept.values())} families by class: "
        f"{dict(kept) or '-'}; dropped by reason: {dict(dropped) or '-'}")
    log(f"library -> {final}; report -> {report}; use it as "
        f"annotate-loci --te-library {final}")


def _append_dfam(species, outdir, prefix, out, names):
    if not which("famdb.py"):
        log("famdb.py is not on PATH; the Dfam merge is skipped", "WARN")
        return 0
    dfa = outdir / f"{prefix}.dfam.fa"
    try:
        _run(["famdb.py", "families", "--curated", "-a", "-ad", "-f", "fasta_name",
              species], stdout_path=str(dfa))
    except (subprocess.CalledProcessError, OSError) as exc:
        log(f"famdb.py failed ({exc}); the Dfam merge is skipped", "WARN")
        dfa.unlink(missing_ok=True)
        return 0
    n = rna = 0
    for head, seq in read_fasta(dfa):
        if not seq:
            continue
        if family_class(family_name(head)) in RNA_CLASSES:
            rna += 1
            continue
        nm = "dfam|" + family_name(head)
        if nm in names:
            continue
        names.add(nm)
        out.write(f">dfam|{head}\n{seq}\n")
        n += 1
    if rna:
        log(f"{rna} Dfam families of RNA classes ({', '.join(RNA_CLASSES)}) dropped")
    if n == 0:
        log(f"famdb.py returned no families for {species!r}; nothing merged", "WARN")
    else:
        log(f"{n} curated Dfam families for {species} appended")
    return n
