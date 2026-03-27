#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
genome-repeat-summary.py
Quantifica a presença de regiões repetitivas no genoma usando RepeatMasker.

Entradas:
  --chrom-sizes   : TSV com duas colunas: chrom \t length
  --rm-gff        : 1+ arquivos GFF (GFF2/GFF3) do RepeatMasker (opcional)
  --rm-out        : 1+ arquivos .out do RepeatMasker (opcional)
  (Passe pelo menos um entre --rm-gff/--rm-out)

Saídas (prefixadas por --out-prefix):
  <prefix>.summary.tsv             -> resumo global (comprimento do genoma, TE bp, fração TE)
  <prefix>.per_chrom.tsv           -> por cromossomo (len, TE bp, fração TE)
  <prefix>.te_by_class.tsv         -> por classe (bp, % do genoma, % da carga TE)
  <prefix>.te_by_family.tsv        -> por família (bp, % do genoma, % da carga TE)
  <prefix>.top_families.tsv        -> top N famílias (métrica escolhida)
  <prefix>.top_families.pdf        -> barplot (PDF) das top famílias
  <prefix>.class_donut.pdf         -> donut (PDF) por classe

Requisitos:
  - bedtools no PATH, pybedtools, pandas
  - (opcional) Rscript com ggplot2, dplyr, scales, viridis (para gerar os PDFs)
"""

import argparse, sys, os, re, glob, shutil, subprocess
import pandas as pd
from collections import defaultdict
try:
    import pybedtools
except Exception:
    sys.stderr.write(
        "[ERRO] pybedtools não importou.\n"
        "Instale com: conda install -c conda-forge -c bioconda bedtools pybedtools pandas\n"
    )
    raise

# --------------------------- utils ---------------------------

def read_chrom_sizes(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"): continue
            parts = line.strip().split()
            if len(parts) < 2: continue
            chrom, length = parts[0], int(float(parts[1]))
            rows.append((chrom, length))
    if not rows:
        raise SystemExit("[ERRO] chrom-sizes vazio.")
    df = pd.DataFrame(rows, columns=["chrom","length"])
    return df

def parse_class_family(cf):
    if not cf or cf == ".": return "NA","NA"
    parts = cf.split("/", 1)
    if len(parts)==2: return (parts[0].strip() or "NA"), (parts[1].strip() or "NA")
    cf = cf.strip()
    return cf, cf

RE_CLASS  = re.compile(r'(?:^|[;\s])class\s*=?\s*"?([^";]+)"?', re.I)
RE_FAMILY = re.compile(r'(?:^|[;\s])family\s*=?\s*"?([^";]+)"?', re.I)
RE_TARGET = re.compile(r'(?:^|[;\s])(?:target|name|repeat)\s*=?\s*"?([^";]+)"?', re.I)

def parse_rm_attrs(attr):
    val_class=val_family=val_target=None
    m=RE_CLASS.search(attr);  val_class  = m.group(1).strip() if m else None
    m=RE_FAMILY.search(attr); val_family = m.group(1).strip() if m else None
    m=RE_TARGET.search(attr); val_target = m.group(1).strip() if m else None
    clazz,family=None,None
    if val_class:
        if "/" in val_class:
            clazz,f2 = val_class.split("/",1); clazz=clazz.strip(); family=(val_family or f2.strip()) or None
        else:
            clazz=val_class.strip(); family=(val_family or clazz) or None
    else:
        if val_target:
            if "/" in val_target:
                clazz,f2 = val_target.split("/",1); clazz=clazz.strip(); family=(val_family or f2.strip() or clazz) or None
            else:
                clazz=val_target.strip(); family=val_family or clazz
        else:
            clazz=(val_family or "NA"); family=(val_family or "NA")
    return (clazz or "NA"), (family or "NA"), (val_target or "NA")

def make_te_bed_from_out(rm_out_path):
    feats=[]
    with open(rm_out_path) as fh:
        for line in fh:
            if not line.strip(): continue
            ls=line.lstrip()
            if ls.startswith("SW") or ls.startswith("score"): continue
            parts=line.strip().split()
            if not parts or not parts[0].isdigit(): continue
            if len(parts) < 15: continue
            chrom = parts[4]
            try:
                qbeg = int(parts[5]); qend=int(parts[6])
            except: continue
            if qbeg>qend: qbeg,qend=qend,qbeg
            start0 = max(0, qbeg-1); end1=qend
            strand = parts[8]
            if strand=="C": strand="-"
            elif strand=="+": strand="+"
            else: strand="."
            cf = parts[10]
            clazz,family = parse_class_family(cf)
            te_name = f"{family}|{clazz}"
            feats.append("\t".join([chrom,str(start0),str(end1),te_name,"0",strand]))
    return pybedtools.BedTool("\n".join(feats), from_string=True) if feats else None

def make_te_bed_from_gff(rm_gff):
    feats=[]
    with open(rm_gff) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"): continue
            parts=line.rstrip("\n").split("\t")
            if len(parts)<9: continue
            chrom,src,ftype,start,end,score,strand,phase,attrs = parts
            try:
                s=int(start)-1; e=int(end)
            except: continue
            clazz,family,_ = parse_rm_attrs(attrs)
            te_name = f"{family}|{clazz}"
            feats.append("\t".join([chrom,str(s),str(e),te_name,"0",strand if strand in "+-" else "."]))
    return pybedtools.BedTool("\n".join(feats), from_string=True) if feats else None

def build_te_bed(gff_list, out_list):
    bt_gffs=[]
    for g in gff_list or []:
        bt=make_te_bed_from_gff(g)
        if bt is not None: bt_gffs.append(bt)
    bt_outs=[]
    for o in out_list or []:
        bt=make_te_bed_from_out(o)
        if bt is not None: bt_outs.append(bt)
    bt_gff=None
    if bt_gffs:
        bt_gff=bt_gffs[0]
        for extra in bt_gffs[1:]: bt_gff=bt_gff.cat(extra, postmerge=False)
        bt_gff=bt_gff.sort()
    bt_out=None
    if bt_outs:
        bt_out=bt_outs[0]
        for extra in bt_outs[1:]: bt_out=bt_out.cat(extra, postmerge=False)
        bt_out=bt_out.sort()
    if bt_gff is not None and bt_out is not None:
        gff_minus = bt_gff.subtract(bt_out)
        return bt_out.cat(gff_minus, postmerge=False).sort()
    return (bt_gff or bt_out).sort() if (bt_gff or bt_out) is not None else None

# --------------------------- plotting (R) ---------------------------

R_BAR_FAM = r"""suppressMessages({
  if (!require("ggplot2")) install.packages("ggplot2", repos="https://cloud.r-project.org", dependencies=TRUE)
  if (!require("dplyr"))   install.packages("dplyr",   repos="https://cloud.r-project.org", dependencies=TRUE)
  if (!require("scales"))  install.packages("scales",  repos="https://cloud.r-project.org", dependencies=TRUE)
  if (!require("viridis")) install.packages("viridis", repos="https://cloud.r-project.org", dependencies=TRUE)
  library(ggplot2); library(dplyr); library(scales); library(viridis)
})
in_tsv <- "__IN__"; out_pdf <- "__OUT__"; ylabel <- "__YL__"; frac <- __FRAC__;
df <- read.delim(in_tsv, sep="\t", header=TRUE, stringsAsFactors=FALSE)
if (nrow(df)==0) quit(save="no", status=0)
df$family <- factor(df$family, levels=df$family[order(df$value)])
p <- ggplot(df, aes(x=family, y=value, fill=family)) +
  geom_col() + coord_flip() +
  scale_fill_viridis_d(guide="none") +
  labs(x="TE family", y=ylabel, title="Top TE families (genome-wide)") +
  theme_minimal(base_size=12) +
  theme(axis.title.y=element_text(face="bold"), axis.title.x=element_text(face="bold"))
if (frac==1) p <- p + scale_y_continuous(labels=scales::percent)
ggsave(out_pdf, p, width=8, height=6, dpi=300, device=cairo_pdf)
"""

R_DONUT_CLASS = r"""suppressMessages({
  if (!require("ggplot2")) install.packages("ggplot2", repos="https://cloud.r-project.org", dependencies=TRUE)
  if (!require("dplyr"))   install.packages("dplyr",   repos="https://cloud.r-project.org", dependencies=TRUE)
  if (!require("scales"))  install.packages("scales",  repos="https://cloud.r-project.org", dependencies=TRUE)
  if (!require("viridis")) install.packages("viridis", repos="https://cloud.r-project.org", dependencies=TRUE)
  library(ggplot2); library(dplyr); library(scales); library(viridis)
})
in_tsv <- "__IN__"; out_pdf <- "__OUT__";
df <- read.delim(in_tsv, sep="\t", header=TRUE, stringsAsFactors=FALSE)
if (nrow(df)==0) quit(save="no", status=0)
df <- df %>% mutate(label=paste0(class, " (", round(100*frac_te,2), "%)"))
p <- ggplot(df, aes(x=3, y=frac_te, fill=reorder(class, -frac_te))) +
  geom_bar(stat="identity", width=0.3, color="white") +
  coord_polar(theta="y", start=0) + theme_void() +
  theme(legend.position="right", legend.title=element_text(face="bold")) +
  geom_text(aes(x=3.5, label=paste0(round(100*frac_te,2), "%")),
            position=position_stack(vjust=0.5), size=5, color="black") +
  scale_fill_viridis_d(name="TE class") + xlim(1,4)
ggsave(out_pdf, p, width=8, height=6, dpi=300, device=cairo_pdf)
"""

def run_r_bar_family(in_tsv, out_pdf, ylabel, is_fraction):
    rscript = shutil.which("Rscript")
    script_path = f"{out_pdf}.make.R"
    txt = (R_BAR_FAM
           .replace("__IN__", os.path.abspath(in_tsv))
           .replace("__OUT__", os.path.abspath(out_pdf))
           .replace("__YL__", ylabel.replace('"','\\"'))
           .replace("__FRAC__", "1" if is_fraction else "0"))
    with open(script_path,"w") as fh: fh.write(txt)
    if rscript:
        try: subprocess.check_call([rscript, script_path])
        except subprocess.CalledProcessError: sys.stderr.write("[WARN] R falhou no barplot famílias.\n")
    else:
        sys.stderr.write(f"[WARN] Rscript não encontrado. Para plotar: Rscript {script_path}\n")

def run_r_donut_class(in_tsv, out_pdf):
    rscript = shutil.which("Rscript")
    script_path = f"{out_pdf}.make.R"
    txt = (R_DONUT_CLASS
           .replace("__IN__", os.path.abspath(in_tsv))
           .replace("__OUT__", os.path.abspath(out_pdf)))
    with open(script_path,"w") as fh: fh.write(txt)
    if rscript:
        try: subprocess.check_call([rscript, script_path])
        except subprocess.CalledProcessError: sys.stderr.write("[WARN] R falhou no donut por classe.\n")
    else:
        sys.stderr.write(f"[WARN] Rscript não encontrado. Para plotar: Rscript {script_path}\n")

# --------------------------- main ---------------------------

def main():
    ap = argparse.ArgumentParser(description="Measure genome-wide repeat content from RepeatMasker.")
    ap.add_argument("--chrom-sizes", required=True,
                    help="TSV with chrom and length (e.g., genome.fa.fai first 2 cols)")
    ap.add_argument("--rm-gff", nargs="+", help="RepeatMasker GFF (GFF2/GFF3)", default=[])
    ap.add_argument("--rm-out", nargs="+", help="RepeatMasker .out tables", default=[])
    ap.add_argument("--top-families", type=int, default=20, help="Top N families for the barplot.")
    ap.add_argument("--family-metric", choices=["bp_sum","frac_genome","frac_te"], default="frac_te",
                    help="Ranking/plot metric: bp_sum | frac_genome | frac_te")
    ap.add_argument("--out-prefix", default="genome_repeat")
    args = ap.parse_args()

    if not args.rm_gff and not args.rm_out:
        sys.stderr.write("[ERRO] passe pelo menos um --rm-gff e/ou --rm-out\n"); sys.exit(2)

    sizes_df = read_chrom_sizes(args.chrom_sizes)
    chrom_lengths = dict(zip(sizes_df["chrom"], sizes_df["length"]))
    genome_len = int(sizes_df["length"].sum())
    sys.stderr.write(f"[INFO] Genome length: {genome_len:,} bp; chroms: {len(chrom_lengths)}\n")

    # build TE bed
    te_bt = build_te_bed(args.rm_gff, args.rm_out)
    if te_bt is None:
        raise SystemExit("[ERRO] não consegui montar BedTool de TE.")
    te_bt = te_bt.sort()

    # per-chrom stats: merge to avoid double-count
    per_chrom = []
    te_total_bp = 0
    for chrom, L in chrom_lengths.items():
        chrom_bt = pybedtools.BedTool("\n".join([f"{chrom}\t0\t{L}\t{chrom}"]), from_string=True)
        te_on_chrom = te_bt.intersect(chrom_bt, wa=True).sort().merge()  # merge on this chrom
        bp = sum(f.end - f.start for f in te_on_chrom)
        te_total_bp += bp
        per_chrom.append((chrom, L, bp, bp/L if L>0 else 0.0))
    per_chrom_df = pd.DataFrame(per_chrom, columns=["chrom","length","te_bp","te_frac"]).sort_values("chrom")

    te_frac_global = te_total_bp / genome_len if genome_len>0 else 0.0
    summary_df = pd.DataFrame([{
        "genome_length_bp": genome_len,
        "te_total_bp": te_total_bp,
        "te_fraction": round(te_frac_global, 6)
    }])

    # class/family aggregation — use original intervals (not merged) for load accounting
    by_family_bp = defaultdict(int)
    by_class_bp  = defaultdict(int)
    for f in te_bt:
        name = f.name or "Unknown|Unknown"
        fam, clazz = (name.split("|",1)+["NA"])[:2] if "|" in name else (name, "NA")
        bp = int(f.end - f.start)
        by_family_bp[fam] += bp
        by_class_bp[clazz] += bp

    fam_df = (pd.DataFrame([{"family":k, "te_bp":v} for k,v in by_family_bp.items()])
              .sort_values("te_bp", ascending=False))
    fam_df["frac_te"]     = fam_df["te_bp"] / te_total_bp if te_total_bp>0 else 0.0
    fam_df["frac_genome"] = fam_df["te_bp"] / genome_len   if genome_len>0   else 0.0

    class_df = (pd.DataFrame([{"class":k, "te_bp":v} for k,v in by_class_bp.items()])
                .sort_values("te_bp", ascending=False))
    class_df["frac_te"]     = class_df["te_bp"] / te_total_bp if te_total_bp>0 else 0.0
    class_df["frac_genome"] = class_df["te_bp"] / genome_len   if genome_len>0   else 0.0

    # write TSVs
    prefix = args.out_prefix
    summary_df.to_csv(f"{prefix}.summary.tsv", sep="\t", index=False)
    per_chrom_df.to_csv(f"{prefix}.per_chrom.tsv", sep="\t", index=False)
    class_df.to_csv(f"{prefix}.te_by_class.tsv", sep="\t", index=False)
    fam_df.to_csv(f"{prefix}.te_by_family.tsv", sep="\t", index=False)
    sys.stderr.write("[OK] TSVs escritos.\n")

    # choose metric & top families table
    if args.family_metric == "bp_sum":
        fam_plot = fam_df.assign(value=fam_df["te_bp"])
        ylabel = "Total TE bp"
        is_frac = False
    elif args.family_metric == "frac_genome":
        fam_plot = fam_df.assign(value=fam_df["frac_genome"])
        ylabel = "Share of genome"
        is_frac = True
    else:  # frac_te
        fam_plot = fam_df.assign(value=fam_df["frac_te"])
        ylabel = "Share of TE load"
        is_frac = True

    fam_top = fam_plot.sort_values("value", ascending=False).head(int(args.top_families))
    fam_top[["family","value","te_bp","frac_te","frac_genome"]].to_csv(f"{prefix}.top_families.tsv", sep="\t", index=False)

    # plots (PDF) via R, se disponível
    run_r_bar_family(f"{prefix}.top_families.tsv", f"{prefix}.top_families.pdf", ylabel, is_frac)
    class_for_plot = class_df[["class","frac_te"]].rename(columns={"class":"class","frac_te":"frac_te"})
    class_for_plot.to_csv(f"{prefix}.class_for_donut.tsv", sep="\t", index=False)
    run_r_donut_class(f"{prefix}.class_for_donut.tsv", f"{prefix}.class_donut.pdf")

    # resumo no stderr
    sys.stderr.write(f"[RESUMO] TE fraction (global): {te_frac_global:.4%}  (TE bp: {te_total_bp:,} / {genome_len:,})\n")
    sys.stderr.write("[OK] Plots finalizados (se Rscript disponível).\n")

if __name__ == "__main__":
    main()
