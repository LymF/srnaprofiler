#!/usr/bin/env python3
"""
weblogo.py
------------------------------------------
Generate WebLogo PDFs for each read length in [--min, --max].
If input is SAM, split by strand (sense vs antisense).
If input is FASTA, generate a single combined logo.
All sequences are converted from DNA (T) to RNA (U) automatically.

NEW:
  --uniq              -> remove duplicated sequences por bucket (tamanho [+ strand])
  --limit-positions N -> trunca as sequências para as primeiras N bases (padrão: 11)

Exemplos:
  python weblogo.py -i sample.sam   -o out/prefix --min 24 --max 28 --uniq --limit-positions 11
  python weblogo.py -i sample.fasta -o out/prefix --min 21 --max 22 --limit-positions 11

Dependências:
  conda install -c bioconda weblogo biopython
"""

import subprocess
import argparse
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple
from Bio import SeqIO

# ----------------------- Helper functions -----------------------

def read_sam_sequences_with_strand(sam_path: Path) -> List[Tuple[str, str]]:
    """
    Extrai (sequência, strand) do SAM.
    Converte T->U para usar alfabeto RNA.
    """
    out = []
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            try:
                flag = int(parts[1])
            except ValueError:
                continue
            seq = parts[9].strip().upper()
            if not seq or seq == "*":
                continue
            seq = seq.replace("T", "U")  # DNA -> RNA
            strand = "antisense" if (flag & 16) else "sense"
            out.append((seq, strand))
    return out

def read_fasta_sequences(fasta_path: Path) -> List[str]:
    """Lê FASTA e converte T->U."""
    return [str(r.seq).upper().replace("T","U") for r in SeqIO.parse(str(fasta_path), "fasta")]

def revcomp_rna(seq: str) -> str:
    """Reverse complement em RNA."""
    comp = str.maketrans("ACGUacguNn", "UGCAugcaNn")
    return seq.translate(comp)[::-1]

def generate_weblogo(input_fasta: Path, output_pdf: Path, title: str):
    """
    Chama o WebLogo (CLI). Sintaxe moderna do --color: COLOR REGEX DESCRIPTION
    """
    cmd = [
        "weblogo",
        "-f", str(input_fasta),
        "-o", str(output_pdf),
        "--format", "pdf",
        "--fineprint", "",
        "--title", title,
        "--color", "green", "A", "Adenine (A)",
        "--color", "blue",  "C", "Cytosine (C)",
        "--color", "orange","G", "Guanine (G)",
        "--color", "red",   "U", "Uracil (U)",
        "--size", "large",
        "--sequence-type", "rna",
        "--resolution", "300"
    ]
    subprocess.run(cmd, check=True)
    print(f"[INFO] WebLogo generated: {output_pdf}")

def uniq_list(seq_list: List[str]) -> List[str]:
    """
    Remove duplicatas preservando a primeira ocorrência (equivalente a awk '!seen[$0]++').
    """
    seen = set()
    out = []
    for s in seq_list:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out

def trim_to_first_n(seq_list: List[str], n: int) -> List[str]:
    """Trunca todas as sequências para as primeiras n bases."""
    if n is None or n <= 0:
        return seq_list
    return [s[:n] if len(s) >= n else s for s in seq_list]

# ----------------------- Main pipeline -----------------------

def main():
    ap = argparse.ArgumentParser(
        description="Generate WebLogo plots por classe de tamanho (SAM ou FASTA)."
    )
    ap.add_argument("-i", "--input", required=True, help="Arquivo SAM ou FASTA")
    ap.add_argument("-o", "--outprefix", required=True, help="Prefixo de saída para arquivos")
    ap.add_argument("--min", type=int, default=20, help="Tamanho mínimo (nt)")
    ap.add_argument("--max", type=int, default=30, help="Tamanho máximo (nt)")
    ap.add_argument("--title", default="Base composition across small RNAs", help="Título base dos plots")
    ap.add_argument("--rc-antisense", action="store_true",
                    help="Faz reverse-complement das reads antisense antes do WebLogo (SAM apenas)")
    ap.add_argument("--uniq", action="store_true",
                    help="Remove sequências duplicadas por bucket (tamanho [+ strand])")
    ap.add_argument("--limit-positions", type=int, default=11,
                    help="Truncar sequências para as primeiras N bases (padrão: 11)")
    args = ap.parse_args()

    in_path = Path(args.input)
    outprefix = Path(args.outprefix)
    outprefix.parent.mkdir(parents=True, exist_ok=True)

    # ----------- Entrada SAM -----------
    if in_path.suffix.lower() == ".sam":
        print(f"[INFO] Reading SAM: {in_path}")
        seqs = read_sam_sequences_with_strand(in_path)
        by_size_strand = defaultdict(lambda: {"sense": [], "antisense": []})
        for seq, strand in seqs:
            L = len(seq)
            if args.min <= L <= args.max:
                by_size_strand[L][strand].append(seq)

        total_logos = 0
        for L in range(args.min, args.max + 1):
            group = by_size_strand.get(L, None)
            if not group:
                print(f"[WARN] No reads of length {L} nt.")
                continue

            sense_reads = group["sense"]
            antisense_reads = group["antisense"]

            # --uniq: deduplicar por bucket
            if args.uniq:
                before_sense = len(sense_reads)
                before_anti  = len(antisense_reads)
                sense_reads = uniq_list(sense_reads)
                antisense_reads = uniq_list(antisense_reads)
                removed_sense = before_sense - len(sense_reads)
                removed_anti  = before_anti  - len(antisense_reads)
                if removed_sense or removed_anti:
                    print(f"[INFO] {L} nt: removed duplicates -> sense -{removed_sense}, antisense -{removed_anti}")

            # RC no antisense se pedido
            if args.rc_antisense and antisense_reads:
                antisense_reads = [revcomp_rna(s) for s in antisense_reads]

            # TRUNCAR para as primeiras N posições
            sense_reads     = trim_to_first_n(sense_reads, args.limit_positions)
            antisense_reads = trim_to_first_n(antisense_reads, args.limit_positions)

            # FASTA combinado
            combined_fasta = outprefix.parent / f"{outprefix.name}.size_{L}.combined.first{args.limit_positions}.fasta"
            with open(combined_fasta, "w") as fh:
                i = 1
                for s in sense_reads:
                    fh.write(f">sense_{i}\n{s}\n"); i += 1
                for s in antisense_reads:
                    fh.write(f">antisense_{i}\n{s}\n"); i += 1

            if len(sense_reads) + len(antisense_reads) >= 2:
                pdf = outprefix.parent / f"{outprefix.name}.size_{L}.combined.first{args.limit_positions}.weblogo.pdf"
                try:
                    generate_weblogo(combined_fasta, pdf, f"{args.title} ({L} nt, combined; first {args.limit_positions} nt)")
                    total_logos += 1
                except subprocess.CalledProcessError as e:
                    print(f"[ERROR] WebLogo failed (combined {L}): {e}")

            # Sense-only
            if len(sense_reads) >= 2:
                sense_fa = outprefix.parent / f"{outprefix.name}.size_{L}.sense.first{args.limit_positions}.fasta"
                with open(sense_fa, "w") as fh:
                    for i, s in enumerate(sense_reads, 1):
                        fh.write(f">sense_{i}\n{s}\n")
                pdf = outprefix.parent / f"{outprefix.name}.size_{L}.sense.first{args.limit_positions}.weblogo.pdf"
                try:
                    generate_weblogo(sense_fa, pdf, f"{args.title} ({L} nt, sense; first {args.limit_positions} nt)")
                    total_logos += 1
                except subprocess.CalledProcessError as e:
                    print(f"[ERROR] WebLogo failed (sense {L}): {e}")
            else:
                print(f"[WARN] Not enough sense reads for {L} nt.")

            # Antisense-only
            if len(antisense_reads) >= 2:
                anti_fa = outprefix.parent / f"{outprefix.name}.size_{L}.antisense.first{args.limit_positions}.fasta"
                with open(anti_fa, "w") as fh:
                    for i, s in enumerate(antisense_reads, 1):
                        fh.write(f">antisense_{i}\n{s}\n")
                pdf = outprefix.parent / f"{outprefix.name}.size_{L}.antisense.first{args.limit_positions}.weblogo.pdf"
                try:
                    generate_weblogo(anti_fa, pdf, f"{args.title} ({L} nt, antisense{' RC' if args.rc_antisense else ''}; first {args.limit_positions} nt)")
                    total_logos += 1
                except subprocess.CalledProcessError as e:
                    print(f"[ERROR] WebLogo failed (antisense {L}): {e}")
            else:
                print(f"[WARN] Not enough antisense reads for {L} nt.")

        if total_logos == 0:
            print("[WARN] No WebLogos generated. Check size range or input.")

    # ----------- Entrada FASTA -----------
    elif in_path.suffix.lower() in [".fa", ".fasta"]:
        print(f"[INFO] Reading FASTA: {in_path}")
        seqs = read_fasta_sequences(in_path)

        # --uniq global (FASTA)
        if args.uniq:
            before = len(seqs)
            seqs = uniq_list(seqs)
            removed = before - len(seqs)
            if removed:
                print(f"[INFO] Removed {removed} duplicate sequences (global FASTA).")

        by_size = defaultdict(list)
        for s in seqs:
            L = len(s)
            if args.min <= L <= args.max:
                by_size[L].append(s)

        total_logos = 0
        for L in range(args.min, args.max + 1):
            reads = by_size.get(L, [])
            if len(reads) < 2:
                print(f"[WARN] Not enough reads of length {L} nt (need ≥2).")
                continue

            # truncar para as primeiras N posições
            reads = trim_to_first_n(reads, args.limit_positions)

            fa = outprefix.parent / f"{outprefix.name}.size_{L}.combined.first{args.limit_positions}.fasta"
            with open(fa, "w") as fh:
                for i, s in enumerate(reads, 1):
                    fh.write(f">seq_{i}\n{s}\n")
            pdf = outprefix.parent / f"{outprefix.name}.size_{L}.combined.first{args.limit_positions}.weblogo.pdf"
            try:
                generate_weblogo(fa, pdf, f"{args.title} ({L} nt; first {args.limit_positions} nt)")
                print(f"[INFO] {L} nt: {len(reads)} reads → {pdf.name}")
                total_logos += 1
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] WebLogo failed ({L}): {e}")
        if total_logos == 0:
            print("[WARN] No WebLogos generated. Check size range or input.")
    else:
        print("[ERROR] Unsupported input format. Provide .sam or .fasta/.fa")

if __name__ == "__main__":
    main()
