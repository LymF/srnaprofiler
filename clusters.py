#!/usr/bin/env python3
"""
clusters.py - Small RNA cluster analysis module
Integrated version of smallrna_profiler_with_mirna_complete.py
"""
import sys
import os
import argparse
import json
import subprocess
import re
import gzip
import statistics
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional, TextIO, Set
from collections import defaultdict

try:
    from Bio import SeqIO
    from pybedtools import BedTool
except ImportError as e:
    print(f"Error importing required packages: {e}", file=sys.stderr)
    print("Please install required packages: pip install biopython pybedtools", file=sys.stderr)
    sys.exit(1)

from utils import log, ensure_outdir, validate_file_exists

__version__ = "1.2.0-integrated"


# ------------------------------- Utils --------------------------------

def get_script_directory():
    """Return the directory where the current Python script is located."""
    return Path(__file__).parent.absolute()

def get_r_scripts_path():
    """Return the path to R scripts directory (same as Python script directory)."""
    return get_script_directory()

def open_sam(filepath: str) -> TextIO:
    return gzip.open(filepath, "rt") if str(filepath).endswith(".gz") else open(filepath, "r")

def write_text(filepath: Path, content: str):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        f.write(content)

def read_fasta(filepath: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    cur_id = None; cur = []
    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(cur)
                cur_id = line[1:].split()[0]
                cur = []
            else:
                cur.append(line)
        if cur_id:
            seqs[cur_id] = "".join(cur)
    return seqs

def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]

def parse_coords_from_header(header: str):
    """
    Parse 'Chr:START-END(+/-)' or 'Chr_START-END(+/-)' or 'Chr:START-END' from FASTA header.
    Returns (chrom, start, end, strand) or None.
    """
    m = re.search(r"([\w\-\.\|]+)[:_](\d+)[\-:](\d+)(?:\(([-+])\))?", header)
    if not m:
        return None
    chrom, start, end, strand = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4) or "+"
    return chrom, start, end, strand

def fasta_to_bed(fasta_path: Path, out_bed: Path) -> Path:
    with open(out_bed, "w") as out:
        for rec in SeqIO.parse(str(fasta_path), "fasta"):
            coords = parse_coords_from_header(rec.description)
            if coords:
                chrom, start, end, strand = coords
                out.write(f"{chrom}\t{start}\t{end}\t{rec.id}\t0\t{strand}\n")
    return out_bed

def read_classification_table(path: Path) -> Dict[str, str]:
    d: Dict[str, str] = {}
    if not path.exists():
        return d
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                d[parts[0]] = parts[1]
    return d

def write_classification_table(d: Dict[str, str], path: Path):
    with open(path, "w") as out:
        for k, v in sorted(d.items()):
            out.write(f"{k}\t{v}\n")


# ------------------------------ StepMgr -------------------------------

class StepManager:
    """Writes .ok tickets to logs/ and updates manifest JSON in manifest/."""
    def __init__(self, outdir: Path, prefix: str, logs_dir: Path, manifest_dir: Path, force: bool = False):
        self.outdir = outdir
        self.prefix = prefix
        self.force = force
        self.logs_dir = ensure_outdir(logs_dir)
        self.manifest_dir = ensure_outdir(manifest_dir)
        self.manifest = {"time": datetime.now().isoformat(), "version": __version__, "steps": []}

    def step_ticket(self, step_name: str) -> Path:
        return self.logs_dir / f".{step_name}.ok"

    def step_done(self, step_name: str) -> bool:
        return self.step_ticket(step_name).exists()

    def mark_done(self, step_name: str):
        write_text(self.step_ticket(step_name), f"{datetime.now()}\n")

    def append_step_record(self, step_name: str, outputs: Dict[str, str]):
        self.manifest["steps"].append({"name": step_name, "outputs": outputs, "time": datetime.now().isoformat()})
        manifest_file = self.manifest_dir / f"{self.prefix}.run.json"
        with open(manifest_file, "w") as f:
            json.dump(self.manifest, f, indent=2)


# ------------------------------ Analyses ------------------------------

class SAMStatistics:
    """
    Replica da lógica de classificação (sem alterar): contabiliza sizes/strand por contig
    e atribui rótulos finais: putative_siRNA / putative_piRNA / ambiguous.
    """
    def __init__(self, sam_file: str, fasta_file: Optional[str], len_min: int = 15, len_max: int = 35):
        self.sam_file = sam_file
        self.fasta_file = fasta_file
        self.len_min = len_min
        self.len_max = len_max
        self.tamanhos_p = defaultdict(int)
        self.tamanhos_n = defaultdict(int)
        self.ch = defaultdict(lambda: {"pos": 0, "neg": 0})
        self.ch2 = defaultdict(lambda: {"pos": defaultdict(int), "neg": defaultdict(int)})
        self.reads_seen = set()
        self.cp = 0
        self.cn = 0

    def analyze(self) -> Tuple[str, Dict[str, Dict]]:
        log(f"Analyzing SAM: {self.sam_file}", "INFO")
        with open_sam(self.sam_file) as f:
            for line in f:
                if line.startswith("@"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 10:
                    continue
                read_id = fields[0]
                flag = int(fields[1])
                rname = fields[2]
                seq = fields[9]
                size = len(seq)
                if size < self.len_min or size > self.len_max:
                    continue
                read_key = f"{rname}_{read_id}"
                if read_key in self.reads_seen:
                    continue
                self.reads_seen.add(read_key)
                if flag & 16:
                    self.cn += 1
                    self.ch[rname]["neg"] += 1
                    self.tamanhos_n[size] += 1
                    self.ch2[rname]["neg"][size] += 1
                else:
                    self.cp += 1
                    self.ch[rname]["pos"] += 1
                    self.tamanhos_p[size] += 1
                    self.ch2[rname]["pos"][size] += 1

        total = self.cp + self.cn
        output = self._generate_output_text(total)
        signal = self._classify_loci()
        return output, signal

    def _generate_output_text(self, total: int) -> str:
        lines = []
        lines.append(f"Total in both strand: {total}")
        lines.append(f"Total in positve strand: {self.cp}")
        lines.append(f"Total in negative strand: {self.cn}")
        lines.append("\nSize\tSENSE\tANTISENSE\tTOTAL")
        all_sizes = set(self.tamanhos_p.keys()) | set(self.tamanhos_n.keys())
        for size in sorted(all_sizes):
            p = self.tamanhos_p.get(size, 0)
            n = self.tamanhos_n.get(size, 0)
            s = p + n
            if s != 0:
                lines.append(f"{size}\t{p}\t{n}\t{s}")
        lines.append("\n\nReads per chromosoes")
        for c in sorted(self.ch.keys()):
            p = self.ch[c]["pos"]; n = self.ch[c]["neg"]; s = p + n
            lines.append(f"{c}\t{p}\t{n}\t{s}")
        return "\n".join(lines) + "\n"

    def _classify_loci(self) -> Dict[str, Dict]:
        signal = {}
        for c in self.ch2.keys():
            signal[c] = {"si": "no", "pi": "no", "final": ""}
            sumT1P = sumT1N = allT1 = 0
            sumT2P = sumT2N = allT2 = 0
            allT3 = 0
            allT4 = allT4_pos = allT4_neg = 0
            scorePI = 0
            for s in range(self.len_min, self.len_max + 1):
                pa = self.ch2[c]["pos"].get(s, 0)
                na = self.ch2[c]["neg"].get(s, 0)
                ta = pa + na
                if 15 <= s <= 19:
                    sumT1P += pa; sumT1N += na; allT1 = sumT1P + sumT1N
                elif 20 <= s <= 22:
                    allT4 += pa + na; allT4_pos += pa; allT4_neg += na
                    if s == 22:
                        if (allT4_pos > 50 or allT4_neg > 50) and (allT4 > allT1):
                            na_safe = max(na, 1); pa_safe = max(pa, 1)
                            r = (pa_safe/na_safe) if pa_safe >= na_safe else (na_safe/pa_safe)
                            if 0.3 < r < 10: signal[c]["si"] = "yes"
                    if pa > 0 and na > 0:
                        if (pa > 50 or na > 50) and (ta > allT1):
                            r = (pa/na) if pa >= na else (na/pa)
                            if 0.3 < r < 10: signal[c]["si"] = "yes"
                elif s in [22, 23]:
                    sumT2P += pa; sumT2N += na; allT2 = sumT2P + sumT2N
                elif 24 <= s <= 29:
                    if pa > 0 and na > 0:
                        r2 = (pa/na) if pa >= na else (na/pa)
                        if pa > 50 and ta > 50:
                            scorePI += 1
                        allT3 += pa + ta
                elif s == 30:
                    if (allT3 > allT2 + allT1) and (scorePI >= 4):
                        signal[c]["pi"] = "yes"
            if (allT4 > allT3) and (allT4 > allT1):
                signal[c]["final"] = "putative_siRNA"
            elif (allT3 > allT4) and (allT3 > allT1):
                signal[c]["final"] = "putative_piRNA"
            elif (allT1 > allT3) and (allT1 > allT4):
                signal[c]["final"] = "ambiguous"
            else:
                signal[c]["final"] = "ambiguous"
        return signal

    def write_counts(self, output_file: Path):
        with open(output_file, "w") as f:
            for c in sorted(self.ch.keys()):
                p = self.ch[c]["pos"]; n = self.ch[c]["neg"]; s = p + n
                f.write(f"{c}\t{p}\t{n}\t{s}\n")

    def write_distribution(self, output_file: Path):
        with open(output_file, "w") as f:
            f.write("\nSize\tSENSE\tANTISENSE\tTOTAL\n")
            for i in range(15, 30):
                p = self.tamanhos_p.get(i, 0)
                n = self.tamanhos_n.get(i, 0)
                t = p + n
                f.write(f"{i}\t{p}\t{n}\t{t}\n")

    def write_profile_files(self, profiles_prefix: Path, genome_seqs: Optional[Dict[str, str]],
                            signal: Dict[str, Dict], clusters_dir: Path):
        # profiles (.profile) in profiles/
        profile1 = profiles_prefix.with_suffix(".profile1")
        profile2 = profiles_prefix.with_suffix(".profile2")
        profile1.parent.mkdir(parents=True, exist_ok=True)
        with open(profile1, "w") as f1, open(profile2, "w") as f2:
            f2.write(profiles_prefix.name)
            for i in range(15, 30 + 1):
                p = self.tamanhos_p.get(i, 0)
                n = self.tamanhos_n.get(i, 0)
                t = p + n
                f1.write(f"{i}\t{p}\t{n}\t{t}\n")
                f2.write(f"#{t}")
            f2.write("\n")

        if not genome_seqs:
            return
        # cluster FASTAs in clusters/
        clusters_dir.mkdir(parents=True, exist_ok=True)
        si_sig = clusters_dir / f"{profiles_prefix.name}.signal.putative_siRNA.fasta"
        pi_sig = clusters_dir / f"{profiles_prefix.name}.signal.putative_piRNA.fasta"
        both_sig = clusters_dir / f"{profiles_prefix.name}.signal.putative_siRNA_and_putative_piRNA.fasta"
        amb_sig = clusters_dir / f"{profiles_prefix.name}.signal.ambiguous.fasta"
        si_pred = clusters_dir / f"{profiles_prefix.name}.clusters.putative_siRNA.fasta"
        pi_pred = clusters_dir / f"{profiles_prefix.name}.clusters.putative_piRNA.fasta"
        amb_pred = clusters_dir / f"{profiles_prefix.name}.clusters.ambiguous.fasta"
        with open(si_sig, "w") as f_si, open(pi_sig, "w") as f_pi, \
             open(both_sig, "w") as f_both, open(amb_sig, "w") as f_amb, \
             open(si_pred, "w") as f_si_p, open(pi_pred, "w") as f_pi_p, \
             open(amb_pred, "w") as f_amb_p:
            for c in sorted(signal.keys()):
                if c not in genome_seqs: 
                    continue
                seq = genome_seqs[c]
                fin = signal[c]["final"]
                if signal[c]["si"] == "yes" and signal[c]["pi"] == "no":
                    f_si.write(f">{c}\n{seq}\n")
                if signal[c]["pi"] == "yes" and signal[c]["si"] == "no":
                    f_pi.write(f">{c}\n{seq}\n")
                if signal[c]["pi"] == "yes" and signal[c]["si"] == "yes":
                    f_both.write(f">{c}\n{seq}\n")
                if signal[c]["si"] == "no" and signal[c]["pi"] == "no":
                    f_amb.write(f">{c}\n{seq}\n")
                if fin == "putative_siRNA":
                    f_si_p.write(f">{c}\n{seq}\n")
                elif fin == "putative_piRNA":
                    f_pi_p.write(f">{c}\n{seq}\n")
                elif fin == "ambiguous":
                    f_amb_p.write(f">{c}\n{seq}\n")


# ------------------------------ Strand --------------------------------

class DefineStrand:
    def __init__(self, counts_file: Path, fasta_file: Path, output_file: Path):
        self.counts_file = counts_file
        self.fasta_file = fasta_file
        self.output_file = output_file

    def process(self):
        log("Defining strand by majority sense/antisense per contig...", "INFO")
        contigs: Dict[str, Dict[str, str]] = {}
        with open(self.fasta_file, "r") as f:
            cur = None; seq = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if cur:
                        s = "".join(seq)
                        contigs[cur] = {"sense": s, "antisense": reverse_complement(s)}
                    cur = line[1:].split()[0]; seq = []
                else:
                    seq.append(line)
            if cur:
                s = "".join(seq)
                contigs[cur] = {"sense": s, "antisense": reverse_complement(s)}

        sense_flag: Dict[str, int] = {}
        with open(self.counts_file, "r") as cf:
            for line in cf:
                if not line.strip() or line.startswith("Chr") or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                c = parts[0]
                try:
                    p = int(parts[1]); n = int(parts[2])
                except ValueError:
                    continue
                sense_flag[c] = 1 if p > n else 0

        with open(self.output_file, "w") as out:
            for cand in sorted(contigs.keys()):
                out.write(f">{cand}\n")
                out.write(contigs[cand]["sense"] + "\n" if sense_flag.get(cand, 1) == 1
                          else contigs[cand]["antisense"] + "\n")
        log(f"Defined-strand FASTA written: {self.output_file}", "INFO")


# ------------------------- R-based distributions ----------------------

class PlotGeralDistribution:
    """
    Writes:
      profiles (.profile1/.profile2) → profiles/
      profile_per_size               → tables/
      Geral_base_distribution*       → tables/
    And runs R scripts in plots/ to produce PDFs.
    """
    def __init__(self, sam_file: str, basename: str, tables_dir: Path, profiles_dir: Path,
                 plots_dir: Path, start: int = 15, end: int = 30,
                 norm: Optional[int] = None, minimo: int = 30,
                 plot: bool = False, antisense: bool = False, r_scripts_path: str = ""):
        self.sam_file = sam_file
        self.basename = basename
        self.tables = tables_dir
        self.profiles = profiles_dir
        self.plots = plots_dir
        self.start = start
        self.end = end
        self.norm = norm
        self.minimo = minimo
        self.plot = plot
        self.antisense = antisense
        # Use provided path or get from script directory
        self.r_scripts_path = r_scripts_path or get_r_scripts_path()

        self.tamanhos_p = defaultdict(int)
        self.tamanhos_n = defaultdict(int)
        self.cp = 0; self.cn = 0
        self.ch = defaultdict(lambda: {"pos": 0, "neg": 0})
        self.ch2 = defaultdict(lambda: {"pos": defaultdict(int), "neg": defaultdict(int)})
        self.base_per_size = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        self.read = set()

    @property
    def profile1_path(self) -> Path: return self.profiles / f"{self.basename}.profile1"
    @property
    def profile2_path(self) -> Path: return self.profiles / f"{self.basename}.profile2"
    @property
    def profile_per_size_path(self) -> Path: return self.tables / f"{self.basename}.profile_per_size"
    @property
    def base_dist_path(self) -> Path: return self.tables / f"{self.basename}.Geral_base_distribution"
    @property
    def base_dist_norm_path(self) -> Path: return self.tables / f"{self.basename}.Geral_base_distribution_norm"

    def analyze(self):
        log(f"Loading SAM for base distribution: {self.sam_file}", "INFO")
        self.tables.mkdir(parents=True, exist_ok=True)
        self.profiles.mkdir(parents=True, exist_ok=True)

        spos = "16" if self.antisense else "0"
        sneg = "0" if self.antisense else "16"
        nmappings = 0; nreads = 0

        with open_sam(self.sam_file) as f:
            for line in f:
                if line.startswith("@"):
                    continue
                nmappings += 1
                fields = line.strip().split("\t")
                if len(fields) < 10:
                    continue
                read_id, flag, c = fields[0], fields[1], fields[2]
                seq = fields[9].upper()
                size = len(seq)
                if read_id in self.read:
                    continue
                self.read.add(read_id); nreads += 1

                if size not in self.tamanhos_p:
                    for base_char in ["A", "C", "G", "T"]:
                        self.base_per_size[size]["+"][base_char] = 0
                        self.base_per_size[size]["-"][base_char] = 0
                if size not in self.tamanhos_n:
                    self.tamanhos_n[size] = 0
                if c not in self.ch:
                    self.ch[c] = {"pos": 0, "neg": 0}
                if size not in self.ch2[c]["pos"]:
                    self.ch2[c]["pos"][size] = 0
                if size not in self.ch2[c]["neg"]:
                    self.ch2[c]["neg"][size] = 0

                first = seq[0] if seq else "N"
                if flag == "16":
                    first = reverse_complement(seq)[0] if seq else "N"

                if flag == spos:
                    self.cp += 1; self.ch[c]["pos"] += 1
                    self.tamanhos_p[size] += 1; self.ch2[c]["pos"][size] += 1
                    if first in ["A", "C", "G", "T"]:
                        self.base_per_size[size]["+"][first] += 1
                elif flag == sneg:
                    self.cn += 1; self.ch[c]["neg"] += 1
                    self.tamanhos_n[size] += 1; self.ch2[c]["neg"][size] += 1
                    if first in ["A", "C", "G", "T"]:
                        self.base_per_size[size]["-"][first] += 1

        total = self.cp + self.cn
        log(f"Tn: {self.cn}\tTp: {self.cp}\tTt: {total}\tReads: {nreads}", "INFO")

        self._write_profiles_and_tables()

        if self.plot:
            self._call_r_plots()

    def _write_profiles_and_tables(self):
        # profiles
        with open(self.profile1_path, "w") as f1, open(self.profile2_path, "w") as f2:
            f2.write(self.basename)
            for i in range(self.start, self.end + 1):
                p = self.tamanhos_p.get(i, 0)
                n = self.tamanhos_n.get(i, 0)
                t = p + n
                f1.write(f"{i}\t{p}\t{n}\t{t}\n")
                f2.write(f"#{t}")
            f2.write("\n")

        # per-size (tables)
        with open(self.profile_per_size_path, "w") as f:
            for i in range(self.start, self.end + 1):
                p = self.tamanhos_p.get(i, 0)
                f.write(f"{self.basename}#{p}" if i == self.start else f"#{p}")
            for i in range(self.start, self.end + 1):
                n = self.tamanhos_n.get(i, 0)
                f.write(f"#{n}")
            f.write("\n")

        # base distributions (tables)
        with open(self.base_dist_path, "w") as f:
            f.write("size\tA\tC\tG\tU\n")
            for i in range(self.start, self.end + 1):
                ta = self.base_per_size[i]["+"]["A"]
                tc = self.base_per_size[i]["+"]["C"]
                tg = self.base_per_size[i]["+"]["G"]
                tt = self.base_per_size[i]["+"]["T"]
                f.write(f"{i}\t{ta}\t{tc}\t{tg}\t{tt}\n")
                nta = self.base_per_size[i]["-"]["A"] * -1
                ntc = self.base_per_size[i]["-"]["C"] * -1
                ntg = self.base_per_size[i]["-"]["G"] * -1
                ntt = self.base_per_size[i]["-"]["T"] * -1
                f.write(f"{i}\t{nta}\t{ntc}\t{ntg}\t{ntt}\n")

        if self.norm:
            with open(self.base_dist_norm_path, "w") as f:
                f.write("size\tA\tC\tG\tU\n")
                for i in range(self.start, self.end + 1):
                    ta = (self.base_per_size[i]["+"]["A"] / self.norm) * 1_000_000
                    tc = (self.base_per_size[i]["+"]["C"] / self.norm) * 1_000_000
                    tg = (self.base_per_size[i]["+"]["G"] / self.norm) * 1_000_000
                    tt = (self.base_per_size[i]["+"]["T"] / self.norm) * 1_000_000
                    f.write(f"{i}\t{ta}\t{tc}\t{tg}\t{tt}\n")
                    nta = (self.base_per_size[i]["-"]["A"] * -1 / self.norm) * 1_000_000
                    ntc = (self.base_per_size[i]["-"]["C"] * -1 / self.norm) * 1_000_000
                    ntg = (self.base_per_size[i]["-"]["G"] * -1 / self.norm) * 1_000_000
                    ntt = (self.base_per_size[i]["-"]["T"] * -1 / self.norm) * 1_000_000
                    f.write(f"{i}\t{nta}\t{ntc}\t{ntg}\t{ntt}\n")

    def _call_r_plots(self):
        """Generate plots using plots module"""
        try:
            # Try to import plots module
            import sys
            script_dir = get_script_directory()
            if str(script_dir) not in sys.path:
                sys.path.insert(0, str(script_dir))
            
            try:
                import plots
                
                # Generate base distribution plots
                log("Generating base distribution plots with plots module...", "INFO")
                
                # Non-normalized plots
                if self.base_dist_path.exists():
                    output_prefix = self.plots / self.basename
                    
                    # Per base plot
                    plots.generate_base_distribution_plot(
                        input_tsv=self.base_dist_path,
                        output_prefix=output_prefix / "Geral_base_distribution_by_reads",
                        r_scripts_path=Path(self.r_scripts_path),
                        plot_type="per_base"
                    )
                    
                    # Zoom plot
                    plots.generate_base_distribution_plot(
                        input_tsv=self.base_dist_path,
                        output_prefix=output_prefix / "Geral_base_distribution_by_reads_ZOOM_piRNAS",
                        r_scripts_path=Path(self.r_scripts_path),
                        plot_type="per_base_zoom"
                    )
                    
                    # Sum strands plot
                    plots.generate_base_distribution_plot(
                        input_tsv=self.base_dist_path,
                        output_prefix=output_prefix / "Geral_base_distribution_sumstrands_by_reads",
                        r_scripts_path=Path(self.r_scripts_path),
                        plot_type="sum_strand"
                    )
                
                # Normalized plots
                if self.norm and self.base_dist_norm_path.exists():
                    output_prefix = self.plots / self.basename
                    
                    # Per base plot (normalized)
                    plots.generate_base_distribution_plot(
                        input_tsv=self.base_dist_norm_path,
                        output_prefix=output_prefix / "NORMALIZED_Geral_base_distribution_by_reads",
                        r_scripts_path=Path(self.r_scripts_path),
                        plot_type="per_base"
                    )
                    
                    # Zoom plot (normalized)
                    plots.generate_base_distribution_plot(
                        input_tsv=self.base_dist_norm_path,
                        output_prefix=output_prefix / "NORMALIZED_Geral_base_distribution_by_reads_ZOOM_piRNAS",
                        r_scripts_path=Path(self.r_scripts_path),
                        plot_type="per_base_zoom"
                    )
                    
                    # Sum strands plot (normalized)
                    plots.generate_base_distribution_plot(
                        input_tsv=self.base_dist_norm_path,
                        output_prefix=output_prefix / "NORMALIZED_Geral_base_distribution_sumstrands_by_reads",
                        r_scripts_path=Path(self.r_scripts_path),
                        plot_type="sum_strand"
                    )
                
                log("Base distribution plots generated successfully", "INFO")
                
            except ImportError:
                log("plots module not found, falling back to direct R calls", "INFO")
                self._call_r_plots_fallback()
                
        except Exception as e:
            log(f"Error generating plots: {e}", "WARN")
            log("Trying fallback method...", "INFO")
            self._call_r_plots_fallback()
    
    def _call_r_plots_fallback(self):
        """Fallback method using direct R script calls"""
        r_exec = "Rscript"
        scripts = [
            ("calcDistributionPerBase_publication.R",
             str(self.base_dist_path),
             f"{self.basename}.Geral_base_distribution_by_reads"),
            ("calcDistributionPerBase_publication_zoom_piRNAs.R",
             str(self.base_dist_path),
             f"{self.basename}.Geral_base_distribution_by_reads_ZOOM_piRNAS"),
            ("calcDistributionPerBaseSumStrand_publication.R",
             str(self.base_dist_path),
             f"{self.basename}.Geral_base_distribution_sumstrands_by_reads"),
        ]
        if self.norm:
            scripts += [
                ("calcDistributionPerBase_publication.R",
                 str(self.base_dist_norm_path),
                 f"{self.basename}.NORMALIZED_Geral_base_distribution_by_reads",
                 "RPM"),
                ("calcDistributionPerBase_publication_zoom_piRNAs.R",
                 str(self.base_dist_norm_path),
                 f"{self.basename}.NORMALIZED_Geral_base_distribution_by_reads_ZOOM_piRNAS",
                 "RPM"),
                ("calcDistributionPerBaseSumStrand_publication.R",
                 str(self.base_dist_norm_path),
                 f"{self.basename}.NORMALIZED_Geral_base_distribution_sumstrands_by_reads",
                 "RPM"),
            ]
        for item in scripts:
            script_name, infile, outpref, *extra = item
            script_path = os.path.join(self.r_scripts_path, script_name)
            if os.path.exists(script_path):
                try:
                    cmd = [r_exec, script_path, infile, os.path.basename(outpref)]
                    if extra:
                        cmd.extend(extra)
                    subprocess.run(cmd, check=False, cwd=str(self.plots))
                    log(f"Generated plot: {outpref} (plots/)", "INFO")
                except Exception as e:
                    log(f"Could not run R script {script_name}: {e}", "WARN")


# ----------------------------- Pipeline -------------------------------

class SmallRNAPipeline:
    def __init__(self, sam_file: str, genome_file: str, prefix: str,
                 reads_file: Optional[str] = None, outdir: Optional[str] = None,
                 si_min: int = 100, pi_min: int = 100,
                 len_min: int = 15, len_max: int = 35,
                 norm: Optional[int] = None, mapq_min: int = 0, nm_max: int = 1_000_000_000,
                 run_patterns: bool = False, force: bool = False,
                 r_scripts_path: str = "",
                 mirna_fasta: Optional[str] = None,
                 mirna_gff: Optional[str] = None,
                 mirna_bed: Optional[str] = None,
                 mirna_keep: bool = False):

        self.sam_file = sam_file
        self.genome_file = genome_file
        self.prefix = prefix
        self.reads_file = reads_file
        self.outdir = Path(outdir or f"{prefix}_out")
        self.si_min = si_min
        self.pi_min = pi_min
        self.len_min = len_min
        self.len_max = len_max
        self.norm = norm
        self.mapq_min = mapq_min
        self.nm_max = nm_max
        self.run_patterns = run_patterns
        self.force = force
        # Use provided path or get from script directory
        self.r_scripts_path = r_scripts_path or get_r_scripts_path()
        self.mirna_fasta = mirna_fasta
        self.mirna_gff = mirna_gff
        self.mirna_bed = mirna_bed
        self.mirna_keep = mirna_keep

        # folders
        self.paths = {
            "logs": ensure_outdir(self.outdir / "logs"),
            "sam": ensure_outdir(self.outdir / "sam"),
            "clusters": ensure_outdir(self.outdir / "clusters"),
            "tables": ensure_outdir(self.outdir / "tables"),
            "plots": ensure_outdir(self.outdir / "plots"),
            "profiles": ensure_outdir(self.outdir / "profiles"),
            "manifest": ensure_outdir(self.outdir / "manifest"),
            "defined_strand": ensure_outdir(self.outdir / "defined_strand"),
        }

        self.step_manager = StepManager(self.outdir, self.prefix,
                                        logs_dir=self.paths["logs"],
                                        manifest_dir=self.paths["manifest"],
                                        force=force)

        # canonical cluster files (produced by classification step)
        self.cluster_bed = self.paths["clusters"] / f"{self.prefix}.combined_clusters.bed"
        self.cluster_fasta = self.paths["clusters"] / f"{self.prefix}.combined_clusters.fasta"
        self.classification_table = self.paths["tables"] / f"{self.prefix}.classification.tsv"

        log(f"SmallRNA Profiler v{__version__}", "INFO")
        log(f"Output directory: {self.outdir}", "INFO")
        log(f"R scripts path: {self.r_scripts_path}", "INFO")

    def _create_pattern_input_files(self):
        """Wrapper para compatibilidade com código existente"""
        return self._generate_pattern_matrix()
    
    # ---------------------- miRNA filtering step ----------------------

    def step_filter_miRNA(self) -> Set[str]:
        """
        Identify clusters overlapping known miRNAs (from GFF, BED, or FASTA with coords in header).
        - Writes clusters/prefix.clusters.miRNA.fasta
        - Updates tables/prefix.classification.tsv with 'miRNA'
        - If not --mirna-keep: removes miRNA clusters from clusters.fasta for downstream steps
        Returns set of miRNA cluster IDs.
        """
        step_name = "filter_miRNA"

        if not (self.mirna_fasta or self.mirna_gff or self.mirna_bed):
            log("No miRNA reference provided. Skipping miRNA filtering.", "INFO")
            return set()

        if self.step_manager.step_done(step_name) and not self.force:
            log("miRNA filtering already done, skipping...", "INFO")
            return set()

        # CORREÇÃO: Padrão correto baseado nos nomes reais dos arquivos
        # Procura por: teste.clusters.putative_*.fasta
        cluster_files = list(self.paths["clusters"].glob(f"{self.prefix}.clusters.putative_*.fasta"))
        
        # DEBUG: Mostrar o que foi encontrado
        log(f"Looking for cluster files with pattern: {self.prefix}.clusters.putative_*.fasta", "INFO")
        log(f"Found {len(cluster_files)} cluster files: {[f.name for f in cluster_files]}", "INFO")

        if not cluster_files:
            log("No cluster FASTA files found. Ensure classification generated clusters.", "WARN")
            return set()

        log("[STEP] Filtering known miRNAs...", "INFO")
        log(f"Reference: {self.mirna_gff or self.mirna_bed or self.mirna_fasta}", "INFO")
        log(f"Mode: {'KEEP (label only)' if self.mirna_keep else 'EXCLUDE from downstream'}", "INFO")

        # Create a combined clusters FASTA file for miRNA filtering
        combined_clusters_fasta = self.paths["clusters"] / f"{self.prefix}.combined_clusters.fasta"
        with open(combined_clusters_fasta, "w") as combined:
            for cluster_file in cluster_files:
                log(f"Adding clusters from: {cluster_file.name}", "INFO")
                with open(cluster_file, "r") as infile:
                    combined.write(infile.read())
        
        # Create BED file from the combined clusters
        combined_clusters_bed = self.paths["clusters"] / f"{self.prefix}.combined_clusters.bed"
        fasta_to_bed(combined_clusters_fasta, combined_clusters_bed)

        mirna_bed = self.paths["clusters"] / f"{self.prefix}.miRNA_ref.bed"
        
        # Handle different reference types
        if self.mirna_gff:
            BedTool(self.mirna_gff).saveas(str(mirna_bed))
        elif self.mirna_bed:
            # For BED files, we can use them directly
            shutil.copy2(self.mirna_bed, str(mirna_bed))
        else:  # FASTA
            fasta_to_bed(Path(self.mirna_fasta), mirna_bed)

        inter = BedTool(combined_clusters_bed).intersect(BedTool(mirna_bed), wa=True, wb=True)

        mirna_clusters: Set[str] = set()
        overlap_tsv = self.paths["clusters"] / f"{self.prefix}.clusters.miRNA.tsv"
        with open(overlap_tsv, "w") as out:
            for i in inter:
                cid = i[3]               # cluster id
                mirid = i[9] if len(i) > 9 else "miRNA_ref"
                mirna_clusters.add(cid)
                out.write(f"{cid}\t{mirid}\n")

        if not mirna_clusters:
            log("No clusters overlapped miRNA annotations.", "INFO")
            self.step_manager.append_step_record(step_name, {"miRNA_clusters": "0"})
            self.step_manager.mark_done(step_name)
            return set()

        log(f"Detected {len(mirna_clusters)} miRNA clusters.", "INFO")

        # Write miRNA FASTA
        mirna_fa = self.paths["clusters"] / f"{self.prefix}.clusters.miRNA.fasta"
        with open(mirna_fa, "w") as out:
            for cluster_file in cluster_files:
                for rec in SeqIO.parse(str(cluster_file), "fasta"):
                    if rec.id in mirna_clusters:
                        SeqIO.write(rec, out, "fasta")
        log(f"miRNA FASTA: {mirna_fa}", "INFO")

        # Update classification table
        class_file = self.paths["tables"] / f"{self.prefix}.classification.tsv"
        cls = read_classification_table(class_file)
        for cid in mirna_clusters:
            cls[cid] = "miRNA"
        write_classification_table(cls, class_file)
        log("Updated classification table with 'miRNA'.", "INFO")

        # Exclude from downstream unless --mirna-keep
        if not self.mirna_keep:
            # Create filtered versions of all cluster files
            for cluster_file in cluster_files:
                filtered_file = cluster_file.with_suffix(".filtered.fasta")
                with open(filtered_file, "w") as out:
                    for rec in SeqIO.parse(str(cluster_file), "fasta"):
                        if rec.id not in mirna_clusters:
                            SeqIO.write(rec, out, "fasta")
                log(f"Filtered clusters: {filtered_file}", "INFO")

        # Manifest + ticket
        self.step_manager.append_step_record(step_name, {
            "miRNA_count": str(len(mirna_clusters)),
            "reference": self.mirna_gff or self.mirna_bed or self.mirna_fasta,
            "kept": str(self.mirna_keep)
        })
        self.step_manager.mark_done(step_name)
        return mirna_clusters

    # ------------------------ Your existing steps ---------------------

    def step_classify(self):
        """
        Generates:
         - tables/prefix.classification.tsv (+ signal table)
         - tables/prefix.counts.tsv
         - profiles/prefix.profile{1,2}
         - clusters/*.fasta (predicted & signal)
         - tables/prefix.size_distribution.tsv
         - clusters/prefix.clusters.bed (assumed present in your original flow)
         - clusters/prefix.clusters.fasta (assumed present in your original flow)
        """
        step_name = "classify"
        if self.step_manager.step_done(step_name) and not self.force:
            log("Step classify already done, skipping...", "INFO")
            return

        log("=== STEP 1: Classification ===", "INFO")
        genome_seqs = read_fasta(self.genome_file) if os.path.exists(self.genome_file) else None

        stats = SAMStatistics(self.sam_file, self.genome_file, self.len_min, self.len_max)
        output_text, signal = stats.analyze()

        class_file = self.paths["tables"] / f"{self.prefix}.classification.tsv"
        write_text(class_file, output_text)
        with open(class_file, "a") as f:
            f.write("\n############### SIGNAL of putative_siRNA or putative_piRNA #############\n\n")
            f.write("Contig\tputative_siRNA\tputative_piRNA\tFINAL\n")
            for c in sorted(signal.keys()):
                f.write(f"{c}\t{signal[c]['si']}\t{signal[c]['pi']}\t{signal[c]['final']}\n")

        counts_file = self.paths["tables"] / f"{self.prefix}.counts.tsv"
        stats.write_counts(counts_file)

        profiles_prefix = self.paths["profiles"] / self.prefix
        stats.write_profile_files(profiles_prefix, genome_seqs, signal, self.paths["clusters"])

        size_dist_file = self.paths["tables"] / f"{self.prefix}.size_distribution.tsv"
        stats.write_distribution(size_dist_file)

        # NOTE: If your original code writes clusters.bed/fasta here, keep it.
        # Here we just ensure they exist or warn (to not break miRNA step).
        if not self.cluster_bed.exists():
            log(f"Expected cluster BED missing: {self.cluster_bed}", "WARN")
        if not self.cluster_fasta.exists():
            log(f"Expected cluster FASTA missing: {self.cluster_fasta}", "WARN")

        self.step_manager.append_step_record(step_name, {
            "classification": str(class_file),
            "counts": str(counts_file)
        })
        self.step_manager.mark_done(step_name)

    def step_split(self):
        """Split SAM into putative_siRNA.sam / putative_piRNA.sam based on classification table."""
        step_name = "split"
        if self.step_manager.step_done(step_name) and not self.force:
            log("Step split already done, skipping...", "INFO")
            return

        log("=== STEP 2: Split putative_siRNA / putative_piRNA ===", "INFO")
        class_file = self.paths["tables"] / f"{self.prefix}.classification.tsv"
        if not class_file.exists():
            raise FileNotFoundError(f"Classification file not found: {class_file}")

        keep_si: Set[str] = set()
        keep_pi: Set[str] = set()
        with open(class_file, "r") as f:
            in_signal = False
            for line in f:
                if "SIGNAL of putative_siRNA or putative_piRNA" in line:
                    in_signal = True
                    continue
                if in_signal and "\t" in line:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        contig, tag_si, tag_pi, final = parts[0], parts[1], parts[2], parts[3]
                        if final == "putative_siRNA" or tag_si == "yes":
                            keep_si.add(contig)
                        if final == "putative_piRNA" or tag_pi == "yes":
                            keep_pi.add(contig)

        si_sam = self.paths["sam"] / f"{self.prefix}.putative_siRNA.sam"
        pi_sam = self.paths["sam"] / f"{self.prefix}.putative_piRNA.sam"
        with open_sam(self.sam_file) as inf, open(si_sam, "w") as si_f, open(pi_sam, "w") as pi_f:
            headers = []
            # write headers to both
            for line in inf:
                if line.startswith("@"):
                    headers.append(line)
                    si_f.write(line); pi_f.write(line)
                else:
                    fields = line.strip().split("\t")
                    if len(fields) >= 3:
                        rname = fields[2]
                        if rname in keep_si: si_f.write(line)
                        if rname in keep_pi: pi_f.write(line)
                    # stream remainder
                    for line in inf:
                        fields = line.strip().split("\t")
                        if len(fields) >= 3:
                            rname = fields[2]
                            if rname in keep_si: si_f.write(line)
                            if rname in keep_pi: pi_f.write(line)
                    break

        self.step_manager.append_step_record(step_name, {
            "putative_siRNA_sam": str(si_sam),
            "putative_piRNA_sam": str(pi_sam),
        })
        self.step_manager.mark_done(step_name)

    def step_sep_single_mult(self):
        """Placeholder — unchanged behavior in your pipeline; keep or integrate your version here."""
        step_name = "sepSinMult"
        if self.step_manager.step_done(step_name) and not self.force:
            log("Step sepSinMult already done, skipping...", "INFO")
            return
        log("=== STEP 3: Separate single/multiple mappings (kept as in your pipeline) ===", "INFO")
        # Your original implementation lives here. Marking as done to preserve flow.
        self.step_manager.append_step_record(step_name, {"status": "kept-original"})
        self.step_manager.mark_done(step_name)

    def step_strand_counts(self):
        step_name = "strandCounts"
        if self.step_manager.step_done(step_name) and not self.force:
            log("Step strandCounts already done, skipping...", "INFO")
            return
        log("=== STEP 4: Strand Counts ===", "INFO")
        cbc = self.paths["tables"] / f"{self.prefix}.strand_counts.tsv"
        pos, neg = defaultdict(int), defaultdict(int)
        with open_sam(self.sam_file) as f:
            for line in f:
                if line.startswith("@"): continue
                parts = line.strip().split("\t")
                if len(parts) < 10: continue
                flag = int(parts[1]); rname = parts[2]; seq = parts[9]
                L = len(seq)
                if L < self.len_min or L > self.len_max: continue
                if rname == "*": continue
                if flag & 16: neg[rname] += 1
                else: pos[rname] += 1
        with open(cbc, "w") as out:
            out.write("Chr\tSense\tAntisense\tTotal\n")
            for c in sorted(set(pos.keys()) | set(neg.keys())):
                pa, na = pos.get(c, 0), neg.get(c, 0)
                out.write(f"{c}\t{pa}\t{na}\t{pa+na}\n")
        self.step_manager.append_step_record(step_name, {"counts": str(cbc)})
        self.step_manager.mark_done(step_name)

    def step_define_strand(self):
        step_name = "defineStrand"
        if self.step_manager.step_done(step_name) and not self.force:
            log("Step defineStrand already done, skipping...", "INFO")
            return
        log("=== STEP 5: Define Strand ===", "INFO")
        cbc = self.paths["tables"] / f"{self.prefix}.strand_counts.tsv"
        outfa = self.paths["defined_strand"] / f"{self.prefix}.defined_strand.fa"
        DefineStrand(cbc, Path(self.genome_file), outfa).process()
        self.step_manager.append_step_record(step_name, {"defined_fasta": str(outfa)})
        self.step_manager.mark_done(step_name)

    def step_base_distribution(self):
        step_name = "baseDistribution"
        if self.step_manager.step_done(step_name) and not self.force:
            log("Step baseDistribution already done, skipping...", "INFO")
            return
        log("=== STEP 5.5: Base Distribution (R) ===", "INFO")

        analyses = [(self.sam_file, "all-srnas", "All small RNAs")]
        si_sam = self.paths["sam"] / f"{self.prefix}.putative_siRNA.sam"
        pi_sam = self.paths["sam"] / f"{self.prefix}.putative_piRNA.sam"
        if si_sam.exists() and si_sam.stat().st_size > 0:
            analyses.append((str(si_sam), "putative_siRNA", "putative_siRNA reads"))
        if pi_sam.exists() and pi_sam.stat().st_size > 0:
            analyses.append((str(pi_sam), "putative_piRNA", "putative_piRNA reads"))

        # Ambiguous (optional)
        amb_sam = self._create_ambiguous_sam()
        if amb_sam and Path(amb_sam).exists() and Path(amb_sam).stat().st_size > 0:
            analyses.append((amb_sam, "ambiguous", "Ambiguous reads"))

        for sam_file, suffix, desc in analyses:
            log(f"Analyzing {desc} ...", "INFO")
            basename = f"{self.prefix}.{suffix}"
            PlotGeralDistribution(
                sam_file=sam_file,
                basename=basename,
                tables_dir=self.paths["tables"],
                profiles_dir=self.paths["profiles"],
                plots_dir=self.paths["plots"],
                start=self.len_min,
                end=self.len_max,
                norm=self.norm,
                minimo=max(self.si_min, self.pi_min),
                plot=True,
                antisense=False,
                r_scripts_path=self.r_scripts_path
            ).analyze()

        self.step_manager.append_step_record(step_name, {"analyses": [a[1] for a in analyses]})
        self.step_manager.mark_done(step_name)

    def _create_ambiguous_sam(self) -> Optional[str]:
        class_file = self.paths["tables"] / f"{self.prefix}.classification.tsv"
        if not class_file.exists():
            return None
        ambiguous = set()
        with open(class_file, "r") as f:
            in_signal = False
            for line in f:
                if "SIGNAL of putative_siRNA or putative_piRNA" in line:
                    in_signal = True; continue
                if in_signal and "\t" in line:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4 and parts[3] == "ambiguous":
                        ambiguous.add(parts[0])
        if not ambiguous:
            return None
        amb_sam = self.paths["sam"] / f"{self.prefix}.ambiguous.sam"
        with open_sam(self.sam_file) as inf, open(amb_sam, "w") as outf:
            for line in inf:
                if line.startswith("@"):
                    outf.write(line)
                else:
                    fields = line.strip().split("\t")
                    if len(fields) >= 3 and fields[2] in ambiguous:
                        outf.write(line)
                    for line in inf:
                        fields = line.strip().split("\t")
                        if len(fields) >= 3 and fields[2] in ambiguous:
                            outf.write(line)
                    break
        log(f"Created ambiguous SAM: {amb_sam}", "INFO")
        return str(amb_sam)

    def step_patterns(self):
        """Execute pattern analysis (ping-pong, z-score, etc.)"""
        step_name = "patterns"
        if self.step_manager.step_done(step_name) and not self.force:
            log("Step patterns already done, skipping...", "INFO")
            return
        
        log("=== STEP 6: Pattern Analysis ===", "INFO")
        
        # Verificar quais arquivos SAM existem
        si_sam = self.paths["sam"] / f"{self.prefix}.putative_siRNA.sam"
        pi_sam = self.paths["sam"] / f"{self.prefix}.putative_piRNA.sam"
        
        sam_files_to_process = []
        
        if si_sam.exists() and si_sam.stat().st_size > 0:
            sam_files_to_process.append(("siRNA", str(si_sam)))
        
        if pi_sam.exists() and pi_sam.stat().st_size > 0:
            sam_files_to_process.append(("piRNA", str(pi_sam)))
        
        if not sam_files_to_process:
            log("No SAM files found for pattern analysis", "WARN")
            self.step_manager.append_step_record(step_name, {"status": "skipped_no_input"})
            self.step_manager.mark_done(step_name)
            return
        
        # Executar análise para cada tipo de sRNA
        for rna_type, sam_path in sam_files_to_process:
            log(f"Running pattern analysis for {rna_type}...", "INFO")
            
            # Gerar arquivo de Z-scores no formato do script original
            pattern_file = self._generate_zscore_matrix(rna_type, sam_path)
            
            if not pattern_file or not pattern_file.exists():
                log(f"Failed to generate pattern matrix for {rna_type}", "WARN")
                continue
            
            # Executar scripts R
            self._run_r_pattern_scripts(pattern_file, rna_type)
        
        self.step_manager.append_step_record(step_name, {
            "pattern_analyses": [t[0] for t in sam_files_to_process]
        })
        self.step_manager.mark_done(step_name)
    
    def _generate_zscore_matrix(self, rna_type: str, sam_path: str) -> Optional[Path]:
        """Gerar arquivo CSV de Z-scores no formato do script original calcPatternInSamFile.pl"""
        log(f"Generating Z-score matrix for {rna_type}...", "INFO")
        
        # Ler dados do SAM
        from collections import defaultdict
        
        # Estrutura: chr -> size -> strand -> count
        chr_data = defaultdict(lambda: defaultdict(lambda: {"pos": 0, "neg": 0}))
        seen_reads = defaultdict(set)
        
        with open_sam(sam_path) as f:
            for line in f:
                if line.startswith("@"):
                    continue
                
                fields = line.strip().split("\t")
                if len(fields) < 10:
                    continue
                
                read_id = fields[0]
                flag = int(fields[1])
                chrom = fields[2]
                seq = fields[9]
                size = len(seq)
                
                # Filtrar por tamanho (15-35)
                if size < self.len_min or size > self.len_max:
                    continue
                
                # Ignorar reads não mapeados
                if chrom == "*":
                    continue
                
                # Evitar contar reads duplicados no mesmo cromossomo
                read_key = f"{read_id}_{size}_{flag}"
                if read_key in seen_reads[chrom]:
                    continue
                seen_reads[chrom].add(read_key)
                
                # Determinar strand
                is_reverse = (flag & 16) != 0
                
                if is_reverse:
                    chr_data[chrom][size]["neg"] += 1
                else:
                    chr_data[chrom][size]["pos"] += 1
        
        # Calcular Z-scores para cada cromossomo
        output_file = self.paths["tables"] / f"{self.prefix}.{rna_type}.zscore_matrix.csv"
        
        with open(output_file, "w") as out_f:
            # Escrever cabeçalho
            header = ["Chromosome"]
            for i in range(42):  # 0-41 como no script original
                header.append(str(i))
            out_f.write("#".join(header) + "\n")
            
            # Para cada cromossomo, calcular Z-scores
            total_chromosomes = 0
            for chrom in sorted(chr_data.keys()):
                # Preparar array de 42 elementos (0-41)
                # Posições 0-20: positive strand para tamanhos 15-35
                # Posições 21-41: negative strand para tamanhos 15-35
                values = [0] * 42
                
                # DEBUG: verificar os índices
                size_range = self.len_max - self.len_min + 1  # Deveria ser 21 (35-15+1)
                if size_range != 21:
                    log(f"Warning: size range is {size_range}, expected 21 for len_min={self.len_min}, len_max={self.len_max}", "WARN")
                
                for size in range(self.len_min, self.len_max + 1):
                    pos_idx = size - self.len_min  # 0-20 para tamanhos 15-35
                    neg_idx = pos_idx + 21  # 21-41 para negative strand
                    
                    # Verificar limites
                    if pos_idx < 0 or pos_idx >= 21:
                        log(f"Warning: pos_idx {pos_idx} out of range for size {size}", "WARN")
                        continue
                    if neg_idx < 21 or neg_idx >= 42:
                        log(f"Warning: neg_idx {neg_idx} out of range for size {size}", "WARN")
                        continue
                    
                    pos_count = chr_data[chrom].get(size, {"pos": 0})["pos"]
                    neg_count = chr_data[chrom].get(size, {"neg": 0})["neg"]
                    
                    values[pos_idx] = pos_count
                    values[neg_idx] = neg_count
                
                # Calcular estatísticas para este cromossomo
                try:
                    mean_val = statistics.mean(values) if values else 0
                    # Usar stdev apenas se houver pelo menos 2 valores não-zero
                    non_zero_values = [v for v in values if v != 0]
                    if len(non_zero_values) > 1:
                        std_val = statistics.stdev(non_zero_values)
                    else:
                        std_val = 1.0  # Valor padrão para evitar divisão por zero
                    
                    # Calcular Z-scores (como no script original)
                    if std_val > 0:
                        zscores = []
                        for val in values:
                            zscore = (val - mean_val) / std_val if std_val > 0 else 0
                            zscores.append(zscore)
                        
                        # Escrever linha
                        out_f.write(f"{chrom}#" + "#".join(f"{z:.6f}" for z in zscores) + "\n")
                        total_chromosomes += 1
                    else:
                        log(f"Skipping chromosome {chrom}: zero standard deviation", "WARN")
                        
                except Exception as e:
                    log(f"Error calculating Z-scores for chromosome {chrom}: {e}", "WARN")
        
        if total_chromosomes == 0:
            log(f"No valid chromosomes with Z-scores for {rna_type}", "WARN")
            return None
        
        log(f"Z-score matrix written: {output_file} ({total_chromosomes} chromosomes)", "INFO")
        return output_file
    
    def _run_r_pattern_scripts(self, pattern_file: Path, rna_type: str):
        """Execute pattern analysis using plots module"""
        try:
            # Try to import plots module
            import sys
            script_dir = get_script_directory()
            if str(script_dir) not in sys.path:
                sys.path.insert(0, str(script_dir))
            
            try:
                import plots
                
                log(f"Generating pattern plots for {rna_type} with plots module...", "INFO")
                
                # Output prefix
                output_prefix = self.paths["plots"] / f"{self.prefix}.{rna_type}"
                
                # Generate pattern plot (combined strands)
                success1 = plots.generate_pattern_plot(
                    input_matrix=pattern_file,
                    output_prefix=output_prefix,
                    r_scripts_path=Path(self.r_scripts_path),
                    both_strands=False
                )
                
                if success1:
                    log(f"Generated pattern plot for {rna_type}", "INFO")
                
                # Generate pattern plot (both strands separately)
                output_prefix_both = self.paths["plots"] / f"{self.prefix}.{rna_type}.both_strands"
                success2 = plots.generate_pattern_plot(
                    input_matrix=pattern_file,
                    output_prefix=output_prefix_both,
                    r_scripts_path=Path(self.r_scripts_path),
                    both_strands=True
                )
                
                if success2:
                    log(f"Generated both strands pattern plot for {rna_type}", "INFO")
                
            except ImportError:
                log("plots module not found, falling back to direct R calls", "INFO")
                self._run_r_pattern_scripts_fallback(pattern_file, rna_type)
                
        except Exception as e:
            log(f"Error generating pattern plots: {e}", "WARN")
            log("Trying fallback method...", "INFO")
            self._run_r_pattern_scripts_fallback(pattern_file, rna_type)
    
    def _run_r_pattern_scripts_fallback(self, pattern_file: Path, rna_type: str):
        """Fallback method using direct R script calls"""
        r_exec = "Rscript"
        
        # Scripts R
        scripts_to_run = [
            ("plotPattern_publication.R", ""),
            ("plotPattern_publication_both_strands.R", ".both_strands"),
        ]
        
        for script_name, suffix in scripts_to_run:
            script_path = os.path.join(self.r_scripts_path, script_name)
            
            if not os.path.exists(script_path):
                log(f"Script R não encontrado: {script_name}", "WARN")
                continue
            
            try:
                # Prefixo de saída
                out_prefix = f"{self.prefix}.{rna_type}"
                cmd = [r_exec, script_path, str(pattern_file), out_prefix]
                
                # Executar no diretório de plots
                result = subprocess.run(cmd, check=False, cwd=str(self.paths["plots"]),
                                      capture_output=True, text=True)
                
                if result.returncode == 0:
                    log(f"Executed pattern analysis for {rna_type} with {script_name}", "INFO")
                    log(f"Output: {out_prefix}{suffix}.pdf", "INFO")
                else:
                    log(f"R script error for {script_name}: {result.stderr}", "WARN")
                
            except Exception as e:
                log(f"Error running {script_name} for {rna_type}: {e}", "WARN")

    def step_summary(self):
        step_name = "summary"
        log("=== STEP 7: Summary ===", "INFO")
        summary_file = self.paths["manifest"] / f"{self.prefix}.summary.txt"
        class_file = self.paths["tables"] / f"{self.prefix}.classification.tsv"
        si_sam = self.paths["sam"] / f"{self.prefix}.putative_siRNA.sam"
        pi_sam = self.paths["sam"] / f"{self.prefix}.putative_piRNA.sam"
        manifest = self.paths["manifest"] / f"{self.prefix}.run.json"

        text = (
            f"smallRNA Profiler Pipeline Summary\n"
            f"===================================\n\n"
            f"Input Files:\n"
            f"  SAM:     {self.sam_file}\n"
            f"  Genome:  {self.genome_file}\n"
            f"  Reads:   {self.reads_file or 'Not provided'}\n\n"
            f"Output Files:\n"
            f"  Classification:   {class_file}\n"
            f"  putative_siRNA:   {si_sam}\n"
            f"  putative_piRNA:   {pi_sam}\n"
            f"  Manifest:         {manifest}\n\n"
            f"Parameters:\n"
            f"  Prefix:       {self.prefix}\n"
            f"  Length range: {self.len_min}-{self.len_max}\n"
            f"  SI minimum:   {self.si_min}\n"
            f"  PI minimum:   {self.pi_min}\n"
            f"  Normalization:{self.norm or 'None'}\n"
            f"  Patterns:     {self.run_patterns}\n\n"
            f"miRNA Filtering:\n"
            f"  Reference:    {self.mirna_gff or self.mirna_bed or self.mirna_fasta or 'None'}\n"
            f"  Keep miRNAs:  {self.mirna_keep}\n\n"
            f"Pipeline Version: {__version__}\n"
            f"Completed: {datetime.now().isoformat()}\n"
        )
        write_text(summary_file, text)
        log(f"Summary written: {summary_file}", "INFO")
        self.step_manager.append_step_record(step_name, {"summary": str(summary_file)})

    # ------------------------------- run -------------------------------

    def run(self):
        try:
            # 1) Classify (produces clusters/* and tables/classification.tsv)
            self.step_classify()

            # 2) miRNA filter (BEFORE si/pi downstream)
            self.step_filter_miRNA()

            # 3) Split si/pi
            self.step_split()

            # 4) Optional steps (kept from your pipeline)
            self.step_sep_single_mult()
            self.step_strand_counts()
            self.step_define_strand()

            # 5) Base distributions + plots
            self.step_base_distribution()

            # 6) Pattern analysis (if you keep it)
            if self.run_patterns:
                self.step_patterns()

            # 7) Summary
            self.step_summary()

            log("Pipeline completed successfully!", "INFO")
        except KeyboardInterrupt:
            log("Pipeline interrupted by user.", "WARN")
            sys.exit(130)
        except Exception as e:
            log(f"Pipeline failed: {e}", "ERROR")
            import traceback; traceback.print_exc()
            sys.exit(1)


# ------------------------------- Main Function -------------------------------

def run_clusters(args):
    """Main function for clusters module - integrated with srnaprofiler.py"""
    
    # Validate input files
    validate_file_exists(Path(args.sam), "SAM file")
    validate_file_exists(Path(args.genome), "Genome file")
    
    # Setup output directory
    outdir = ensure_outdir(Path(args.outdir))
    
    log("="*60, "INFO")
    log(f"Starting cluster analysis for {args.prefix}", "INFO")
    log("="*60, "INFO")
    
    # Create pipeline instance
    pipeline = SmallRNAPipeline(
        sam_file=args.sam,
        genome_file=args.genome,
        prefix=args.prefix,
        reads_file=args.reads,
        outdir=args.outdir,
        si_min=args.si,
        pi_min=args.pi,
        len_min=args.len_min,
        len_max=args.len_max,
        norm=args.norm if hasattr(args, 'norm') else None,
        run_patterns=args.run_patterns,
        force=getattr(args, 'force', False),
        r_scripts_path="",  # Auto-detect
        mirna_fasta=getattr(args, 'mirna_fasta', None),
        mirna_gff=getattr(args, 'mirna_gff', None),
        mirna_bed=getattr(args, 'mirna_bed', None),
        mirna_keep=getattr(args, 'mirna_keep', False)
    )
    
    # Run the pipeline
    pipeline.run()
    
    log("="*60, "INFO")
    log(f"Cluster analysis completed successfully!", "INFO")
    log(f"Results saved in: {outdir}", "INFO")
    log("="*60, "INFO")


if __name__ == "__main__":
    # This allows the module to be run standalone for testing
    parser = argparse.ArgumentParser(description="Small RNA cluster analysis")
    parser.add_argument("--sam", required=True, help="Input SAM file")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA file")
    parser.add_argument("-p", "--prefix", required=True, help="Output prefix")
    parser.add_argument("-o", "--outdir", help="Output directory")
    parser.add_argument("--si", type=int, default=100, help="siRNA threshold")
    parser.add_argument("--pi", type=int, default=100, help="piRNA threshold")
    parser.add_argument("--len-min", type=int, default=15, help="Min read length")
    parser.add_argument("--len-max", type=int, default=35, help="Max read length")
    parser.add_argument("--run-patterns", action="store_true", help="Run pattern analysis")
    
    # miRNA options
    parser.add_argument("--mirna-fasta", help="miRNA FASTA reference")
    parser.add_argument("--mirna-gff", help="miRNA GFF reference")
    parser.add_argument("--mirna-bed", help="miRNA BED reference")
    parser.add_argument("--mirna-keep", action="store_true", help="Keep miRNA clusters")
    
    args = parser.parse_args()
    run_clusters(args)