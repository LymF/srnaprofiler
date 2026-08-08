#!/usr/bin/env bash
# Build a bowtie index of structural RNA, to filter those reads before clustering.
#
# tRNA, rRNA, snoRNA and snRNA fragments are abundant in any small RNA library --
# 36.5% of the read mass in our Arabidopsis run, 25.7% in the Drosophila one --
# and none of them is siRNA or piRNA. Removing the reads is better than rejecting
# the loci afterwards: the locus never exists, so its abundance is never wrong.
#
# Sequences come from Rfam, which covers every species. That is the point: a
# non-model genome has no structural RNA annotation to extract from, while tRNA
# and rRNA are conserved enough for a 20-30 nt read to match across species at
# one mismatch. There is no single FASTA of all of Rfam -- it ships 4227
# per-family files -- and we do not want all of them anyway.
#
# Rfam accessions are permanent -- RF00005 is tRNA in every release -- and the
# path uses CURRENT, so no version is pinned here. What can still break is the
# EBI path itself; the script fails loudly rather than building a partial index.
#
# Usage:  bash extras/fetch_ncrna_db.sh <outdir>
# Then:   srnaprofiler mapping ... --ncrna-db <outdir>/ncrna

set -euo pipefail

OUT="${1:?usage: fetch_ncrna_db.sh <outdir>}"
BASE="https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files"

# The families that dominate small RNA libraries. tRNA and rRNA first: they are
# both the most abundant contaminant and the most conserved, so they transfer to
# a distant organism. The snRNA entries are less conserved across clades -- keep
# them, but do not expect the same reach in a non-model species.
FAMILIES=(
  RF00005  # tRNA
  RF00001  # 5S ribosomal RNA
  RF00002  # 5.8S ribosomal RNA
  RF01960  # eukaryotic small subunit rRNA
  RF02543  # eukaryotic large subunit rRNA
  RF00004  # U2 spliceosomal RNA
  RF00003  # U1 spliceosomal RNA
)

command -v bowtie-build >/dev/null || { echo "bowtie-build not in PATH" >&2; exit 1; }

mkdir -p "$OUT"
: > "$OUT/ncrna.fa"

# Every family must arrive with sequences in it. A truncated download would
# otherwise build a thinner index in silence, which is the same silent false
# negative this filter exists to avoid.
for rf in "${FAMILIES[@]}"; do
  echo "fetching $rf" >&2
  curl -sSL --fail -m 600 "$BASE/$rf.fa.gz" | gunzip -c > "$OUT/.$rf.fa"
  n=$(grep -c '^>' "$OUT/.$rf.fa" || true)
  if [ "$n" -lt 100 ]; then
    echo "ERROR: $rf returned $n sequences, expected hundreds or more." >&2
    echo "  Either the download was truncated, or Rfam moved the path." >&2
    echo "  Check $BASE/$rf.fa.gz by hand before trusting this index." >&2
    rm -f "$OUT/.$rf.fa"
    exit 1
  fi
  echo "  $rf: $n sequences" >&2
  cat "$OUT/.$rf.fa" >> "$OUT/ncrna.fa"
  rm -f "$OUT/.$rf.fa"
done

echo "building bowtie index" >&2
bowtie-build --threads "${THREADS:-4}" "$OUT/ncrna.fa" "$OUT/ncrna" >/dev/null

n=$(grep -c '^>' "$OUT/ncrna.fa")
echo "done: $n sequences in $OUT/ncrna.fa, index at $OUT/ncrna" >&2
