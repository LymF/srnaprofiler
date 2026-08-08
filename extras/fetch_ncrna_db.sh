#!/usr/bin/env bash
# Build a bowtie index of structural RNA (rRNA, tRNA, snRNA, snoRNA) from Rfam,
# to filter those reads before clustering. miRNA families are excluded: the tool
# calls miRNA as a class.
#
# Usage:  bash extras/fetch_ncrna_db.sh <outdir>
# Then:   srnaprofiler mapping ... --ncrna-db <outdir>/ncrna

set -euo pipefail

OUT="${1:?usage: fetch_ncrna_db.sh <outdir>}"
BASE="https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT"

command -v bowtie-build >/dev/null || { echo "bowtie-build not in PATH" >&2; exit 1; }

mkdir -p "$OUT"

# Select on Rfam's own type field, so a family added later is picked up. Matching
# the whole line instead would take in whatever merely mentions the words in its
# description -- RNase P, tmRNA, telomerase, group I/II introns, viral tRNA-like
# 3' UTR elements (63 families) -- and would drop the five Giardia snoRNAs whose
# names contain "miR".
echo "selecting families" >&2
curl -sSL --fail -m 300 "$BASE/database_files/family.txt.gz" | gunzip -c \
  | LC_ALL=C awk -F'\t' '
      { type = "" ; for (i = 1; i <= NF; i++) if ($i ~ /^Gene;/) { type = $i; break } }
      type ~ /rRNA|tRNA|snRNA|snoRNA/ && $1 ~ /^RF[0-9]+$/ { print $1 }
    ' | sort -u > "$OUT/families.txt"

n_fam=$(wc -l < "$OUT/families.txt")
if [ "$n_fam" -lt 200 ]; then
  echo "ERROR: $n_fam families selected, expected several hundred." >&2
  echo "  Check $BASE/database_files/family.txt.gz" >&2
  exit 1
fi
echo "  $n_fam families" >&2

: > "$OUT/ncrna.fa"
i=0
while read -r rf; do
  i=$((i + 1))
  printf '\r  fetching %d/%d (%s)    ' "$i" "$n_fam" "$rf" >&2
  if ! curl -sSL --fail -m 300 "$BASE/fasta_files/$rf.fa.gz" 2>/dev/null | gunzip -c >> "$OUT/ncrna.fa"; then
    printf '\r  %s: no FASTA, skipped                    \n' "$rf" >&2
  fi
done < "$OUT/families.txt"
printf '\n' >&2

# A truncated run would build a thin index without saying so.
n_seq=$(grep -c '^>' "$OUT/ncrna.fa")
if [ "$n_seq" -lt 100000 ]; then
  echo "ERROR: $n_seq sequences, expected millions. Downloads were truncated." >&2
  exit 1
fi

echo "building index over $n_seq sequences" >&2
bowtie-build --threads "${THREADS:-4}" "$OUT/ncrna.fa" "$OUT/ncrna" >/dev/null
echo "done: index at $OUT/ncrna" >&2
