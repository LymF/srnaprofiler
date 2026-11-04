#!/usr/bin/env bash
# Download one single-end run from NCBI SRA and write <OUTDIR>/<ACC>.fastq.gz.
#
# Usage: download_reads_sra.sh ACCESSION OUTDIR [THREADS]
#
# prefetch is md5-verified and resumable; vdb-validate checks the .sra before conversion;
# the fastq.gz appears atomically (<ACC>.fastq.gz.part, then mv), so an existing
# <ACC>.fastq.gz is always complete and is skipped. Needs sra-tools and pigz
# (see benchmark/environment.yml). Paired-end runs (more than one fastq) are refused.
set -u
if [ $# -lt 2 ]; then
  echo "usage: $0 ACCESSION OUTDIR [THREADS]" >&2
  exit 2
fi
acc=$1
R=$2
threads=${3:-6}
mkdir -p "$R"
[ -s "$R/$acc.fastq.gz" ] && { echo "skip $acc"; exit 0; }
T=$R/.sra_$acc
mkdir -p "$T"
cd "$T" || exit 1
for i in 1 2 3 4 5; do
  prefetch "$acc" -O "$T" --max-size 100g >> prefetch.log 2>&1 && break
  echo "retry $i" >> prefetch.log
done
sra=$(ls "$T/$acc/$acc.sra" 2>/dev/null) || { echo "FAIL $acc: prefetch"; exit 1; }
vdb-validate "$sra" >> validate.log 2>&1 || { echo "FAIL $acc: vdb-validate"; exit 1; }
fasterq-dump "$sra" -O "$T" -e "$threads" --split-files -t "$T/tmp" >> dump.log 2>&1 \
  || { echo "FAIL $acc: fasterq-dump"; exit 1; }
n=$(ls "$T"/*.fastq | wc -l)
[ "$n" = 1 ] || { echo "FAIL $acc: $n fastq files (paired?)"; exit 1; }
pigz -p "$threads" -c "$T"/*.fastq > "$R/$acc.fastq.gz.part" \
  && mv "$R/$acc.fastq.gz.part" "$R/$acc.fastq.gz" \
  && rm -rf "$T" \
  && echo "ok $acc $(du -h "$R/$acc.fastq.gz" | cut -f1)"
