#!/usr/bin/env bash
# Re-run classify on every manifest library, reusing qc/mapping/mirna outputs.
# benchmark run skips steps whose outputs are newer than their inputs (mtimes only),
# so a pipeline code change needs the clusters table deleted to force classify.
# Usage: rerun_classify.sh <repo> <results_dir> <status_file> [parallel=14] [threads=4]
# Requires PYTHON (the srnaprofiler env python) and BENCHMARK_DATA (the shared data
# dir with genomes/ and reads/); without BENCHMARK_DATA the checkout's own
# data/benchmark is used and genomes/reads are re-downloaded.
set -u
REPO=$1; RES=$2; STATUS=$3; PAR=${4:-14}; THR=${5:-4}
PY=${PYTHON:?set PYTHON to the srnaprofiler env python}
: "${BENCHMARK_DATA:?set BENCHMARK_DATA to the shared benchmark data dir}"
export BENCHMARK_DATA SRNAPROFILER_PYTHON=$PY PATH="$(dirname "$PY"):$PATH"
cd "$REPO" || exit 1
: > "$STATUS"
tail -n +2 benchmark/manifest.tsv | cut -f1 | while read -r acc; do
  rm -f "$RES/$acc/work/$acc.clusters.tsv"
  echo "$acc"
done | xargs -P "$PAR" -I{} bash -c \
  "$PY -m benchmark run --library {} --results-dir $RES --threads $THR > $RES/../rerun_{}.log 2>&1; echo {} exit=\$? >> $STATUS"
echo ALL_DONE >> "$STATUS"
