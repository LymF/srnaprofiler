#!/bin/bash

# Usage: ./run_rnafold.sh input.fasta output_dir

FASTA=$1
OUTDIR=$2

mkdir -p "$OUTDIR"

# Run RNAfold sequence by sequence, to keep header names
awk '/^>/ {if (seq) {print header"\n"seq; seq=""} header=$0; next} {seq=seq$0} END {print header"\n"seq}' "$FASTA" | \
while read -r header; do
    read -r seq
    name=$(echo "$header" | sed 's/^>//; s/[^A-Za-z0-9_.-]/_/g')  # clean name

    echo "$header" >> "$OUTDIR/rnafold_results.txt"
    echo "$seq"    >> "$OUTDIR/rnafold_results.txt"

    echo "$seq" | RNAfold -p --MEA --noPS > "$OUTDIR/${name}_results.txt"

    # RNAfold generates rna.ps and dot.ps → rename them
    if [[ -f rna.ps ]]; then
        mv rna.ps "$OUTDIR/${name}_ss.ps"
        ps2pdf "$OUTDIR/${name}_ss.ps" "$OUTDIR/${name}_ss.pdf"
    fi
    if [[ -f dot.ps ]]; then
        mv dot.ps "$OUTDIR/${name}_dp.ps"
        ps2pdf "$OUTDIR/${name}_dp.ps" "$OUTDIR/${name}_dp.pdf"
    fi

    # append results
    cat "$OUTDIR/${name}_results.txt" >> "$OUTDIR/rnafold_results.txt"
    echo "" >> "$OUTDIR/rnafold_results.txt"
done

echo "✅ RNAfold analysis complete."
echo "Results saved in: $OUTDIR"
