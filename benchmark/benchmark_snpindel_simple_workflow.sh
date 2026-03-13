#!/bin/bash

# Benchmark script for HG002 with DeepVariant, FreeBayes, HaplotypeCaller, Manta, Strelka (relative paths)
set -euo pipefail

export HGREF="references/Homo_sapiens_assembly38.fasta"
REF="HG002"
BED="references/${REF}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
TRUTH="references/${REF}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
SDF="references/grch38.sdf"
THREADS=16

RESULTSDIR="results/nf-germline-short-read-variant-calling"
mkdir -p "$RESULTSDIR"

TOOLS=(freebayes deepvariant)
VCF_FREEBAYES="pipelines/nf-germline-short-read-variant-calling/variant_calling/freebayes/HG002.vcf.gz"
VCF_DEEPVARIANT="pipelines/nf-germline-short-read-variant-calling/variant_calling/deepvariant/HG002.vcf.gz"

# Path map
declare -A VCF_PATHS=(
  [freebayes]="$VCF_FREEBAYES"
  [deepvariant]="$VCF_DEEPVARIANT"
)

for TOOL in "${TOOLS[@]}"; do
  TOOLDIR="$RESULTSDIR/$TOOL"
  mkdir -p "$TOOLDIR"
  VCF="${VCF_PATHS[$TOOL]}"
  LOG="$TOOLDIR/${REF}_${TOOL}.log"
  SUMMARY="$TOOLDIR/${REF}_${TOOL}.summary.csv"
  FORMATTED="$TOOLDIR/${REF}_${TOOL}.formatted.summary.csv"
  OUTPUT="$TOOLDIR/${REF}_${TOOL}"

  # If the tool output folder already exists with summary file, continue from AWK formatting (line 55)
  if [[ -f "${OUTPUT}.summary.csv" ]]; then
    echo "Output for $TOOL already exists, formatting summary."
    awk -v d="${TOOL}_${REF}" -F',' 'FNR==1{a="run"} FNR>1{a=d} {print $0",\t"a}' "${OUTPUT}.summary.csv" > "$FORMATTED"
    echo "Completed $TOOL (using existing output)"
    continue
  fi

  if [[ ! -f "$VCF" ]]; then
    echo "VCF file not found for $TOOL: $VCF" >&2
    continue
  fi

  echo "Benchmarking $TOOL for $REF ..."
  hap.py "$TRUTH" "$VCF" \
    -o "$OUTPUT" \
    -V --engine=vcfeval --threads "$THREADS" --engine-vcfeval-template "$SDF" \
    -f "$BED" \
    --logfile "$LOG" \
    --scratch-prefix .

  if [[ -f "${OUTPUT}.summary.csv" ]]; then
    awk -v d="${TOOL}_${REF}" -F',' 'FNR==1{a="run"} FNR>1{a=d} {print $0",\t"a}' "${OUTPUT}.summary.csv" > "$FORMATTED"
  else
    echo "Warning: ${OUTPUT}.summary.csv not found for $TOOL" >&2
  fi

  echo "Completed $TOOL"
done

MERGED="$RESULTSDIR/merged_${REF}_benchmark.csv"
echo "Merging summaries to $MERGED ..."
awk '(NR == 1) || (FNR > 1)' "$RESULTSDIR"/*/*.formatted.summary.csv > "$MERGED"
echo "All done. Results for $REF in $MERGED"