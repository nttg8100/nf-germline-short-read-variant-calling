#!/bin/bash

# Benchmark script for HG002 with DeepVariant, FreeBayes, HaplotypeCaller, Manta, Strelka (relative paths)
set -euo pipefail

HGREF="references/Homo_sapiens_assembly38.fasta"
REF="HG002"
BED="references/${REF}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
TRUTH="references/${REF}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
SDF="references/grch38.sdf"
THREADS=16

RESULTSDIR="results"
mkdir -p "$RESULTSDIR"

TOOLS=(deepvariant freebayes haplotypecaller manta strelka)

VCF_DEEPVARIANT="pipelines/sarek/variant_calling/deepvariant/HG002/HG002.deepvariant.vcf.gz"
VCF_FREEBAYES="pipelines/sarek/variant_calling/freebayes/HG002/HG002.freebayes.vcf.gz"
VCF_HAPLOTYPECALLER="pipelines/sarek/variant_calling/haplotypecaller/HG002/HG002.haplotypecaller.vcf.gz"
VCF_MANTA="pipelines/sarek/variant_calling/manta/HG002/HG002.manta.diploid_sv.vcf.gz"
VCF_STRELKA="pipelines/sarek/variant_calling/strelka/HG002/HG002.strelka.variants.vcf.gz"

# Path map
declare -A VCF_PATHS=(
  [deepvariant]="$VCF_DEEPVARIANT"
  [freebayes]="$VCF_FREEBAYES"
  [haplotypecaller]="$VCF_HAPLOTYPECALLER"
  [manta]="$VCF_MANTA"
  [strelka]="$VCF_STRELKA"
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