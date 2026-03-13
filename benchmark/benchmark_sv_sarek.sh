#!/bin/bash

# Benchmark script for HG002 SV variants using Truvari (Sarek pipeline)
# Converts query VCFs from hg38 to hg19 to match truth set coordinates
set -euo pipefail

HGREF="references/human_g1k_v37_decoy.fasta"
REF="HG002"

# Truth set in hg19 (original coordinates)
SV_TRUTH="references/HG002_SVs_Tier1_v0.6.vcf.gz"
SV_BED="references/HG002_SVs_Tier1_v0.6.bed"

# Chain file for hg38 to hg19 conversion
CHAIN_FILE="references/hg38ToHg19.over.chain.gz"
RESULTSDIR="results/sarek/sv"
mkdir -p "$RESULTSDIR"

# Define SV callers available from Sarek pipeline
TOOLS=(manta)

# VCF paths for SV callers from Sarek (hg38 coordinates)
VCF_MANTA="pipelines/sarek/variant_calling/manta/HG002/HG002.manta.diploid_sv.vcf.gz"

# Path map
declare -A VCF_PATHS=(
  [manta]="$VCF_MANTA"
)

for TOOL in "${TOOLS[@]}"; do
  TOOLDIR="$RESULTSDIR/$TOOL"
  mkdir -p "$TOOLDIR"
  VCF_HG38="${VCF_PATHS[$TOOL]}"
  VCF_HG19="$TOOLDIR/${REF}_${TOOL}_hg19.vcf.gz"
  LOG="$TOOLDIR/${REF}_${TOOL}.log"
  OUTPUT="$TOOLDIR/${REF}_${TOOL}"
  TRUVARI_OUTPUT="$OUTPUT.truvari"
  SUMMARY="$TOOLDIR/${REF}_${TOOL}.summary.txt"
  CONVERSION_LOG="$TOOLDIR/${REF}_${TOOL}_conversion.log"

  ##### Convert VCF from hg38 to hg19 using CrossMap
  echo "Converting $TOOL VCF from hg38 to hg19..."
  # Normalize the VCF before conversion to ensure consistent representation of variants
  bcftools norm -m -any -o "$TOOLDIR/${REF}_${TOOL}_hg38_normalized.vcf" "$VCF_HG38"

  CrossMap vcf \
    "$CHAIN_FILE" \
    "$TOOLDIR/${REF}_${TOOL}_hg38_normalized.vcf" \
    "$HGREF" \
    "$TOOLDIR/${REF}_${TOOL}_hg19_raw.vcf" \
    > "$CONVERSION_LOG" 2>&1 || {
    echo "Warning: CrossMap conversion failed for $TOOL. Check log: $CONVERSION_LOG" >&2
    continue
  }

  ##### Remove chr, sort and index the converted VCF
  # Remove chr prefix from the VCF (for consistency with truth set)
  echo "Processing chromosomes in $TOOL VCF..."
  # Remove 'chr' prefix if present, otherwise add it (to ensure consistency with truth set)
  sed -e 's/chr//g' "$TOOLDIR/${REF}_${TOOL}_hg19_raw.vcf" > "$TOOLDIR/${REF}_${TOOL}_hg19_nochr.vcf"
  # Sort the VCF
  echo "Sorting $TOOL VCF..." 
  bcftools sort -o $TOOLDIR/${REF}_${TOOL}_hg19.vcf "$TOOLDIR/${REF}_${TOOL}_hg19_nochr.vcf"


  # Index the final VCF
  bgzip -f $TOOLDIR/${REF}_${TOOL}_hg19.vcf
  tabix -p vcf "$VCF_HG19"

  # Clean up intermediate files
  rm -f "$TOOLDIR/${REF}_${TOOL}_hg19_raw.vcf" "$TOOLDIR/${REF}_${TOOL}_hg19_nochr.vcf"

  #### Benchmarking with Truvari
  # Run Truvari benchmark using pixi environment
  # Note: Truvari may need reference validation, but we use -f to provide the reference
  echo "Benchmarking $TOOL for $REF (SV variants) ..."
  truvari bench \
    -b "$SV_TRUTH" \
    -c "$VCF_HG19" \
    -o "$TRUVARI_OUTPUT" \
    -f "$HGREF" \
    --includebed "$SV_BED" \
    > "$LOG" 2>&1 || {
    echo "Warning: Truvari benchmarking failed for $TOOL. Check log: $LOG" >&2
    # Don't exit on failure, just continue with next tool
  }

  #### Extract summary statistics
  if [ -f "$TRUVARI_OUTPUT/summary.json" ]; then
    echo "Extracting summary from $TOOL Truvari output..."
    # Convert JSON summary to readable format
    python3 << PYTHON_EOF
import json
import sys

tool = "$TOOL"
ref = "$REF"
summary_file = "$TRUVARI_OUTPUT/summary.json"
output_file = "$SUMMARY"

try:
    with open(summary_file) as f:
        summary = json.load(f)
    
    with open(output_file, 'w') as f:
        f.write(f"SV Benchmarking Summary: {ref} - {tool}\\n")
        f.write("=" * 60 + "\\n\\n")
        f.write(f"TP (True Positives):       {summary['TP-base']}\\n")
        f.write(f"FP (False Positives):      {summary['FP']}\\n")
        f.write(f"FN (False Negatives):      {summary['FN']}\\n")
        f.write(f"Precision:                 {summary['precision']:.4f}\\n")
        f.write(f"Recall:                    {summary['recall']:.4f}\\n")
        f.write(f"F1 Score:                  {summary['f1']:.4f}\\n")
        f.write(f"Genotype Concordance:      {summary['gt_concordance']:.4f}\\n")
        f.write(f"\\nBase VCF Count:            {summary['base cnt']}\\n")
        f.write(f"Comparison VCF Count:      {summary['comp cnt']}\\n")
    
    print(f"Summary extracted to {output_file}")
except Exception as e:
    print(f"Error extracting summary: {e}", file=sys.stderr)
    sys.exit(1)
PYTHON_EOF
  fi

done
echo "All SV benchmarking complete. Results in $RESULTSDIR"
