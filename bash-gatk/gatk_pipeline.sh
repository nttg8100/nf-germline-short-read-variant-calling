#!/bin/bash
#
# GATK Variant Calling Pipeline - Complete Best Practices (10 Steps)
# Based on NGS101 tutorial: https://ngs101.com/how-to-analyze-whole-genome-sequencing-data-for-absolute-beginners-part-1-from-raw-reads-to-high-quality-variants-using-gatk/
# Adapted to use nf-core test data with hard filtering (Option 2)
#
# This pipeline implements the complete GATK workflow with:
# - Adapter trimming and quality filtering (Trim Galore)
# - Read alignment with BWA-MEM
# - Duplicate marking (GATK MarkDuplicates)
# - Base quality score recalibration (BQSR)
# - Alignment quality metrics (GATK CollectMetrics)
# - Variant calling with HaplotypeCaller (GVCF mode)

set -euo pipefail

# Configuration
SAMPLE="${1:-sample1}"  # Accept sample name as first argument, default to sample1
THREADS=2
DATA_DIR="data"
REF_DIR="reference"
RESULTS_DIR="results"

# Reference files
REFERENCE="${REF_DIR}/genome.fasta"
DBSNP="${REF_DIR}/dbsnp_146.hg38.vcf.gz"
KNOWN_INDELS="${REF_DIR}/mills_and_1000G.indels.vcf.gz"

# Input FASTQ files
FASTQ_R1="${DATA_DIR}/${SAMPLE}_R1.fastq.gz"
FASTQ_R2="${DATA_DIR}/${SAMPLE}_R2.fastq.gz"

# Output directories
QC_DIR="${RESULTS_DIR}/qc/${SAMPLE}"
TRIMMED_DIR="${RESULTS_DIR}/trimmed/${SAMPLE}"
ALIGNED_DIR="${RESULTS_DIR}/aligned/${SAMPLE}"
VAR_DIR="${RESULTS_DIR}/var/${SAMPLE}"

# Output files
TRIMMED_R1="${TRIMMED_DIR}/${SAMPLE}_R1_val_1.fq.gz"
TRIMMED_R2="${TRIMMED_DIR}/${SAMPLE}_R2_val_2.fq.gz"
ALIGNED_BAM="${ALIGNED_DIR}/${SAMPLE}.bam"
SORTED_BAM="${ALIGNED_DIR}/${SAMPLE}.sorted.bam"
DEDUP_BAM="${ALIGNED_DIR}/${SAMPLE}_marked_duplicates.bam"
RECAL_TABLE="${ALIGNED_DIR}/${SAMPLE}_recal_data.table"
RECAL_BAM="${ALIGNED_DIR}/${SAMPLE}_recalibrated.bam"
GVCF="${VAR_DIR}/${SAMPLE}.g.vcf.gz"
RAW_VCF="${VAR_DIR}/${SAMPLE}_raw_variants.vcf.gz"
METRICS="${ALIGNED_DIR}/${SAMPLE}_duplicate_metrics.txt"

# Check tools are available
echo "Checking required tools..."
echo "✓ FastQC:" 
fastqc --version
echo "✓ Trim Galore:" 
trim_galore --version
echo "✓ BWA:" 
bwa || true  # BWA doesn't have a --version flag, so we ignore the error
echo "✓ Samtools:" 
samtools --version
echo "✓ GATK:" 
gatk --version
echo "✓ R (for insert size histogram):" 
R --version
echo "All required tools are available."

# Create output directories
mkdir -p ${QC_DIR} ${TRIMMED_DIR} ${ALIGNED_DIR} ${VAR_DIR}

echo "=========================================="
echo "GATK Variant Calling Pipeline - Complete"
echo "Sample: ${SAMPLE}"
echo "Based on NGS101 Best Practices"
echo "=========================================="
echo

# Step 1: Quality Control with FastQC
echo "[$(date)] Step 1: Running FastQC on raw reads..."
fastqc -o ${QC_DIR} -t ${THREADS} ${FASTQ_R1} ${FASTQ_R2}
echo "[$(date)] FastQC completed"
echo

# Step 2: Adapter Trimming and Quality Filtering
echo "[$(date)] Step 2: Adapter trimming with Trim Galore..."
trim_galore \
    --paired \
    --quality 20 \
    --length 50 \
    --fastqc \
    --output_dir ${TRIMMED_DIR} \
    ${FASTQ_R1} ${FASTQ_R2}
echo "[$(date)] Adapter trimming completed"
echo

# Step 3: Read Alignment with BWA-MEM
echo "[$(date)] Step 3: Aligning reads with BWA-MEM..."
bwa mem \
    -t ${THREADS} \
    -M \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}_lib" \
    ${REFERENCE} \
    ${TRIMMED_R1} ${TRIMMED_R2} \
    | samtools view -Sb - > ${ALIGNED_BAM}
echo "[$(date)] Alignment completed"
echo

# Step 4: Sort and Index BAM file
echo "[$(date)] Step 4: Sorting BAM file..."
samtools sort -@ ${THREADS} -o ${SORTED_BAM} ${ALIGNED_BAM}
samtools index ${SORTED_BAM}
echo "[$(date)] Sorting completed"
echo

# Step 5: Mark Duplicates
echo "[$(date)] Step 5: Marking duplicates with GATK..."
gatk MarkDuplicates \
    -I ${SORTED_BAM} \
    -O ${DEDUP_BAM} \
    -M ${METRICS} 
echo "[$(date)] Mark duplicates completed"
echo

# Step 6: Base Quality Score Recalibration (BQSR) - Generate table
echo "[$(date)] Step 6: Generating BQSR recalibration table..."
gatk BaseRecalibrator \
    -I ${DEDUP_BAM} \
    -R ${REFERENCE} \
    --known-sites ${DBSNP} \
    --known-sites ${KNOWN_INDELS} \
    -O ${RECAL_TABLE}
echo "[$(date)] BQSR table generated"
echo

# Step 7: Apply BQSR
echo "[$(date)] Step 7: Applying BQSR..."
gatk ApplyBQSR \
    -I ${DEDUP_BAM} \
    -R ${REFERENCE} \
    --bqsr-recal-file ${RECAL_TABLE} \
    -O ${RECAL_BAM}
echo "[$(date)] BQSR applied"
echo

# Step 8: Alignment Quality Assessment
echo "[$(date)] Step 8: Collecting alignment quality metrics..."
gatk CollectAlignmentSummaryMetrics \
    -R ${REFERENCE} \
    -I ${RECAL_BAM} \
    -O ${QC_DIR}/${SAMPLE}_alignment_summary.txt

# Collect insert size metrics with histogram (requires R)
gatk CollectInsertSizeMetrics \
    -I ${RECAL_BAM} \
    -O ${QC_DIR}/${SAMPLE}_insert_size_metrics.txt \
    -H ${QC_DIR}/${SAMPLE}_insert_size_histogram.pdf
echo "[$(date)] Alignment QC completed"
echo

# Step 9: Variant Calling with HaplotypeCaller (GVCF mode)
echo "[$(date)] Step 9: Calling variants with HaplotypeCaller in GVCF mode..."
gatk HaplotypeCaller \
    -R ${REFERENCE} \
    -I ${RECAL_BAM} \
    -O ${GVCF} \
    -ERC GVCF \
    --dbsnp ${DBSNP}
echo "[$(date)] GVCF generation completed"
echo

# Step 10: Genotype GVCFs
echo "[$(date)] Step 10: Genotyping GVCF to VCF..."
gatk GenotypeGVCFs \
    -R ${REFERENCE} \
    -V ${GVCF} \
    -O ${RAW_VCF}
echo "[$(date)] Genotyping completed"
echo
