# Nextflow GATK Variant Calling Pipeline

This directory contains the **Nextflow DSL2 implementation** of the 16-step GATK germline variant calling workflow. This is a production-ready, scalable pipeline with full container support and multi-sample parallelization.

## What's New (February 2026)

This implementation includes significant optimizations over the baseline bash workflow:

- ✅ **fastp** replaces FastQC + Trim Galore (1.5-2x faster, single-pass QC/trimming)
- ✅ **BWA-MEM2** replaces BWA-MEM (2-3x faster alignment with AVX2 optimization)
- ✅ **GATK Spark** for 3 steps: MarkDuplicates, BaseRecalibrator, ApplyBQSR (1.3-1.5x faster)
- ✅ **Multi-lane support** - Automatically process and merge multiple sequencing lanes per sample
- ✅ **Value channel reuse** - Eliminates ~25 duplicate `.collect()` calls for reference files
- ✅ **Centralized configuration** - All `publishDir` settings in `conf/modules.config`

## Purpose

✅ **Production workloads** - Process 100+ samples in parallel  
✅ **HPC/Cloud deployment** - Run on SLURM, AWS Batch, Google Cloud  
✅ **Container reproducibility** - Singularity/Docker for full portability  
✅ **Resume capability** - Restart from failures with `-resume`  
✅ **Automatic parallelization** - 3-7x faster than serial bash execution

## Features

- **20 modular processes** - One DSL2 module per logical step
- **Multi-sample, multi-lane support** - CSV samplesheet with lane-level processing
- **BWA-MEM2** - Faster alignment with AVX2 optimization (2-3x speedup over BWA)
- **fastp** - All-in-one QC, trimming, and filtering (replaces FastQC + Trim Galore)
- **GATK Spark** - Multi-threaded BQSR processing (MarkDuplicates, BaseRecalibrator, ApplyBQSR)
- **Centralized configuration** - All publishDir settings in `conf/modules.config`
- **Container strategy** - BioContainers for reproducibility
- **Flexible profiles** - Test (single sample) / Full (multi-sample) / SLURM (HPC)

## Quick Start

### Prerequisites

From the **repository root**, ensure:

```bash
# Tools installed via Pixi
pixi install
pixi run nextflow -version

# Singularity available for containers
singularity --version  # or apptainer --version
```
### 1. Run Examples Datasets

```bash
cd workflows/nextflow

# Run with test profile (uses sample1 from data/)
nextflow run main.nf -profile singularity,test -resume
```

**Runtime**: ~3m 38s | **Outputs**: 93 files in `results_nextflow/`
**Multi-lane support**: The pipeline automatically processes each lane independently through Steps 1-3 (fastp, BWA-MEM2, sort), then merges all lanes per sample at Step 4 (samtools merge). Subsequent steps operate on merged per-sample BAMs.

## Complete Pipeline Overview

This pipeline implements the **GATK germline variant calling best practices** (16 steps) with key optimizations:

### Pre-processing (Steps 1-8)
1. **fastp** - All-in-one QC, adapter trimming, and filtering (replaces FastQC + Trim Galore)
2. **BWA-MEM2** - Read alignment with AVX2 optimization (2-3x faster than BWA)
3. **SAMtools sort** - Sort BAM by coordinate (per-lane processing)
4. **SAMtools merge** - Merge multi-lane BAMs per sample
5. **GATK Spark MarkDuplicates** - Mark PCR/optical duplicates (multi-threaded)
6. **GATK Spark BaseRecalibrator** - Model systematic errors (BQSR, multi-threaded)
7. **GATK Spark ApplyBQSR** - Apply recalibration (multi-threaded)
8. **GATK CollectMetrics** - Alignment QC + insert size histogram

### Variant Calling (Steps 9-10)
9. **GATK HaplotypeCaller** - Call variants (GVCF mode)
10. **GATK GenotypeGVCFs** - Joint genotyping across samples

### Variant Filtering (Steps 11-13)
11. **SelectVariants + VariantFiltration (SNPs)** - Filter SNPs (QUAL, QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum)
12. **SelectVariants + VariantFiltration (Indels)** - Filter indels (QUAL, QD, FS, ReadPosRankSum)
13. **GATK MergeVcfs** - Merge filtered SNPs + indels

### Annotation & Statistics (Steps 14-16)
14. **SnpEff** - Functional annotation (gene, transcript, effect)
15. **bcftools stats** - Variant statistics (raw/filtered counts)
16. **Visualization** - BED file + bedGraph coverage track

## Pipeline Architecture

### Workflow Structure

The main workflow (`main.nf`) orchestrates 16 GATK steps across 20 Nextflow processes:

```
main.nf
├── Step 1: QC, Trimming, and Filtering (parallel across all lane-level FASTQs)
│   └── FASTP(reads_ch) → Trimmed FASTQ + HTML/JSON reports
├── Step 2: Read Alignment (parallel across all lanes)
│   └── BWA_MEM2(trimmed_reads, reference) → Per-lane BAM
├── Step 3: Sort BAM files (parallel)
│   └── SAMTOOLS_SORT(bam) → Sorted BAM per lane
├── Step 4: Merge lanes per sample
│   └── SAMTOOLS_MERGE(sorted_bams.groupBy(sample)) → Per-sample BAM
├── Step 5: Mark Duplicates (parallel across samples)
│   └── GATKSPARK_MARKDUPLICATES(merged_bam) → Dedup BAM
├── Step 6-7: Base Quality Score Recalibration (parallel across samples)
│   ├── GATKSPARK_BASERECALIBRATOR(dedup_bam, ref, known_sites) → Recal table
│   └── GATKSPARK_APPLYBQSR(dedup_bam, recal_table, ref) → Recalibrated BAM
├── Step 8: Alignment Quality Assessment (parallel)
│   └── GATK_COLLECTMETRICS(recal_bam, ref) → QC metrics
├── Step 9: Variant Calling (parallel across samples)
│   └── GATK_HAPLOTYPECALLER(recal_bam, ref) → GVCF
├── Step 10: Joint Genotyping (single job, all samples)
│   └── GATK_GENOTYPEGVCFS(all_gvcfs.collect(), ref) → Raw VCF
├── Step 11-12: Filtering (parallel branches: SNPs and Indels)
│   ├── GATK_SELECTVARIANTS_SNP + GATK_VARIANTFILTRATION_SNP
│   └── GATK_SELECTVARIANTS_INDEL + GATK_VARIANTFILTRATION_INDEL
├── Step 13: Merge Filtered Variants
│   └── GATK_MERGEVCFS(snp_vcf, indel_vcf) → Filtered VCF
├── Step 14: Functional Annotation
│   └── SNPEFF(merged_vcf) → Annotated VCF
├── Step 15: Variant Statistics
│   └── BCFTOOLS_STATS(filtered_vcf) → Variant stats
└── Step 16: Visualization
    ├── BCFTOOLS_QUERY(merged_vcf) → BED
    └── BEDTOOLS_GENOMECOV(recal_bam) → bedGraph
```

**Key Optimizations:**
- **Lane-aware processing**: Supports multiple lanes per sample, automatically merges at Step 4
- **GATK Spark parallelization**: Steps 5-7 use multi-threaded Spark implementations
- **Value channel reuse**: Reference files collected once and reused (eliminates duplicate `.collect()` calls)

### Modules Directory

```
modules/
├── fastp.nf                           # Step 1: QC, trimming, filtering (replaces FastQC + Trim Galore)
├── bwa_mem2.nf                        # Step 2: Read alignment with BWA-MEM2
├── samtools_sort.nf                   # Step 3: BAM sorting
├── samtools_merge.nf                  # Step 4: Merge per-lane BAMs by sample
├── gatkspark_markduplicates.nf        # Step 5: PCR duplicate marking (Spark)
├── gatkspark_baserecalibrator.nf      # Step 6: BQSR table generation (Spark)
├── gatkspark_applybqsr.nf             # Step 7: Apply recalibration (Spark)
├── gatk_collectmetrics.nf             # Step 8: Alignment QC
├── gatk_haplotypecaller.nf            # Step 9: Variant calling (GVCF mode)
├── gatk_genotypegvcfs.nf              # Step 10: Joint genotyping
├── gatk_selectvariants_snp.nf         # Step 11a: Extract SNPs
├── gatk_variantfiltration_snp.nf      # Step 11b: Filter SNPs
├── gatk_selectvariants_indel.nf       # Step 12a: Extract indels
├── gatk_variantfiltration_indel.nf    # Step 12b: Filter indels
├── gatk_mergevcfs.nf                  # Step 13: Merge filtered VCFs
├── snpeff.nf                          # Step 14: Functional annotation
├── bcftools_stats.nf                  # Step 15: Variant statistics
├── bcftools_query.nf                  # Step 16a: VCF to BED conversion
├── bedtools_genomecov.nf              # Step 16b: Coverage bedGraph
└── save_reference.nf                  # Utility: Save reference files from S3
```

Each module follows **nf-core DSL2 conventions**:
- `tag` directive for sample tracking
- `publishDir` configurations centralized in `conf/modules.config`
- `container` directive for Singularity/Docker images
- `input/output/script` blocks with clear semantics

## Configuration

The `nextflow.config` file defines profiles and resource allocation:

### Profiles

**test** - Single sample, quick validation
```groovy
params {
    input = null  // Uses default sample1 from data/
    outdir = 'results_nextflow'
    reference = '../../reference/genome.fasta'
    // ... more params
}
```

**full** - Multi-sample production
```groovy
params {
    input = 'samplesheet.csv'  // CSV with 100+ samples
    outdir = 'results_production'
    // ... more params
}
```

**singularity** - Container execution (recommended)
```groovy
singularity {
    enabled = true
    autoMounts = true
    cacheDir = "$HOME/.singularity/cache"
}
```

**slurm** - HPC cluster (coming in Part 3)
```groovy
process {
    executor = 'slurm'
    queue = 'normal'
}
```

### Custom Parameters

Override defaults via command line:

```bash
# Change output directory
nextflow run main.nf -profile singularity,test --outdir my_results
```

## Output Structure

After running the pipeline, results are organized by lowercases of process name


## References

- **Main repository**: [../../README.md](../../README.md) - Project overview, installation, quick start
- **Bash implementation**: [../bash/README.md](../bash/README.md) - Educational single-sample workflow
- **Blog post**: [Part 2 - Migrating GATK Bash to Nextflow with MD5 Validation](https://riverxdata.github.io/river-docs/blog/gatk-bash-nextflow-migration-md5-validation-part2)
- **Nextflow docs**: https://www.nextflow.io/docs/latest/
- **nf-core best practices**: https://nf-co.re/docs/contributing/guidelines
- **BioContainers registry**: https://biocontainers.pro/

## Support

For issues or questions:
- Open an issue on [GitHub](https://github.com/nttg8100/variant-calling-gatk-pipeline-best-practice-from-scratch/issues)
- Check troubleshooting guide: `../../docs/TROUBLESHOOTING.md`
- Read the blog post: https://riverxdata.github.io/river-docs/