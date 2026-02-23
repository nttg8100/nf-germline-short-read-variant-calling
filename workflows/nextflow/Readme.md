# Nextflow GATK Variant Calling Pipeline

This directory contains the **Nextflow DSL2 implementation** of the 16-step GATK germline variant calling workflow. This is a production-ready, scalable pipeline with full container support and multi-sample parallelization.

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

### 1. Download Test Data

```bash
# From repository root
cd /path/to/variant-calling-gatk-pipeline-best-practice-from-scratch
pixi run bash scripts/download_data.sh
```

### 2. Run Single Sample (Test Profile)

```bash
cd workflows/nextflow

# Run with test profile (uses sample1 from data/)
nextflow run main.nf -profile singularity,test -resume
```

**Runtime**: ~3m 38s | **Outputs**: 93 files in `results_nextflow/`

### 3. Run Multiple Samples (Samplesheet)

```bash
# Edit samplesheet.csv to add/remove samples and lanes
cat samplesheet.csv
# sample,lane,fastq_1,fastq_2
# sample1,L001,../../data/sample1_L001_R1.fastq.gz,../../data/sample1_L001_R2.fastq.gz
# sample1,L002,../../data/sample1_L002_R1.fastq.gz,../../data/sample1_L002_R2.fastq.gz
# sample2,L001,../../data/sample2_L001_R1.fastq.gz,../../data/sample2_L001_R2.fastq.gz

# Run with samplesheet
nextflow run main.nf -profile singularity --input samplesheet.csv -resume
```

**Runtime**: 3m 26s for 3 samples (3.2x faster than bash!)

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
    cpus = 4
    memory = '16 GB'
}
```

### Custom Parameters

Override defaults via command line:

```bash
# Change output directory
nextflow run main.nf -profile singularity,test --outdir my_results

# Use custom reference genome
nextflow run main.nf -profile singularity,test --reference /path/to/genome.fasta

# Adjust thread count
nextflow run main.nf -profile singularity,test --threads 8
```

## Output Structure

After running the pipeline, results are organized by sample:

```
results_nextflow/
├── trimmed/
│   ├── sample1/
│   │   ├── sample1_L001_trimmed_1.fastq.gz + _2.fastq.gz
│   │   ├── sample1_L002_trimmed_1.fastq.gz + _2.fastq.gz
│   │   └── sample1_L001_fastp.html/json (QC reports)
│   ├── sample2/
│   └── sample3/
├── aligned/
│   ├── sample1/
│   │   ├── sample1_L001_aligned.bam (per-lane)
│   │   └── sample1_L002_aligned.bam
│   ├── sample2/
│   └── sample3/
├── sorted/
│   ├── sample1/
│   │   ├── sample1_L001_sorted.bam
│   │   └── sample1_L002_sorted.bam
│   ├── sample2/
│   └── sample3/
├── merged/
│   ├── sample1/sample1_merged.bam (lanes merged)
│   ├── sample2/
│   └── sample3/
├── markdup/
│   ├── sample1/sample1_dedup.bam + .bai
│   ├── sample2/
│   └── sample3/
├── bqsr/
│   ├── sample1/sample1_recal.table
│   ├── sample2/
│   └── sample3/
├── recalibrated/
│   ├── sample1/
│   │   ├── sample1_recal.bam + .bai
│   │   └── sample1_coverage.bedgraph
│   ├── sample2/
│   └── sample3/
├── qc/
│   ├── sample1/
│   │   ├── sample1_alignment_summary.txt
│   │   ├── sample1_insert_size_metrics.txt
│   │   └── sample1_insert_size_histogram.pdf
│   ├── sample2/
│   └── sample3/
└── variants/
    ├── sample1/
    │   ├── sample1.g.vcf.gz + .tbi
    │   ├── sample1_raw.vcf.gz + .tbi
    │   ├── sample1_raw_snps.vcf.gz + .tbi
    │   ├── sample1_raw_indels.vcf.gz + .tbi
    │   ├── sample1_filtered_snps.vcf.gz + .tbi
    │   ├── sample1_filtered_indels.vcf.gz + .tbi
    │   ├── sample1_filtered.vcf.gz + .tbi
    │   ├── sample1_annotated.vcf
    │   ├── stats/
    │   │   ├── sample1_variant_stats.txt
    │   │   ├── sample1_snp_count.txt
    │   │   └── sample1_indel_count.txt
    │   └── sample1_variants.bed
    ├── sample2/
    └── sample3/
```

## Container Strategy

### BioContainers (Most Tools)

```groovy
container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
container 'quay.io/biocontainers/bwa:0.7.18--he4a0461_2'
container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'
container 'quay.io/biocontainers/bcftools:1.17--haef29d1_0'
```

### Broad Institute GATK (Step 8 - CollectMetrics)

**Issue**: Standard BioContainers GATK lacks R runtime for PDF histogram generation.

**Solution**: Use Broad's official image:
```groovy
container 'broadinstitute/gatk:4.4.0.0'  // Includes R 4.2.x + plotting libraries
```

**Details**: See [Part 2 blog post, Section 10.3, Challenge #5](https://riverxdata.github.io/river-docs/blog/gatk-bash-nextflow-migration-md5-validation-part2)

### Container Cache

Singularity caches images in `$HOME/.singularity/cache` to avoid re-downloading:

```bash
# Check cached images
ls -lh $HOME/.singularity/cache/

# Clear cache if needed (forces re-download)
rm -rf $HOME/.singularity/cache/*
```

## Performance Benchmarks

### Single Sample (Test Profile)

| Metric | Value |
|--------|-------|
| **Total Runtime** | 3m 38s |
| **Processes** | 19 |
| **CPU hours** | 0.5 |
| **Output files** | 93 |
| **Cache hit rate** | 74.9% (with `-resume`) |

### Multi-Sample (3 Samples)

| Metric | Value |
|--------|-------|
| **Total Runtime** | 3m 26s |
| **Processes** | 57 (19 × 3 samples) |
| **CPU hours** | 1.3 |
| **Output files** | 279 (93 × 3) |
| **Speedup vs Bash** | **3.2x faster** (11 min → 3.4 min) |

### Parallelization Analysis

| Step | Parallelization | Notes |
|------|----------------|-------|
| Steps 1-9 | ✅ Per-sample | All samples run simultaneously |
| Step 10 | ❌ Single job | Joint genotyping requires all GVCFs |
| Steps 11-12 | ✅ Parallel branches | SNPs and indels processed concurrently |
| Steps 13-16 | ✅ Per-sample | Final outputs generated in parallel |

**Why not 3x speedup for 3 samples?**
- Step 10 (GenotypeGVCFs) is a serial bottleneck (must wait for all samples)
- Container startup overhead (~10-15s per process)
- Shared resource contention on test machine

**On HPC with more resources:**
- 10 samples: ~5-10 minutes (8-10x speedup)
- 100 samples: ~2-3 hours (40-50x speedup)

## Troubleshooting

### Issue 1: GATK CollectMetrics fails with "RScript not found"

**Error**:
```
A USER ERROR has occurred: RScript not found in environment path
```

**Cause**: Standard GATK container lacks R runtime for PDF histogram generation.

**Solution**: `modules/gatk_collectmetrics.nf` uses `broadinstitute/gatk:4.4.0.0` which includes R.

**Verify**:
```bash
# Check module container directive
grep "container" modules/gatk_collectmetrics.nf
# Should show: container 'broadinstitute/gatk:4.4.0.0'
```

### Issue 2: Work directory filling up disk space

**Cause**: Nextflow stores intermediate files in `work/` directory.

**Solution**: Clean up after successful runs:
```bash
# Remove work directory (keeps results, loses resume capability)
rm -rf work/

# Or use Nextflow's cleanup command
nextflow clean -f
```

**Note**: Only clean `work/` after confirming results are correct!

### Issue 3: Singularity image download fails

**Error**:
```
FATAL:   Unable to handle docker://quay.io/biocontainers/...
```

**Cause**: Network issues or Singularity version incompatibility.

**Solutions**:
```bash
# Check Singularity version (requires >=3.5)
singularity --version

# Test manual image pull
singularity pull docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0

# Clear Singularity cache and retry
rm -rf $HOME/.singularity/cache/*
nextflow run main.nf -profile singularity,test -resume
```

### Issue 4: "No such file or directory" for BWA index files

**Cause**: BWA requires all 5 index files (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) staged in work directory.

**Solution**: `main.nf` explicitly channels all index files:
```groovy
bwa_index_ch = Channel.of([
    file("${params.reference}.amb"),
    file("${params.reference}.ann"),
    file("${params.reference}.bwt"),
    file("${params.reference}.pac"),
    file("${params.reference}.sa")
]).collect()
```

**Verify**: Ensure all index files exist in `reference/`:
```bash
ls -lh ../../reference/genome.fasta*
```

### Issue 5: Zero variants in output VCF

**Cause**: Test data is chr22 subset with low coverage.

**Expected**: This is normal for test data. Both bash and Nextflow should report 0 variants.

**Validation**:
```bash
# Compare variant counts between bash and Nextflow
cat ../../workflows/bash/results/variants/sample1_filtered_snp_count.txt
cat results_nextflow/variants/sample1/sample1_filtered_snp_count.txt
# Both should show: 0
```

### Issue 6: Process crashes with "OutOfMemoryError"

**Cause**: JVM heap size too small for GATK processes.

**Solution**: Increase memory in process directives:
```groovy
// In main.nf or nextflow.config
process {
    withName: GATK_HAPLOTYPECALLER {
        memory = '16 GB'  // Increase from default 8 GB
    }
}
```

## Validation Against Bash

To verify scientific equivalence between Bash and Nextflow implementations:

```bash
# Run bash pipeline
cd ../bash
pixi run bash gatk_pipeline.sh 2>&1 | tee bash.log

# Run nextflow pipeline
cd ../nextflow
nextflow run main.nf -profile singularity,test -resume

# Compare variant statistics (should match exactly)
diff ../bash/results/variants/sample1_variant_stats_filtered.txt \
     results_nextflow/variants/sample1/sample1_variant_stats_filtered.txt

# Expected: SNP/indel counts identical
# MD5 checksums will differ (timestamps in headers - this is normal!)
```

**Key validation metrics:**
- ✅ Raw SNP count: Should match
- ✅ Raw indel count: Should match
- ✅ Filtered SNP count (PASS): Should match
- ✅ Filtered indel count (PASS): Should match
- ⚠️ MD5 checksums: Will differ (expected due to timestamps)

## Advanced Usage

### Resume from Failures

```bash
# Pipeline fails at step 10
nextflow run main.nf -profile singularity,test

# Fix issue, then resume (skips completed steps)
nextflow run main.nf -profile singularity,test -resume
```

### Custom Samplesheet Format

```csv
sample_id,fastq_1,fastq_2,population,sequencing_center
sample1,/data/s1_R1.fq.gz,/data/s1_R2.fq.gz,EUR,BGI
sample2,/data/s2_R1.fq.gz,/data/s2_R2.fq.gz,AFR,Illumina
sample3,/data/s3_R1.fq.gz,/data/s3_R2.fq.gz,EAS,BGI
```

Modify `main.nf` to parse additional columns:
```groovy
Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row ->
        def meta = [
            id: row.sample_id,
            population: row.population,
            center: row.sequencing_center
        ]
        def reads = [file(row.fastq_1), file(row.fastq_2)]
        return [meta, reads]
    }
```

### Run Specific Steps Only

```bash
# Run only QC steps (Steps 1-2)
nextflow run main.nf -profile singularity,test -entry QC_ONLY

# Run only variant calling (Steps 9-13)
nextflow run main.nf -profile singularity,test -entry VARIANT_CALLING
```

(Requires defining multiple `workflow` entries in `main.nf`)

## Next Steps (Part 3)

- [ ] **SLURM integration** - Deploy on HPC clusters
- [ ] **1000 Genomes validation** - Test with real 30x coverage data
- [ ] **MultiQC reporting** - Aggregate QC across 100+ samples
- [ ] **Cloud deployment** - AWS Batch, Google Cloud Life Sciences
- [ ] **nf-test framework** - Unit/integration testing
- [ ] **nf-core compliance** - Align with best practices

## References

- **Blog post**: [Part 2 - Migrating GATK Bash to Nextflow with MD5 Validation](https://riverxdata.github.io/river-docs/blog/gatk-bash-nextflow-migration-md5-validation-part2)
- **Nextflow docs**: https://www.nextflow.io/docs/latest/
- **nf-core best practices**: https://nf-co.re/docs/contributing/guidelines
- **BioContainers registry**: https://biocontainers.pro/
- **Repository root README**: `../../README.md`

## Support

For issues or questions:
- Open an issue on [GitHub](https://github.com/nttg8100/variant-calling-gatk-pipeline-best-practice-from-scratch/issues)
- Check troubleshooting guide: `../../docs/TROUBLESHOOTING.md`
- Read the blog post: https://riverxdata.github.io/river-docs/

---

**Status**: ✅ Production Ready (16/16 steps implemented)  
**Last Updated**: February 2026
