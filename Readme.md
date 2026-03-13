# Nextflow Short Reads Germline Variant Calling

This repository implements a comprehensive Nextflow pipeline for germline short-read variant calling, supporting multiple variant callers (DeepVariant, GATK HaplotypeCaller, FreeBayes) with integrated quality control and annotation.

## Primary Use Case

**Primary support**: **30X whole-genome sequencing (WGS) Illumina short reads with GRCh38 (hg38) alignment**

The pipeline is optimized for the following workflow:

- **Input**: Illumina short-read FASTQ files (~30X coverage)
- **Quality Control**: FASTP for read-level quality filtering
- **Alignment**: BWA-MEM2 alignment to GRCh38 (hg38) reference genome
- **Preprocessing**: Skipped by default (can be enabled if needed)
- **Variant Calling**: **DeepVariant** (default) for high-accuracy germline variant detection
- **Quality Reports**: 
  - Alignment quality metrics (samtools stats, FastQC, MultiQC)
  - Variant call quality reports
- **Variant Annotation**: 
  - **SnpEff** (annotation database)
  - **VEP** (Ensembl Variant Effect Predictor)
- **Output**: Annotated VCF files with comprehensive quality metrics

### SNP and INDEL Benchmarking Performance

This pipeline is specifically optimized for **SNP and small INDEL detection** using **DeepVariant**. Benchmark results against HG002 (Genome in a Bottle) demonstrate competitive performance compared to nf-core/sarek:

#### SNP Detection Performance (HG002)

| Caller | Pipeline | Recall | Precision | F1 Score | TP | FN |
|--------|----------|--------|-----------|----------|-----|-----|
| **DeepVariant** | **nf-germline-short-read** | **99.39%** | **99.82%** | **99.60%** | 3,344,672 | 20,455 |
| GATK HaplotypeCaller | nf-core/sarek | 99.26% | 99.16% | 99.21% | 3,340,110 | 25,017 |
| FreeBayes | nf-core/sarek | 99.39% | 96.29% | 97.81% | 3,344,473 | 20,654 |

#### INDEL Detection Performance (HG002)

| Caller | Pipeline | Recall | Precision | F1 Score | TP | FN |
|--------|----------|--------|-----------|----------|--------|---------|
| **DeepVariant** | **nf-germline-short-read** | **98.78%** | **99.38%** | **99.08%** | 519,079 | 6,390 |
| GATK HaplotypeCaller | nf-core/sarek | 98.37% | 98.83% | 98.60% | 516,926 | 8,543 |
| FreeBayes | nf-core/sarek | 95.99% | 96.40% | 96.20% | 504,432 | 21,037 |

**Key Results**: DeepVariant in nf-germline-short-read achieves the highest SNP accuracy (99.60% F1) and excellent INDEL accuracy (99.08% F1), outperforming GATK HaplotypeCaller and FreeBayes from nf-core/sarek on the same HG002 test dataset.

### Configuration for Primary Use Case

```bash
# Default configuration uses:
# - DeepVariant as variant caller
# - HG38/GRCh38 reference genome
# - FASTP quality filtering
# - Preprocessing skipped (skip_preprocessing: true)
# - SnpEff + VEP annotation enabled

pixi run nextflow run main.nf -profile docker -resume
```

## Pipeline Architecture

![Nextflow Germline Variant Calling Pipeline](docs/nf-germline-pipeline.png)

For a detailed breakdown of the pipeline architecture, tool versions, parameters, and usage examples, see [PIPELINE_ARCHITECTURE.md](docs/PIPELINE_ARCHITECTURE.md).

## Quick Start

### 1. Prepare a Samplesheet

Create a CSV samplesheet with your input FASTQ files. The samplesheet format is:

```csv
sample,lane,fastq_1,fastq_2
HG002,1,/path/to/HG002_R1.fastq.gz,/path/to/HG002_R2.fastq.gz
HG003,1,/path/to/HG003_R1.fastq.gz,/path/to/HG003_R2.fastq.gz
HG004,1,/path/to/HG004_R1.fastq.gz,/path/to/HG004_R2.fastq.gz
```

**Columns**:
- `sample`: Sample identifier (used for output file naming)
- `lane`: Sequencing lane (if multiple lanes, create separate rows per lane)
- `fastq_1`: Path to first-in-pair (R1) FASTQ file (gzipped)
- `fastq_2`: Path to second-in-pair (R2) FASTQ file (gzipped)

### 2. Run the Pipeline

```bash
pixi run nextflow run main.nf \
  --input samplesheet.csv \
  --profile docker \
  -resume
```

For test mode with sample data:

```bash
pixi run nextflow run main.nf -profile docker,test -resume 
```

### 3. View Results

Output files will be generated in the `results/` directory:
- Aligned BAM files
- Variant calls (VCF format)
- Annotated VCF files (SnpEff + VEP)
- Quality control reports (MultiQC)

For more advanced usage and configuration options, see the [Pipeline Architecture](docs/PIPELINE_ARCHITECTURE.md) documentation.

## Key Features

- **Multiple Variant Callers**: DeepVariant (default), GATK HaplotypeCaller, FreeBayes
- **Quality Control**: FastQC, MultiQC, samtools stats
- **Variant Annotation**: SnpEff, VEP
- **Structural Variants**: Manta SV calling (with benchmarking support)
- **Flexible Configuration**: Container support (Docker/Singularity), multiple profiles
- **Benchmarking Tools**: Integrated Truvari for SV benchmarking

## Benchmarking Results

This simple, streamlined workflow achieves performance comparable to **nf-core/sarek** for short-read germline variant calling:

### SNP/INDEL Variant Calling (Sarek Comparison)
- **DeepVariant**: High sensitivity and precision for SNP/INDEL detection
- **GATK HaplotypeCaller**: Robust multi-sample calling capability
- **FreeBayes**: Excellent for population-level variant discovery

### Structural Variant Calling (Manta)
- **HG002 Benchmarking Results**:
  - Sensitivity: 7.88% (1,082 TP out of 13,732 variants)
  - Precision: 36.8% (1,082 TP out of 2,940 calls)
  - Genotype Concordance: 90.39%

For detailed benchmarking methodology, results, and comparisons, see [Benchmarking Guide](benchmark/Readme.md).

## Documentation

- [Pipeline Architecture](docs/architecture.md) - Detailed documentation with tool versions and parameters
- [Benchmarking Guide](benchmark/Readme.md) - SV and SNP/INDEL benchmarking instructions with Truvari metrics

## License
MIT ([LICENSE](LICENSE))