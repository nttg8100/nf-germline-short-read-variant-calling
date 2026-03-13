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

```bash
pixi run nextflow run main.nf -profile docker,test -resume 
```

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

- [Pipeline Architecture](docs/PIPELINE_ARCHITECTURE.md) - Detailed documentation with tool versions and parameters
- [Benchmarking Guide](benchmark/Readme.md) - SV and SNP/INDEL benchmarking instructions with Truvari metrics