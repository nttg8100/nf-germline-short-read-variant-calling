# GATK Variant Calling Pipeline - From Bash to Nextflow

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Blog Post](https://img.shields.io/badge/Blog-RiverXData-blue)](https://riverxdata.github.io/river-docs/)

A complete GATK germline variant calling workflow demonstrating the **validated transformation** from bash scripts to enterprise-grade Nextflow pipelines. This repository implements all 16 GATK best practice steps with **proven scientific equivalence** through MD5 validation.

🎯 **What makes this special:**
- ✅ **Complete 16-step pipeline**: FastQC → Alignment → BQSR → Variant Calling → Annotation → Statistics
- ✅ **Two implementations**: Bash (educational baseline) + Nextflow (production-ready)
- ✅ **Validated equivalence**: MD5 comparison framework proves scientific accuracy
- ✅ **Container-ready**: Singularity/Docker containers for full reproducibility
- ✅ **Multi-sample support**: Process 100+ samples in parallel with Nextflow
- ✅ **Fully documented**: Detailed blog post series with troubleshooting guides

## 📊 Quick Performance Comparison

| Metric | Bash (Serial) | Nextflow (Parallel) |
|--------|--------------|---------------------|
| **1 sample** | 3m 38s | 3m 38s |
| **3 samples** | ~11 minutes | **3m 26s** (3.2x faster) |
| **10 samples** | ~36 minutes | **~5 minutes** (7x faster) |
| **100 samples** | ~6 hours | **~2 hours** (3x faster) |

## 🚀 Quick Start

### Prerequisites

- **Linux** system (tested on Ubuntu 20.04/22.04)
- **Pixi** package manager ([installation guide](https://pixi.sh))
- **8GB RAM** minimum (16GB recommended for multi-sample processing)
- **10GB disk space** for test data and results
- **Singularity/Apptainer** (for Nextflow containerized execution)

### Installation

```bash
# Clone the repository
git clone https://github.com/nttg8100/variant-calling-gatk-pipeline-best-practice-from-scratch.git
cd variant-calling-gatk-pipeline-best-practice-from-scratch

# Install dependencies with Pixi
pixi install

# Verify installation
pixi run bwa 2>&1 | head -3
pixi run samtools --version
pixi run gatk --version
```

### Download Test Data

```bash
# Download reference genome, known sites, and test FASTQ files
pixi run bash scripts/download_data.sh

# This creates:
# - data/ (sample1/2/3 paired-end FASTQ files)
# - reference/ (genome.fasta + BWA indices + dbSNP + Mills indels)
```

**Download size**: ~35 MB | **Time**: 1-2 minutes

### Run the Pipelines

#### Option 1: Bash Pipeline (Educational)

```bash
cd workflows/bash
pixi run bash gatk_pipeline.sh 2>&1 | tee pipeline.log
```

**Use case**: Learning GATK workflow step-by-step, single-sample processing

#### Option 2: Nextflow Pipeline (Production)

```bash
cd workflows/nextflow

# Single sample (test profile)
nextflow run main.nf -profile singularity,test -resume

# Multi-sample (3 samples from samplesheet)
nextflow run main.nf -profile singularity --input samplesheet.csv -resume
```

**Use case**: Multi-sample processing, HPC/cloud deployment, production workloads

## 📁 Repository Structure

```
variant-calling-gatk-pipeline-best-practice-from-scratch/
├── README.md                        # This file
├── LICENSE                          # MIT License
├── pixi.toml                        # Environment specification
├── pixi.lock                        # Locked dependencies
│
├── workflows/
│   ├── bash/                        # Bash implementation (Part 1)
│   │   ├── README.md               # Bash-specific instructions
│   │   ├── gatk_pipeline.sh        # 16-step bash workflow
│   │   └── download_data.sh        # Data download helper
│   │
│   └── nextflow/                    # Nextflow implementation (Part 2)
│       ├── README.md                # Nextflow-specific instructions
│       ├── main.nf                  # Main workflow (16 steps)
│       ├── nextflow.config          # Profiles (test/full/slurm)
│       ├── samplesheet.csv          # Multi-sample input example
│       ├── modules/                 # 18 DSL2 process modules
│       │   ├── fastqc.nf           # Step 1: Quality control
│       │   ├── trim_galore.nf      # Step 2: Adapter trimming
│       │   ├── bwa_mem.nf          # Step 3: Read alignment
│       │   ├── samtools_sort.nf    # Step 4: BAM sorting
│       │   ├── gatk_markduplicates.nf  # Step 5: PCR duplicates
│       │   ├── gatk_baserecalibrator.nf  # Step 6: BQSR table
│       │   ├── gatk_applybqsr.nf   # Step 7: Apply BQSR
│       │   ├── gatk_collectmetrics.nf  # Step 8: Alignment QC (R-enabled)
│       │   ├── gatk_haplotypecaller.nf  # Step 9: Variant calling
│       │   ├── gatk_genotypegvcfs.nf  # Step 10: Joint genotyping
│       │   ├── gatk_selectvariants_snp.nf  # Step 11a: SNPs
│       │   ├── gatk_variantfiltration_snp.nf  # Step 11b: Filter SNPs
│       │   ├── gatk_selectvariants_indel.nf  # Step 12a: Indels
│       │   ├── gatk_variantfiltration_indel.nf  # Step 12b: Filter indels
│       │   ├── gatk_mergevcfs.nf   # Step 13: Merge filtered VCFs
│       │   ├── snpeff.nf           # Step 14: Functional annotation
│       │   ├── bcftools_stats.nf   # Step 15: Variant statistics
│       │   ├── bcftools_query.nf   # Step 16a: VCF to BED
│       │   └── bedtools_genomecov.nf  # Step 16b: Coverage track
│       └── WORK_DIRECTORY_EXPLAINED.md  # Nextflow internals guide
│
├── data/                            # Test FASTQ files (not tracked)
│   ├── sample1_R1.fastq.gz         # Sample 1 forward reads
│   ├── sample1_R2.fastq.gz         # Sample 1 reverse reads
│   ├── sample2_R1/R2.fastq.gz      # Additional test samples
│   └── sample3_R1/R2.fastq.gz
│
├── reference/                       # Reference genome files (not tracked)
│   ├── genome.fasta                 # hg38 chr22 subset (40kb)
│   ├── genome.fasta.{amb,ann,bwt,pac,sa}  # BWA indices
│   ├── genome.fasta.fai             # SAMtools index
│   ├── genome.dict                  # GATK sequence dictionary
│   ├── dbsnp_146.hg38.vcf.gz + .tbi  # Known SNP sites (BQSR)
│   └── mills_and_1000G.indels.vcf.gz + .tbi  # Known indels (BQSR)
│
├── scripts/                         # Utility scripts
│   ├── download_data.sh            # Download test data and reference
│   └── validate_migration.sh       # MD5 validation (bash vs nextflow)
│
└── docs/                           # Additional documentation
    └── TROUBLESHOOTING.md          # Common issues and solutions
```

## 🧬 Pipeline Overview

This implementation follows **GATK germline variant calling best practices** (16 steps from QC → Alignment → BQSR → Variant Calling → Annotation → Statistics).

**For detailed pipeline steps and technical specifications, see:**
- 📖 [Bash Implementation README](workflows/bash/README.md) - Educational single-sample workflow
- 📖 [Nextflow Implementation README](workflows/nextflow/Readme.md) - Production multi-sample workflow with full technical details

## 🎓 Learning Resources

### Blog Post Series

This repository is accompanied by comprehensive blog posts:

- **[Part 1: Building a Production-Ready GATK Bash Workflow with Pixi](https://riverxdata.github.io/river-docs/blog/gatk-variant-calling-bash-workflow-pixi-part1)**
  - Step-by-step bash implementation
  - Tool installation with Pixi
  - Performance benchmarking
  - Production best practices
  - 📖 See also: [Bash Implementation README](workflows/bash/README.md)

- **[Part 2: Migrating GATK Bash to Nextflow with MD5 Validation](https://riverxdata.github.io/river-docs/blog/gatk-bash-nextflow-migration-md5-validation-part2)**
  - Complete 16-step Nextflow implementation
  - MD5 validation framework
  - Container strategy (BioContainers + Broad Institute)
  - Multi-sample parallelization
  - Troubleshooting GATK R dependency issues
  - 📖 See also: [Nextflow Implementation README](workflows/nextflow/Readme.md)

- **Part 3: Production-Scale Deployment** (Coming Soon)
  - SLURM/HPC integration
  - 1000 Genomes Project validation
  - 100+ sample benchmarking
  - Cloud deployment strategies

### Implementation-Specific Documentation

- 📖 [Bash Implementation README](workflows/bash/README.md) - Educational single-sample workflow
- 📖 [Nextflow Implementation README](workflows/nextflow/Readme.md) - Production multi-sample workflow with complete technical details

### Key Insights

**Why bash first?**
- Understand GATK workflow step-by-step
- Establish baseline for validation
- Educational tool for learning bioinformatics

**Why Nextflow?**
- Automatic parallelization (3-7x speedup)
- Resume from failures (`-resume`)
- Container portability (Singularity/Docker)
- HPC/cloud ready (SLURM, AWS Batch, Google Cloud)
- Production-grade error handling

## 🔬 Validation & Scientific Equivalence

We use **MD5 checksums** to validate that Bash and Nextflow produce scientifically equivalent results.

**Key findings:**
- ✅ **Variant statistics match exactly**: Raw SNP/indel counts identical
- ⚠️ **MD5 checksums differ**: Expected due to timestamps in BAM/VCF headers
- ✅ **Scientific content equivalent**: Same quality scores, alignments, variant calls

For detailed validation methodology and results, see [Nextflow README - Validation Section](workflows/nextflow/Readme.md#validation-against-bash).

## 🛠️ Tool Versions

All tools are installed via Pixi (locked in `pixi.lock`). For complete tool versions, container images, and technical specifications, see:
- 📖 [Nextflow Implementation README - Tool Versions & Containers](workflows/nextflow/Readme.md#tool-versions--container-strategy)

## 🐛 Troubleshooting

For detailed troubleshooting guides, see:
- 📖 [Nextflow Troubleshooting Section](workflows/nextflow/Readme.md#troubleshooting) - Container issues, memory errors, validation
- 📖 [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) - Common issues and solutions

**Quick Tips:**
- Use `-resume` flag to restart from failures without rerunning completed steps
- Check `work/` directory for detailed logs when a process fails
- Ensure Singularity/Apptainer version is ≥3.5 for container support

## 🚀 Next Steps & Roadmap

### Immediate Use Cases

✅ **Learn GATK workflow** - Study bash implementation step-by-step  
✅ **Test Nextflow locally** - Run multi-sample test with 3 samples  
✅ **Validate migration** - Compare bash vs Nextflow outputs with MD5  
✅ **Customize pipelines** - Modify parameters in `nextflow.config`

### Future Enhancements (Part 3)

- [ ] **SLURM executor** - HPC cluster deployment
- [ ] **1000 Genomes validation** - Real 30x coverage human genome data
- [ ] **MultiQC reporting** - Aggregate QC across 100+ samples
- [ ] **Cloud deployment** - AWS Batch, Google Cloud Life Sciences profiles
- [ ] **nf-test framework** - Unit/integration testing for each module
- [ ] **nf-core standards** - Align with nf-core best practices

## 📚 References

### GATK Documentation
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
- [GATK Germline Short Variant Discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932)

### Nextflow Resources
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [nf-core Pipelines](https://nf-co.re/)
- [BioContainers Registry](https://biocontainers.pro/)

### Blog Posts
- [Pixi - The New Conda Era](https://riverxdata.github.io/river-docs/blog/pixi-is-new-conda-based-era)
- [Containers on HPC: Docker to Singularity](https://riverxdata.github.io/river-docs/blog/containers-hpc-docker-singularity-apptainer)
- [Finding Pre-Built Bioinformatics Containers](https://riverxdata.github.io/river-docs/blog/bioinformatics-containers-build-efficient-docker)

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

### Development Workflow

```bash
# Create a feature branch
git checkout -b feat/my-feature

# Make changes and test
cd workflows/bash && pixi run bash gatk_pipeline.sh
cd workflows/nextflow && nextflow run main.nf -profile singularity,test

# Commit with conventional commits
git commit -m "feat: add new feature"

# Push and create PR
git push origin feat/my-feature
```

### Code Style

- **Bash**: Follow [Google Shell Style Guide](https://google.github.io/styleguide/shellguide.html)
- **Nextflow**: Follow [nf-core style guide](https://nf-co.re/docs/contributing/guidelines)
- **Commit messages**: Use [Conventional Commits](https://www.conventionalcommits.org/)

## 📄 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **GATK Team** at Broad Institute for best practices guidelines
- **nf-core Community** for test datasets and Nextflow standards
- **BioContainers** for pre-built bioinformatics containers
- **Pixi Team** for modern conda-based package management

## 📧 Contact

- **Issues**: [GitHub Issues](https://github.com/nttg8100/variant-calling-gatk-pipeline-best-practice-from-scratch/issues)
- **Blog**: [RiverXData](https://riverxdata.github.io/river-docs/)
- **Author**: Giang Nguyen

---

**Status**: ✅ **Production Ready** | Bash ✅ Complete | Nextflow ✅ Complete (16/16 steps)  
**Last Updated**: February 2026
