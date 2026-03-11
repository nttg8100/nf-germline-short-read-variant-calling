# cloning repo
git clone https://github.com/nttg8100/nf-germline-short-read-variant-calling.git -b 0.4.0

# Freebayes recommendation: 
# 1. fastp: quality control
# 2. bwa-mem2: read alignment  
# 3. sambamba: PCR remove duplicate
# 4. freebayes: Variant calling
nextflow run nf-germline-short-read-variant-calling \
    -profile docker -c customized_labels.config \
    --genome GRCh38 \
    --input input.csv \
    --variant_caller "freebayes" \
    --preprocessor "sambamba" \
    --output variant_calling


# Deepvariant recommendation: 
# 1. fastp: quality control
# 2. bwa-mem2: read alignment  
# 3. deepvariant: variant calling 
nextflow run nf-germline-short-read-variant-calling \
    -profile docker -c customized_labels.config \
    --genome GRCh38 \
    --input input.csv \
    --variant_caller "freebayes" \
    --preprocessor "sambamba" \
    --output variant_calling
