scratch/data/nf-germline-short-read-variant-calling
nextflow run main.nf \
    -profile docker -c /scratch/data/nf-germline-short-read-variant-calling/benchmark/pipelines/nf-germline-short-read-variant-calling/customized_labels.config \
    --genome GRCh38 \
    --input /scratch/data/nf-germline-short-read-variant-calling/benchmark/pipelines/nf-germline-short-read-variant-calling/input.csv \
    --variant_caller "freebayes" \
    --preprocessor "sambamba" \
    --output variant_calling
