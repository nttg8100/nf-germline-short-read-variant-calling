# raw fastq
samtools view -b HG002_merged.bam chr22:17000000-20000000 > subset.bam
samtools sort -n -o subset_name.bam subset.bam
samtools fastq \
  -1 HG002_subset_R1.fastq.gz \
  -2 HG002_subset_R2.fastq.gz \
  -0 /dev/null \
  -s /dev/null \
  -n subset_name.bam

# subset genome
# download only chr22
aws s3 ls --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/Chromosomes/chr22.fasta .
# index
samtools faidx chr22.fasta
# get fasta dict
picard CreateSequenceDictionary \
    R=chr22.fasta \
    O=chr22.dict