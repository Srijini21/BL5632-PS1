## Snakefile 

from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3=S3RemoteProvider()

SAMPLES=["SRR2589044"]

rule all:
    input:S3.remote(expand("s3://ps1-bl5632/results/final-variants/{SRR}_final_variants.vcf",SRR=SAMPLES))

rule bowtie:
    input:
        fq1=S3.remote(expand("s3://ps1-bl5632/results/trimmed/{SRR}_1.trim.fastq.gz",SRR=SAMPLES)),
        fq2=S3.remote(expand("s3://ps1-bl5632/results/trimmed/{SRR}_2.trim.fastq.gz",SRR=SAMPLES)),
    output:S3.remote(expand("s3://ps1-bl5632/results/sam/{SRR}.sam",SRR=SAMPLES))
    singularity:"docker://biocontainers/bowtie2:v2.4.1_cv1"
    shell:"bowtie2 -x genome/Ec606 --very-fast --threads 4 -1 {input.fq1} -2 {input.fq2} -S {output}"

rule bam:
    input:S3.remote(expand("s3://ps1-bl5632/results/sam/{SRR}.sam",SRR=SAMPLES))
    output:S3.remote(expand("s3://ps1-bl5632/results/bam/{SRR}.bam",SRR=SAMPLES))
    singularity:"docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:"samtools view -S -b {input} > {output}"

rule sort:
    input:S3.remote(expand("s3://ps1-bl5632/results/bam/{SRR}.bam",SRR=SAMPLES))
    output:S3.remote(expand("s3://ps1-bl5632/results/sorted-bam/{SRR}-sorted.bam",SRR=SAMPLES))
    singularity:"docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:"samtools sort -o {output} {input}"

rule index:
    input:S3.remote(expand("s3://ps1-bl5632/results/sorted-bam/{SRR}-sorted.bam",SRR=SAMPLES))
    output:S3.remote(expand("s3://ps1-bl5632/results/sorted-bam/{SRR}-sorted.bam.bai",SRR=SAMPLES))
    singularity:"docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:"samtools index {input}"

rule mpileup:
    input:
        genome=S3.remote(expand("s3://ps1-bl5632/data-genome/ecoli_rel606.fasta",SRR=SAMPLES)),
        bam=S3.remote(expand("s3://ps1-bl5632/results/sorted-bam/{SRR}-sorted.bam",SRR=SAMPLES)),
        bai=S3.remote(expand("s3://ps1-bl5632/results/sorted-bam/{SRR}-sorted.bam.bai",SRR=SAMPLES))
    output:S3.remote(expand("s3://ps1-bl5632/results/bcf/{SRR}.bcf",SRR=SAMPLES))
    singularity:"docker://biocontainers/bcftools:v1.9-1-deb_cv1"
    shell:"bcftools mpileup -O b -o {output} -f {input.genome} {input.bam}"

rule variant_calling:
   input:S3.remote(expand("s3://ps1-bl5632/results/bcf/{SRR}.bcf",SRR=SAMPLES))
   output:S3.remote(expand("s3://ps1-bl5632/results/vcf/{SRR}-variants.vcf",SRR=SAMPLES))
   singularity:"docker://biocontainers/bcftools:v1.9-1-deb_cv1"
   shell:"bcftools call --ploidy 1 -m -v {input} > {output}"

rule final_variants:
   input:S3.remote(expand("s3://ps1-bl5632/results/vcf/{SRR}-variants.vcf",SRR=SAMPLES))
   output:S3.remote(expand("s3://ps1-bl5632/results/final-variants/{SRR}_final_variants.vcf",SRR=SAMPLES))
   singularity:"docker://migbro/vcfutils:latest"
   shell:"vcfutils.pl varFilter {input} > {output}"
