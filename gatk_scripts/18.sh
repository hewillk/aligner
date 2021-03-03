#!/usr/bin/bash

reference=/mnt/hewill/ref/hs37d5.fa
gatk=~/software/gatk-4.1.8.1/gatk
sample=$1

chrom=(
18:1-40000000
18:40000001-78077248
)
for i in ${chrom[@]}; do
    (time $gatk HaplotypeCaller \
     -R $reference \
     -I ${sample}.18.bam \
     -L $i \
     --read-filter ProperlyPairedReadFilter \
     --read-filter NotSupplementaryAlignmentReadFilter \
     --allow-non-unique-kmers-in-ref true \
     --min-dangling-branch-length 10000000 \
     --base-quality-score-threshold 6 \
     -O output/${sample}.${i}.vcf.gz \
     && echo "** ${sample}.${i}.vcf.gz done **") \
     > log/${sample}.${i}.log 2>&1 &
done && wait