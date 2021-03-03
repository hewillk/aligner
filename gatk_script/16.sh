#!/usr/bin/bash

reference=/mnt/hewill/ref/hs37d5.fa
gatk=~/software/gatk-4.1.8.1/gatk
sample=$1

chrom=(
16:1-25000000
16:25000001-50000000
16:50000001-90354753
)
for i in ${chrom[@]}; do
    (time $gatk HaplotypeCaller \
     -R $reference \
     -I ${sample}.16.bam \
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