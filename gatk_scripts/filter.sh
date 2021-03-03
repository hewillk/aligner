#!/usr/bin/bash

reference=/mnt/hewill/ref/hs37d5.fa
gatk=~/software/gatk-4.1.8.1/gatk
sample=$1

(time $gatk VariantFiltration \
 -R $reference \
 -V output/$sample.vcf.gz \
 -O output/$sample.filtered.vcf.gz \
 --filter-expression "QD < 2.0" \
 --filter-name "QD20" \
 --filter-expression "FS > 60.0" \
 --filter-name "FS60" \
 && echo "** ${sample} vcf filter done **") \
> log/${sample}.filter.log 2>&1