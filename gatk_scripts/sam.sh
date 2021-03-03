#!/usr/bin/bash

samtools=~/software/samtools-1.10/samtools
sample=$1

chrom=(1 2 3 4 5 6 7 X 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
for i in ${chrom[@]}; do
    $samtools sort -@ 80 -l 1 --no-PG $sample.$i.sam > $sample.$i.bam
    $samtools index -@ 80 $sample.$i.bam
    rm $sample.$i.sam &
    bash $i.sh $sample &
done && wait
