#!/bin/bash
dir="."
for sam in $(find $dir -name "*.sam"|sort); do
bam=${sam//.sam/.bam};
##testing
#echo ${sam}
#echo ${bam}
samtools view -S -b ${sam} > ${bam}
bamtools sort -in ${bam} -out ${bam}
bamtools index -in ${bam}
done
rm "$dir"/*.sam
##samtools sort somehow cannot be run parallelly at the same time