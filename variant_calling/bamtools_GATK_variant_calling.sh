#!/bin/bash
dir="."
for sample_name in $(find $dir -name "*.bam"|cut -d "_" -f 1,2|uniq); do
ls $sample_name*bam > $sample_name.bamlist
bamtools merge -list $sample_name.bamlist -out $sample_name.bam
bamtools sort -in $sample_name.bam -out $sample_name.bam
bamtools index -in $sample_name.bam
java -jar /opt/app/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../RHD.fa -I $sample_name.bam -nt 10 -o ${sample_name}_GATK.vcf -glm BOTH -dt NONE &
#testing
# echo ${sample_name}.bam
# echo $sample_name.bamlist
done
rm "$dir"/*.bamlist
#rm "$dir"/"$sample_name"_L*.bam

