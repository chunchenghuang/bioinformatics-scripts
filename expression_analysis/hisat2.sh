#!/bin/bash
#environment setting
dir="."
ref_dir="./grcm38_snp_tran"
index=$ref_dir/"genome_snp_tran"
#cpu for one pair-end reads
#you need to consider maximum cpu that could be used parallely
cpu="25"

#loops
for list_R1 in $(find $dir -name "*R1_inserted.fastq"|sort); do
number=$(ls $list_R1 | wc -l)
list_R2=${list_R1//_R1/_R2};
output=${list_R1//_R1_inserted.fastq/};
##testing
##echo ${list_R1};
##echo ${list_R2};
##echo ${output};
##echo ${index_dir}zebrafinch_index 
##echo ${output}.sam
##echo ${output}.txt
hisat2 -p $cpu -x ${index} -1 ${list_R1} -2 ${list_R2} --summary-file ${output}_mice.txt --new-summary |\
samtools sort -@ 25 -O BAM > ${output}_mice.bam
#[ $(jobs|wc -l) -ge $(expr $(nproc) / $cpu) ] && wait
done
#wait

# for sample_name in $(find $dir -name "*.bam"|cut -d "_" -f 1 | sort | uniq); do
# sample_name=${sample_name//.bam/}
#     if [ ! -e "$sample_name.bam" ]; then
#         if [ $(ls $sample_name*.bam|wc -l) -gt 1 ]; then
#             ls $sample_name*.bam > $sample_name.bamlist
#             bamtools merge -list $sample_name.bamlist -out $sample_name.bam
#             rm $sample_name.bamlist
#         else
#             cp $sample_name*.bam $sample_name.bam
#         fi
#         bamtools sort -in $sample_name.bam -out $sample_name.bam
#         bamtools index -in $sample_name.bam
#         rm $sample_name*_L*.bam
#     fi
# done

# java -jar /opt/AlignerBoost/AlignerBoost.jar run filterPE -in S37-1.bam -out S37-1_AB.bam

# bamtools sort -mem 51200 -in S37-1_AB.bam -out S37-1_AB.bam

# bamtools index -in S37-1_AB.bam