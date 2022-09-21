#!/bin/bash
dir="."
Quality="10"
for bam_file in $(find $dir -name "*.bam"); do
sample_name=${bam_file//.bam/}
samtools depth -aa -d 0 -Q $Quality $sample_name.bam > ${sample_name}_${Quality}_DP.txt
awk '{print $1"\t"$2-1"\t"$2"\t"$3}' ${sample_name}_${Quality}_DP.txt > ${sample_name}_${Quality}_DP.bed 
#bedtools intersect -a $sample_name_$Quality_DP.bed -b IDT_exome_hg19.bed > M158-1_S1_Q30_intersect.bed
#awk '{print $1"\t"$3"\t"$4}' M158-1_S1_Q30_intersect.bed > M158-1_S1_Q30_intersect_single.bed
#awk '{if($4 < 20){print $0}}' M158-1_S1_Q30_intersect.bed > M158-1_S1_Q30_intersect_coverage_20.bed
awk '{print $4}' ${sample_name}_${Quality}_DP.bed > ${sample_name}_${Quality}_DP.tmp
done
first_sample=$(ls *_DP.bed |head -n 1)
awk '{print $1"\t"$2"\t"$3}' $dir/$first_sample > table_row.header
sample_name=$(ls *_DP.tmp|cut -d "_" -f 1| tr "\n" "\t")
echo -e "Region\tStart\tEnd\t$sample_name" > table_column.header
sample_list=$(ls *_DP.tmp|tr "\n" " ")
paste table_row.header $sample_list > cov_table.tmp
cat table_column.header cov_table.tmp > cov_table.txt
rm *.tmp
rm *.header