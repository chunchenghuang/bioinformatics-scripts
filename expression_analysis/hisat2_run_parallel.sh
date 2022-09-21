#!/bin/bash
dir="/home/cch/zebrafinch_brain_pre_3/trimmomatic_paired/"
index_dir="/home/cch/zebrafinch_brain_pre/zebrafinch_HISAT2/"
output_dir="/home/cch/zebrafinch_brain_pre_3/HISAT2_trim_paired/"
for list_R1 in $(find $dir -name "*L001_R1_paired.fastq.gz"|sort); do
list_R2=${list_R1//_R1/_R2};
output=${list_R1//$dir/$output_dir};
output=${output//_R1_paired.fastq.gz/_paired_HISAT2};
##testing
##echo ${list_R1};
##echo ${list_R2};
##echo ${output};
##echo ${index_dir}zebrafinch_index 
##echo ${output}.sam
##echo ${output}.txt
hisat2 -p 12 --dta -x ${index_dir}zebrafinch_index -1 ${list_R1} -2 ${list_R2} -S ${output}.sam --summary-file ${output}.txt --new-summary
done&
for list_R1 in $(find $dir -name "*L002_R1_paired.fastq.gz"|sort); do
list_R2=${list_R1//_R1/_R2};
output=${list_R1//$dir/$output_dir};
output=${output//_R1_paired.fastq.gz/_paired_HISAT2};
##testing
##echo ${list_R1};
##echo ${list_R2};
##echo ${output};
##echo ${index_dir}zebrafinch_index 
##echo ${output}.sam
##echo ${output}.txt
hisat2 -p 12 --dta -x ${index_dir}zebrafinch_index -1 ${list_R1} -2 ${list_R2} -S ${output}.sam --summary-file ${output}.txt --new-summary
done&
for list_R1 in $(find $dir -name "*L002_2_R1_paired.fastq.gz"|sort); do
list_R2=${list_R1//_R1/_R2};
output=${list_R1//$dir/$output_dir};
output=${output//_R1_paired.fastq.gz/_paired_HISAT2};
##testing
##echo ${list_R1};
##echo ${list_R2};
##echo ${output};
##echo ${index_dir}zebrafinch_index 
##echo ${output}.sam
##echo ${output}.txt
hisat2 -p 12 --dta -x ${index_dir}zebrafinch_index -1 ${list_R1} -2 ${list_R2} -S ${output}.sam --summary-file ${output}.txt --new-summary
done