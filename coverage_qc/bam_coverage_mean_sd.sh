#!/bin/bash
#the maximum cpu number this script uses = (maximum number of processor -2)
dir="."
Quality="10"
bed_file="IDT_exome_hg19.bed"
for chr in $(cut -f 1 $bed_file | sort -k1,1 -V -s | uniq); do
grep -w $chr $bed_file > IDT_exome_hg19_${chr}.bed
done
for bam_file in $(find $dir -name "*.bam"); do
sample_name=${bam_file//.bam/}
samtools depth -aa -b $bed_file -d 0 -Q $Quality $sample_name.bam > ${sample_name}_${Quality}_DP.txt
for chr in $(cut -f 1 ${sample_name}_${Quality}_DP.txt | sort -k1,1 -V -s | uniq); do
grep -w $chr ${sample_name}_${Quality}_DP.txt > ${sample_name}_${Quality}_DP_${chr}.tmp
while IFS=$'\t' read -r -a line;
do
max=${line[2]}
min=$(expr ${line[1]} + 1)
region=$(expr ${line[2]} - ${line[1]})
chr=${line[0]}
grep -m 1 -P "$chr\t$min" -A $(($region-1)) ${sample_name}_${Quality}_DP_${chr}.tmp | \
awk -v min="$min" -v max="$max" -v chr="$chr" -v region=$region \
'{sum+=$3;sumsq+=$3*$3}END{print chr"\t"min-1"\t"max"\t",sum/region"\t"sqrt((sumsq-(sum^2/region))/(region-1))}' &
[ $(jobs|wc -l) -ge $(expr $(nproc) - 2) ] && wait
done < IDT_exome_hg19_${chr}.bed > ${sample_name}_${Quality}_DP_${chr}_mean_sd.tmp
wait
done
cat ${sample_name}_${Quality}_DP_chr*_mean_sd.tmp | sort -k1,1V -k2,2n -k3,3n > ${sample_name}_${Quality}_DP_mean_sd.bed
awk '{print $4"\t"$5}' ${sample_name}_${Quality}_DP_mean_sd.bed > ${sample_name}_${Quality}_DP_mean_sd.tmp
done
rm "$dir"/*_DP_chr*
first_sample=$(ls *_DP_mean_sd.bed|sort|head -n 1)
awk '{print $1"\t"$2"\t"$3}' $dir/$first_sample > table_mean_row.header
sample_name_mean_sd=$(ls *_DP_mean_sd.tmp|sort|awk -F "_" '{print $1"_mean""\t"$1"_sd"}'|tr "\n" "\t")
echo -e "Region\tStart\tEnd\t$sample_name_mean_sd" > table_mean_column.header
sample_list=$(ls *_DP_mean_sd.tmp|sort|tr "\n" " ")
paste table_mean_row.header $sample_list > cov_mean_sd_table.tmp
cat table_mean_column.header cov_mean_table.tmp > cov_mean_sd_table.txt
rm "$dir"/*.tmp
rm "$dir"/*.header