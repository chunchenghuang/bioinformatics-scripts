#!/bin/bash
#for line in $(cat filename); do
#    read input
#    echo $input;
#done
#or use while loop
cp all_ratio_table.txt region_all_ratio_table.txt
while IFS=$'\t' read -r -a line;
do
max=$(expr ${line[3]} - 25599028)
min=$(expr ${line[2]} - 25599028)
region=${line[0]}
awk -v min="$min" -v max="$max" -v region="$region" '{if($2 <= max && $2 >= min){$1 = region}{print $0}}' region_all_ratio_table.txt | tr -s " " "\t" > region_all_ratio_table.tmp
rm region_all_ratio_table.txt
mv region_all_ratio_table.tmp region_all_ratio_table.txt
#echo $region $min $max
done < RHD_exons.bed #input here