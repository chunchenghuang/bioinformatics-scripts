awk '{if($3<200){print $0}}' CHD2_S1_10_DP.txt > CHD2_S1_10_DP_less_200.txt
fake_last=$(tail -1 CHD2_S1_10_DP_less_200.txt | awk '{print $1"\t"$2+2"\t"$3}')
echo "$fake_last" | cat CHD2_S1_10_DP_less_200.txt - > CHD2_S1_10_DP_less_200_new.tmp
awk 'NR == 1 {start = prev = $2; chr = $1; next}
{if ($2!=prev+1) {print chr"\t"start"\t"(prev==start ? start:prev);start=$2}
{prev = $2 ; chr = $1}}' CHD2_S1_10_DP_less_200_new.tmp > region.tmp
while IFS=$'\t' read -r -a line;
do
max=${line[2]}
min=${line[1]}
region=$(expr ${line[2]} - ${line[1]} + 1)
chr=${line[0]}
grep -m 1 -P "$chr\t$min" -A $(($region-1)) CHD2_S1_10_DP_less_200.txt | \
awk -v min="$min" -v max="$max" -v chr="$chr" -v region=$region \
'{sum+=$3}END{print chr"\t"min-1"\t"max"\t",sum/region}' 
done < region.tmp > CHD2_S1_10_DP_less_200_final.txt
