#!/bin/bash
while IFS=$'\t' read -r -a line;
do
 seq ${line[1]} ${line[2]} > POS.tmp
 repeat=$(expr ${line[2]} - ${line[1]} + 1)
 echo ${line[0]}| awk -v n="$repeat" '{for(i=0;i<n;i++)print}' > chr.tmp
 samtools faidx Homo_sapiens_assembly38.fasta ${line[0]}:${line[1]}-${line[2]} | sed 1d | tr -d "\n" | \
 while read -n1 i;do echo "$i"; done > ref.tmp
 paste chr.tmp POS.tmp ref.tmp
#echo $region $min $max
done < IDT_INHDIS_HG38.bed > chr_pos_ref.tmp
number=$(wc -l chr_pos_ref.tmp | awk '{print $1}')
echo '.' | awk -v n="$number" '{for(i=0;i<n;i++)print}' > id.tmp
grep "#" Pompe.vcf > Pompe.header
sed '/#/d' Pompe.vcf > Pompe.table
awk '{print $5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' Pompe.table | head -n 1 | \
awk -v n="$number" '{for(i=0;i<n;i++)print}' | \
awk '{{$1="G"}print}' | tr " " "\t" > rep_G.table
paste chr_pos_ref.tmp id.tmp rep_G.table > full_G.tmp
awk '{print $1"\t"$2"\t"$4"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' full_G.tmp > reorder_full_G.tmp
awk '{if($4!=$5){print $0}}' reorder_full_G.tmp > reorder_noGG_full.tmp
cat Pompe.header reorder_noGG_full.tmp > Pompe_fake_G.vcf
