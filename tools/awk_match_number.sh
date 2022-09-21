dir="."
bedfile1="../RHD_common_snp.bed"
bedfile2="../RHCE_common_snp.bed"
gene_starting_position="25599028"
cat $bedfile1 $bedfile2 > all.bed
awk -v a="$gene_starting_position" '{print "\t"$3-a"\t"}' all.bed > common_snp_position
#snp's bp position will be at your column three of your bed file
for vcf_files in $(find $dir -name "*GATK.vcf"|sort); do
clean_vcf=${vcf_files//.vcf/_clean.vcf}
sample_common_snp=${vcf_files//.vcf/_csnp.vcf}
grep -f common_snp_position $vcf_files > $sample_common_snp
grep -vf common_snp_position $vcf_files > $clean_vcf
done
