#!bin/bash
dir="."
ref_dir="/home/public/hg38_1000G"
ref_file=${ref_dir}/"Homo_sapiens_assembly38.fasta"
bed=$dir/"PGK1.bed"
#gene name
gene_name="PGK1"

for bam in $(find $dir -name "*.bam"|sort|cut -d "_" -f1| sed 's/.bam//' | uniq); do
    samtools depth -aa -b $bed -d 0 -Q 10 ${bam}_${gene_name}_realigned_mod.bam > ${bam}_${gene_name}_10_DP.txt
    bgzip -c ${bam}_${gene_name}_realigned_mod_markduplicates_bqsr_raw.vcf > ${bam}_${gene_name}_realigned_mod_markduplicates_bqsr_raw.vcf.gz
    tabix -p vcf ${bam}_${gene_name}_realigned_mod_markduplicates_bqsr_raw.vcf.gz
    cat $ref_file | bcftools consensus ${bam}_${gene_name}_realigned_mod_markduplicates_bqsr_raw.vcf.gz > ${bam}_${gene_name}_consensus.fa
    python3 /home/cch/scripts/cnv_calculation.py ${bam}_${gene_name}_10_DP.txt ${bam}_${gene_name}_consensus.fa $bed
done