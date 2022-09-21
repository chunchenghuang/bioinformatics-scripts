dir="."

for vcf in $(find $dir -name "*raw.vcf"|sort); do
    output=${vcf//_raw.vcf/_raw_mod.vcf}
    grep -v "^#" "$vcf" | awk -vOFS="\t" '{$1="chr6";$2=$2+32038415-1;print $0}'>vcf.temp
    grep "^#" "$vcf" > header.temp
    cat header.temp vcf.temp > $output
done

# for bam in $(find $dir -maxdepth 1 -name "*.bam"|sort|cut -d "_" -f1| sed 's/.bam//' | uniq); do
#     samtools view -h ${bam}.bam | grep "^@" > header.tmp
#     samtools view -h ${bam}_markduplicates.bam | grep -v "^@" | awk -vOFS="\t" '{if($3=="CYP21A2"){$3="chr6";$4=$4+32038415-1}{print $0}}'| awk -vOFS="\t" '{if($7=="="){$8=$8+32038415-1}{print $0}}' > sam.tmp
#     cat header.tmp sam.tmp | samtools view -@ 20 -S -b - > ${bam}_markduplicates_mod.bam
#     bamtools sort -in ${bam}_markduplicates_mod.bam -out ${bam}_markduplicates_mod.bam
#     bamtools index -in ${bam}_markduplicates_mod.bam
# done