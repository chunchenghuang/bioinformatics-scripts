#!bin/bash
source /opt/miniconda2/etc/profile.d/conda.sh
conda activate gatk

dir="."
ref_dir="/home/public/hg38_1000G"
ref_file=${ref_dir}/"Homo_sapiens_assembly38.fasta"
#known snp site files
omni=${ref_dir}/"1000G_omni2.5.hg38.vcf.gz"
dbsnp=${ref_dir}/"dbsnp_146.hg38.vcf.gz"
hapmap=${ref_dir}/"hapmap_3.3.hg38.vcf.gz"
indels=${ref_dir}/"Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
confidence=${ref_dir}/"1000G_phase1.snps.high_confidence.hg38.vcf.gz"
germline=${ref_dir}/"af-only-gnomad.hg38.vcf.gz"
#target region
regions=${dir}/"gene_list.bed"
#gene name
gene_name="PGK1"
#gene range
gene_range="chrX:78104340-78125830"


#reference index from alignment
if [ ! -e "$ref_file.sa" ]; then
    bwa index $ref_file
fi

if [ ! -e "$ref_file.fai" ]; then
    samtools faidx $ref_file
fi

ref_file_name=${ref_file//.fa/}
if [ ! -e "$ref_file_name.dict" ]; then
    java -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar CreateSequenceDictionary \
        --REFERENCE $ref_file --OUTPUT $ref_file_name.dict
fi

for list_R1 in $(find $dir -name "*R1*.fastq.gz"|sort); do
    list_R1=${list_R1//.fastq.gz/}
    list_R2=${list_R1//_R1/_R2}
    bam=${list_R1//_R1*/}
    sample_name=$(echo "$list_R1"|awk -F "/" '{print $NF}'|cut -d "_" -f 1)
    group_id=$(zcat $list_R1.fastq.gz | head -n 1 | tr -d "@" | cut -d ":" -f 3,4 | tr ":" ".")
    library=$(zcat $list_R1.fastq.gz | head -n 1 | awk -F ":" '{print $NF}')
    header=$(echo "@RG\tID:"$group_id"\tSM:"${sample_name}"\tPL:ILLUMINA\tLB:"$library)
    #if samples are grouped, input in tID
    if [ ! -e "${sample_name}.sam" ] && [ ! -e "${sample_name}.bam" ]; then
        bwa mem -M -t 20 -R $header $ref_file ${list_R1}.fastq.gz ${list_R2}.fastq.gz |\
        samtools view -@ 20 -S -b - > ${bam}.bam
    fi
done

for sample_name in $(find $dir -name "*.bam"|awk -F "/" '{print $NF}'|cut -d "_" -f 1 | sort | uniq); do
sample_name=${sample_name//L*.bam/}
    if [ ! -e "$sample_name.bam" ]; then
        if [ $(ls $sample_name*.bam|wc -l) -gt 1 ]; then
            ls $sample_name*.bam > $sample_name.bamlist
            bamtools merge -list $sample_name.bamlist -out $sample_name.bam
        else
            cp $sample_name*.bam $sample_name.bam
        fi
        bamtools sort -mem 80000 -in $sample_name.bam -out $sample_name.bam
        bamtools index -in ${sample_name}.bam
        rm $sample_name*_L*.bam
    fi
done

ref_file=${ref_dir}/"PGK1.fa"
#reference index from alignment
if [ ! -e "$ref_file.sa" ]; then
    bwa index $ref_file
fi

if [ ! -e "$ref_file.fai" ]; then
    samtools faidx $ref_file
fi

ref_file_name=${ref_file//.fa/}
if [ ! -e "$ref_file_name.dict" ]; then
    java -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar CreateSequenceDictionary \
        --REFERENCE $ref_file --OUTPUT $ref_file_name.dict
fi

for list_R1 in $(find $dir -name "*R1*.fastq.gz"|sort); do
    list_R1=${list_R1//.fastq.gz/}
    list_R2=${list_R1//_R1/_R2}
    sample_name=$(echo "$list_R1"|awk -F "/" '{print $NF}'|cut -d "_" -f 1)
    group_id=$(cat ${list_R1}_${gene_name}.fastq | head -n 1 | tr -d "@" | cut -d ":" -f 3,4 | tr ":" ".")
    library=$(cat ${list_R1}_${gene_name}.fastq | head -n 1 | awk -F ":" '{print $NF}')
    header=$(echo "@RG\tID:"$group_id"\tSM:"${sample_name}"\tPL:ILLUMINA\tLB:"$library)
    samtools view ${sample_name}.bam $gene_range | awk '{print $1}' > ${sample_name}_${gene_name}.list
    seqtk subseq $list_R1.fastq.gz ${sample_name}_${gene_name}.list 1> ${list_R1}_${gene_name}.fastq
    seqtk subseq $list_R2.fastq.gz ${sample_name}_${gene_name}.list 1> ${list_R2}_${gene_name}.fastq
    bwa mem -M -t 20 -R $header $ref_file ${list_R1}_${gene_name}.fastq ${list_R2}_${gene_name}.fastq |\
        samtools view -@ 20 -S -b - > ${sample_name}_${gene_name}_realigned.bam
    bamtools sort -in ${sample_name}_${gene_name}_realigned.bam -out ${sample_name}_${gene_name}_realigned.bam
    bamtools index -in ${sample_name}_${gene_name}_realigned.bam
done

##modification

for bam in $(find $dir -name "*_realigned.bam"|sort|cut -d "_" -f1| sed 's/.bam//' | uniq); do
    samtools view -h $bam.bam | grep "^@" > header.tmp
    samtools view -h ${bam}_${gene_name}_realigned.bam | grep -v "^@" | awk -vOFS="\t" '{if($3=="PGK1"){$3="chrX";$4=$4+78104340-1}{print $0}}'| awk -vOFS="\t" '{if($7=="="){$8=$8+78104340-1}{print $0}}' > sam.tmp
    cat header.tmp sam.tmp | samtools view -@ 20 -S -b - > ${bam}_${gene_name}_realigned_mod.bam
    bamtools sort -in ${bam}_${gene_name}_realigned_mod.bam -out ${bam}_${gene_name}_realigned_mod.bam
    bamtools index -in ${bam}_${gene_name}_realigned_mod.bam
done

ref_file=${ref_dir}/"Homo_sapiens_assembly38.fasta"
regions=${dir}/"CYP21A2.bed"

#call variants
for sample_name in $(find $dir -name "*realigned_mod.bam"); do
    sample_name=${sample_name//.bam/}
    sample_name=${sample_name//_markduplicates*/}
    #GATK mark duplicate and bamtools index
    if [ ! -e "${sample_name}_markduplicates.bam" ]; then
        java -Xmx100g -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar MarkDuplicates \
            -I $sample_name.bam \
            -O ${sample_name}_markduplicates.bam \
            -M ${sample_name}_markduplicates.metrics
        bamtools index -in ${sample_name}_markduplicates.bam
    fi
    #GATK base recalibrator
    if [ ! -e "${sample_name}.realn.recal" ]; then
        java -Xmx100g -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar BaseRecalibrator \
            -R $ref_file \
            -I ${sample_name}_markduplicates.bam \
            -L $regions \
            --known-sites $omni \
            --known-sites $confidence \
            --known-sites $dbsnp \
            --known-sites $hapmap \
            --known-sites $indels \
            -O ${sample_name}.realn.recal
    fi
    #apply recalibration
    if [ ! -e "${sample_name}_markduplicates_bqsr.bam" ]; then
        java -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar ApplyBQSR \
            -R $ref_file \
            -I ${sample_name}_markduplicates.bam \
            -L $regions \
            -bqsr ${sample_name}.realn.recal \
            -O ${sample_name}_markduplicates_bqsr.bam
    fi
done

#calling germline mutants for single unmatched sample
for samples in $(find $dir -name "*_markduplicates_bqsr.bam");do
    samples=${samples//.bam/}
    #haplotype caller
    if [ ! -e "${samples}_raw.vcf" ]; then
        java -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar HaplotypeCaller \
            --native-pair-hmm-threads 10 \
            --input $samples.bam \
            -L $regions \
            -R $ref_file \
            --output ${samples}_raw.vcf
    fi
done

conda deactivate

#HAPCUT2 for haplotyping

for bam in $(find $dir -name "*realigned_mod.bam");do
    bam=${bam//.bam/}
    extractHAIRS --mmq 10 --bam $bam.bam --VCF ${bam}_markduplicates_bqsr_raw.vcf --out ${bam}_fragment_file
    HAPCUT2 --fragments ${bam}_fragment_file --VCF ${bam}_markduplicates_bqsr_raw.vcf --output ${bam}_haplotype_output_file
done