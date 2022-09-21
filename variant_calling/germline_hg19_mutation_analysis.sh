#!bin/bash

#for haplotype caller you need conda gatk python paqcakges installed 
#environment
source /opt/miniconda2/etc/profile.d/conda.sh
conda activate gatk
dir="."
genome="hg19"

#hg19
ref_dir="/home/public/hg19_1000G"
ref_file=${ref_dir}/"ucsc.hg19.fasta"
#known snp site files
omni=${ref_dir}/"1000G_omni2.5.hg19.sites.vcf.gz"
dbsnp=${ref_dir}/"dbsnp_138.hg19.vcf.gz"
hapmap=${ref_dir}/"hapmap_3.3.hg19.sites.vcf.gz"
indels=${ref_dir}/"Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
confidence=${ref_dir}/"1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz"
germline=${ref_dir}/"af-only-gnomad.raw.sites.vcf"
#target region
regions=${dir}/"xgen-exome-research-panel-v2-targets-hg19_20.bed"


#operate under the target directory
cd $dir

for f in $(find $dir -name "*.bed"); do
    if [ ! -e "$f" ]; then
        echo "You need a bed file to run this script " >&2
        exit 1
    else
        break
    fi
done

for f in $(find $dir -name "*R1*.fastq.gz"|cut -d "_" -f 1); do
    f=${f//$dir\//}
    if [ -e "${f}_final.hg19_multianno.txt.intervar" ]; then
        echo "WARNING:Sample ${f}'s final output already exists!
Remove ${f}'s fastq.gz file if you do not wish to re-do the intermediate files."
    fi
done

#reference index from alignment
if [ ! -e "$ref_file.sa" ]; then
    bwa index $ref_file
fi

if [ ! -e "$ref_file.fai" ]; then
    samtools faidx $ref_file
fi

ref_file_name=${ref_file//.fasta/}
if [ ! -e "$ref_file_name.dict" ]; then
    java -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar CreateSequenceDictionary \
        --REFERENCE $ref_file --OUTPUT $ref_file_name.dict
fi

#trimmomatic
for list_R1 in $(find $dir -name "*R1_*.fastq.gz"|sort); do
    list_R1=${list_R1//_paired.fastq.gz/}
    list_R1=${list_R1//_unpaired.fastq.gz/}
    list_R1=${list_R1//.fastq.gz/}
    list_R2=${list_R1//_R1/_R2}
    if [ ! -e "${list_R1}_paired.fastq.gz" ]; then
        java -jar /opt/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 \
            ${list_R1}.fastq.gz ${list_R2}.fastq.gz \
            ${list_R1}_paired.fastq.gz ${list_R1}_unpaired.fastq.gz \
            ${list_R2}_paired.fastq.gz ${list_R2}_unpaired.fastq.gz \
            ILLUMINACLIP:/opt/Trimmomatic-0.38/adapters/TruSeq.fa:2:30:10:1:true \
            LEADING:30 TRAILING:30 MINLEN:100 SLIDINGWINDOW:4:20
    rm "$dir"/*_unpaired.fastq.gz
    fi
done

#bwa alignment, merge and sort bam file
for list_R1 in $(find $dir -name "*R1*.fastq.gz"|sort); do
    list_R1=${list_R1//.fastq.gz/}
    list_R2=${list_R1//_R1/_R2}
    bam=${list_R1//_R*/}
    sample_name=$(echo "$list_R1"|awk -F "/" '{print $NF}'|cut -d "_" -f 1)
    group_id=$(zcat $list_R1.fastq.gz | head -n 1 | tr -d "@" | cut -d ":" -f 3,4 | tr ":" ".")
    library=$(zcat $list_R1.fastq.gz | head -n 1 | awk -F ":" '{print $NF}')
    header=$(echo "@RG\tID:"$group_id"\tSM:"${sample_name}"\tPL:ILLUMINA\tLB:"$library)
    #if samples are grouped, input in tID
    if [ ! -e "${sample_name}.sam" ] && [ ! -e "${sample_name}.bam" ]; then
        bwa mem -M -t 20 -R $header $ref_file ${list_R1}.fastq.gz ${list_R2}.fastq.gz |\
        samtools view -@ 20 -S -b - > $bam.bam
    fi
done

#merge, sort bam files
for sample_name in $(find $dir -name "*.bam"|awk -F "/" '{print $NF}'|cut -d "_" -f 1| sort | uniq); do
    sample_name=${sample_name//.bam/}
    if [ ! -e "$sample_name.bam" ]; then
        if [ $(ls $sample_name*.bam|wc -l) -gt 1 ]; then
            ls $sample_name*.bam > $sample_name.bamlist
            bamtools merge -list $sample_name.bamlist -out $sample_name.bam
            rm $sample_name.bamlist
        else
            cp $sample_name*.bam $sample_name.bam
        fi
        bamtools sort -in $sample_name.bam -out $sample_name.bam
        bamtools index -in $sample_name.bam
        rm $sample_name*_L*.bam
        rm *.sam
    fi
done


for sample_name in $(find $dir -name "*.bam"); do
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
    #CNNScore
    if [ ! -e "${samples}_raw_CNNScore.vcf" ]; then
        java -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar CNNScoreVariants \
        -I $samples.bam \
        -V ${samples}_raw.vcf \
        -R $ref_file \
        -O ${samples}_raw_CNNScore.vcf \
        -tensor-type read_tensor
    fi
    #FilterVariantTranches
    if [ ! -e "${samples}_filtered.vcf" ]; then
        java -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar FilterVariantTranches \
        -V ${samples}_raw_CNNScore.vcf \
        --resource $hapmap \
        --resource $indels \
        --info-key CNN_2D \
        --snp-tranche 99.9 --snp-tranche 99.95 \
        --indel-tranche 99.0 --indel-tranche 99.4 \
        -O ${samples}_filtered.vcf
    fi
done

#filter PASS calls
for vcfs in $(find $dir -name "*_filtered.vcf"); do
    vcfs=${vcfs//.vcf/}
    if [ ! -e "${vcfs}_pass.vcf" ]; then
        java -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar SelectVariants \
            -V $vcfs.vcf \
            --exclude-filtered \
            -O ${vcfs}_pass.vcf
    fi
done

#not yet download full dataset for hg19
# #ANNOVAR and Intervar variant annotation
# for pass_vcfs in $(find $dir -name "*_pass.vcf");do
#     pass_vcfs=${pass_vcfs//.vcf/}
#     out=${pass_vcfs//_markduplicates*/}
#     out=${out//$dir\//}
#     #ANNOVAR format conversion
#     if [ ! -e "$pass_vcfs.avinput" ]; then
#         /home/hao/Downloads/annovar/convert2annovar.pl \
#             -format vcf4 $pass_vcfs.vcf \
#             -outfile $pass_vcfs.avinput \
#             -allsample -withfreq -include
#     fi
#     #ANNOVAR
#     if [ ! -e "${out}_final.hg19_multianno.txt" ]; then
#         /home/hao/Downloads/annovar/table_annovar.pl $pass_vcfs.avinput \
#             /home/hao/Downloads/annovar/humandb \
#             --buildver hg19 \
#             --remove \
#             --outfile ${out}_final \
#             --protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,avsnp150,dbnsfp33a,dbnsfp35c,clinvar_20190305,exac03,cosmic70,gnomad_genome,gnomad30_genome,dbscsnv11,dbnsfp31a_interpro,rmsk,ensGene,knownGene \
#             --operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,r,g,g \
#             --nastring . \
#             --otherinfo
#     fi
#     #Intervar
#     if [ ! -e "${out}_final.hg19_multianno.txt.intervar" ]; then
#         /home/hao/Downloads/InterVar-master/Intervar.py \
#             --table_annovar=/home/hao/Downloads/annovar/table_annovar.pl \
#             --convert2annovar=/home/hao/Downloads/annovar/convert2annovar.pl \
#             --annotate_variation=/home/hao/Downloads/annovar/annotate_variation.pl \
#             --database_locat=/home/hao/Downloads/annovar/humandb \
#             --database_intervar=/home/hao/Downloads/InterVar-master/intervardb \
#             --skip_annovar \
#             -b hg19 -i $pass_vcfs.avinput --input_type=AVinput -o ${out}_final
#     fi
#     #final column re-order of the annovar and intervar file
#     if [ ! -e "${out}_final_hg19_selected.intervar" ]; then
#         awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$192"\t"$197}' ${out}_final.hg19_multianno.txt |\
#         awk '($6~/,/){print $0}' > ${out}_multi_pheno.tmp
#         awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$192"\t"$197}' ${out}_final.hg19_multianno.txt |\
#         awk '($6!~/,/){print $0}' > ${out}_single_pheno.tmp
#         sed 's/:/\t/g;s/,/\t/g' ${out}_multi_pheno.tmp | awk 'NR%2' |\
#         awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10/($9+$10+$11)}' > ${out}_odd_ratio.tmp
#         sed 's/:/\t/g;s/,/\t/g' ${out}_multi_pheno.tmp | awk '!(NR%2)' |\
#         awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11/($9+$10+$11)}' > ${out}_even_ratio.tmp
#         sed 's/:/\t/g;s/,/\t/g' ${out}_single_pheno.tmp | sed 1d |\
#         awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$9/($8+$9)}' > ${out}_single_ratio.tmp
#         echo -e "$(head -1 ${out}_single_pheno.tmp|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}')\tAllele_ratio" |\
#         sed 's/\t/$/;s/\t/$/;s/\t/$/;s/\t/$/' > ${out}_ratio_all_wanted_header.tmp
#         cat ${out}_odd_ratio.tmp ${out}_even_ratio.tmp ${out}_single_ratio.tmp |\
#         sed 's/\t/$/;s/\t/$/;s/\t/$/;s/\t/$/' | sort -k1,1V > ${out}_ratio_all_wanted_noheader.tmp
#         awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$89"\t"$154}' ${out}_final.hg19_multianno.txt |\
#         sed 's/\t/$/;s/\t/$/;s/\t/$/;s/\t/$/' > ${out}_annovar_wanted.tmp
#         head -1 ${out}_annovar_wanted.tmp > ${out}_annovar_wanted_header.tmp
#         sed 1d ${out}_annovar_wanted.tmp | sort -k1,1V > ${out}_annovar_wanted_noheader.tmp
#         cut --complement -f 16,18,19,23,24,25,26,27,28,34 ${out}_final.hg19_multianno.txt.intervar | tr -d "#" |\
#         sed 's/\t/$/;s/\t/$/;s/\t/$/;s/\t/$/' > ${out}_intervar_wanted.tmp
#         head -1 ${out}_intervar_wanted.tmp > ${out}_intervar_wanted_header.tmp
#         sed 1d ${out}_intervar_wanted.tmp | sed 's/^/chr/' | sort -k1,1V > ${out}_intervar_wanted_noheader.tmp
#         join -t $'\t' -1 1 -2 1 ${out}_ratio_all_wanted_header.tmp ${out}_annovar_wanted_header.tmp > ${out}_semi_final_header.tmp
#         join -t $'\t' -1 1 -2 1 ${out}_semi_final_header.tmp ${out}_intervar_wanted_header.tmp > ${out}_final_header.tmp
#         join -t $'\t' -1 1 -2 1 ${out}_ratio_all_wanted_noheader.tmp ${out}_annovar_wanted_noheader.tmp >  ${out}_semi_final_noheader.tmp
#         join -t $'\t' -1 1 -2 1 ${out}_semi_final_noheader.tmp ${out}_intervar_wanted_noheader.tmp > ${out}_final_noheader.tmp
#         cat ${out}_final_header.tmp ${out}_final_noheader.tmp | tr "$" "\t" |\
#         sort -k1,1V -k2,2n -k3,3n -k4 -k5 > ${out}_final_hg19_selected.intervar
#         rm "$dir"/*.tmp
#     else
#         echo "All outputs exists!"
#     fi
# done
conda deactivate