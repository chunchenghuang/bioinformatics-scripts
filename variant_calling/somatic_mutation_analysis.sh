#!bin/bash

#environment
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
exac=${ref_dir}/"small_exac_common_3.hg38.vcf.gz"
#target region
regions=${dir}/"hg38_IDT_temp_110failed.bed"
#sample list and treatments
sample_list=${dir}/"sample_list.txt"

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
    if [ -e "${f}_final.hg38_multianno.txt.intervar" ]; then
        echo "WARNING:Sample ${f}'s final output already exists!
Remove ${f}'s fastq.gz file if you do not wish to re-do the intermediate files."
    fi
done

# if [ ! -e "$sample_list" ]; then
#     echo "must provide a tab-delimited sample_list.txt in the following format.
# format: <sample_name>   <normal/tumor>
# for example:
# sample_1    normal
# sample_2    tumor
# sample_3    tumor
# sample_4    normal" >&2
#     exit 1
# fi

#reference index from alignment
if [ ! -e "$ref_file.sa" ]  then
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
for list_R1 in $(find $dir -name "*R1*.fastq.gz"|sort); do
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

#bwa alignment
for list_R1 in $(find $dir -name "*R1*_paired.fastq.gz"|sort); do
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
        samtools view -@ 20 -S -b - > $bam.bam
    fi
done


#merge, sort bam files
for sample_name in $(find $dir -name "*.bam"|awk -F "/" '{print $NF}'|cut -d "_" -f 1 | sort | uniq); do
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
    sample_name=${sample_name//$dir\//}
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

#calling somatic mutants for single unmatched tumor sample
for tumor_samples in $(find $dir -name "*_markduplicates_bqsr.bam");do
tumor_samples=${tumor_samples//.bam/}
    #Mutect2
    if [ ! -e "${tumor_samples}_raw.vcf" ]; then
        java -Xmx100g -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar Mutect2 \
            --native-pair-hmm-threads 10 \
            -I $tumor_samples.bam \
            -tumor $tumor_samples \
            -L $regions \
            -R $ref_file \
            --germline-resource $germline \
            --af-of-alleles-not-in-resource 0.0000312 \
            --f1r2-tar-gz ${tumor_samples}_f1r2.tar.gz \
            --output ${tumor_samples}_raw.vcf
            #-pon pon.vcf.gz \
    fi
    #LearnOrientationModel
    if [ ! -e "${tumor_samples}_read-orientation-model.tar.gz" ]; then
        java -Xmx100g -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar LearnReadOrientationModel \
            -I ${tumor_samples}_f1r2.tar.gz \
            -O ${tumor_samples}_read-orientation-model.tar.gz
    fi
    #GetPileupSummaries
    if [ ! -e "${tumor_samples}_getpileupsummaries.table" ]; then
        java -Xmx100g -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar GetPileupSummaries \
            -I $tumor_samples.bam \
            -V $exac \
            -L $exac \
            -O ${tumor_samples}_getpileupsummaries.table
    fi
    #CalculateContamination
    if [ ! -e "${tumor_samples}_contamination.table" ]; then
        java -Xmx100g -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar CalculateContamination \
            -I ${tumor_samples}_getpileupsummaries.table \
            -tumor-segmentation ${tumor_samples}_segments.table \
            -O ${tumor_samples}_contamination.table
    fi
    #FilterMutectCalls
    if [ ! -e "${tumor_samples}_filtered.vcf" ]; then
        java -Xmx100g -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar FilterMutectCalls \
            -R $ref_file \
            --stats \
            -V ${tumor_samples}_raw.vcf \
            --tumor-segmentation ${tumor_samples}_segments.table \
            --contamination-table ${tumor_samples}_contamination.table \
            --ob-priors ${tumor_samples}_read-orientation-model.tar.gz \
            -O ${tumor_samples}_filtered.vcf
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

#ANNOVAR and Intervar variant annotation
for pass_vcfs in $(find $dir -name "*_pass.vcf");do
    pass_vcfs=${pass_vcfs//.vcf/}
    out=${pass_vcfs//_markduplicates*/}
    out=${out//$dir\//}
    #ANNOVAR format conversion
    if [ ! -e "$pass_vcfs.avinput" ]; then
        /home/hao/Downloads/annovar/convert2annovar.pl \
            -format vcf4 $pass_vcfs.vcf \
            -outfile $pass_vcfs.avinput \
            -allsample -withfreq -include
    fi
    #ANNOVAR
    if [ ! -e "${out}_final.hg38_multianno.txt" ]; then
        /home/hao/Downloads/annovar/table_annovar.pl $pass_vcfs.avinput \
            /home/hao/Downloads/annovar/humandb \
            --buildver hg38 \
            --remove \
            --outfile ${out}_final \
            --protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,avsnp150,dbnsfp33a,dbnsfp35c,clinvar_20190305,exac03,cosmic70,gnomad_genome,gnomad30_genome,dbscsnv11,dbnsfp31a_interpro,rmsk,ensGene,knownGene \
            --operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,r,g,g \
            --nastring . \
            --otherinfo
    fi
    #Intervar
    if [ ! -e "${out}_final.hg38_multianno.txt.intervar" ]; then
        /home/hao/Downloads/InterVar-master/Intervar.py \
            --table_annovar=/home/hao/Downloads/annovar/table_annovar.pl \
            --convert2annovar=/home/hao/Downloads/annovar/convert2annovar.pl \
            --annotate_variation=/home/hao/Downloads/annovar/annotate_variation.pl \
            --database_locat=/home/hao/Downloads/annovar/humandb \
            --database_intervar=/home/hao/Downloads/InterVar-master/intervardb \
            --skip_annovar \
            -b hg38 -i $pass_vcfs.avinput --input_type=AVinput -o ${out}_final
    fi
    #final column re-order of the annovar and intervar file
    if [ ! -e "${out}_final_hg38_selected.intervar" ]; then
        awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$192"\t"$197}' ${out}_final.hg38_multianno.txt |\
        awk '($6~/,/){print $0}' > ${out}_multi_pheno.tmp
        awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$192"\t"$197}' ${out}_final.hg38_multianno.txt |\
        awk '($6!~/,/){print $0}' > ${out}_single_pheno.tmp
        sed 's/:/\t/g;s/,/\t/g' ${out}_multi_pheno.tmp | awk 'NR%2' |\
        awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10/($9+$10+$11)}' > ${out}_odd_ratio.tmp
        sed 's/:/\t/g;s/,/\t/g' ${out}_multi_pheno.tmp | awk '!(NR%2)' |\
        awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$11/($9+$10+$11)}' > ${out}_even_ratio.tmp
        sed 's/:/\t/g;s/,/\t/g' ${out}_single_pheno.tmp | sed 1d |\
        awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$9/($8+$9)}' > ${out}_single_ratio.tmp
        echo -e "$(head -1 ${out}_single_pheno.tmp|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}')\tAllele_ratio" |\
        sed 's/\t/$/;s/\t/$/;s/\t/$/;s/\t/$/' > ${out}_ratio_all_wanted_header.tmp
        cat ${out}_odd_ratio.tmp ${out}_even_ratio.tmp ${out}_single_ratio.tmp |\
        sed 's/\t/$/;s/\t/$/;s/\t/$/;s/\t/$/' | sort -k1,1V > ${out}_ratio_all_wanted_noheader.tmp
        awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$89"\t"$154}' ${out}_final.hg38_multianno.txt |\
        sed 's/\t/$/;s/\t/$/;s/\t/$/;s/\t/$/' > ${out}_annovar_wanted.tmp
        head -1 ${out}_annovar_wanted.tmp > ${out}_annovar_wanted_header.tmp
        sed 1d ${out}_annovar_wanted.tmp | sort -k1,1V > ${out}_annovar_wanted_noheader.tmp
        cut --complement -f 16,18,19,23,24,25,26,27,28,34 ${out}_final.hg38_multianno.txt.intervar | tr -d "#" |\
        sed 's/\t/$/;s/\t/$/;s/\t/$/;s/\t/$/' > ${out}_intervar_wanted.tmp
        head -1 ${out}_intervar_wanted.tmp > ${out}_intervar_wanted_header.tmp
        sed 1d ${out}_intervar_wanted.tmp | sed 's/^/chr/' | sort -k1,1V > ${out}_intervar_wanted_noheader.tmp
        join -t $'\t' -1 1 -2 1 ${out}_ratio_all_wanted_header.tmp ${out}_annovar_wanted_header.tmp > ${out}_semi_final_header.tmp
        join -t $'\t' -1 1 -2 1 ${out}_semi_final_header.tmp ${out}_intervar_wanted_header.tmp > ${out}_final_header.tmp
        join -t $'\t' -1 1 -2 1 ${out}_ratio_all_wanted_noheader.tmp ${out}_annovar_wanted_noheader.tmp >  ${out}_semi_final_noheader.tmp
        join -t $'\t' -1 1 -2 1 ${out}_semi_final_noheader.tmp ${out}_intervar_wanted_noheader.tmp > ${out}_final_noheader.tmp
        cat ${out}_final_header.tmp ${out}_final_noheader.tmp | tr "$" "\t" |\
        sort -k1,1V -k2,2n -k3,3n -k4 -k5 > ${out}_final_hg38_selected.intervar
        rm "$dir"/*.tmp
    else
        echo "All outputs exists!"
    fi
done