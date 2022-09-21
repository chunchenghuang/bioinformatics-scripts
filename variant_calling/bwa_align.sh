#!bin/bash
dir="."
ref_dir="/home/cch/hg38_reference/UCSC_hg38"
ref_file=${ref_dir}/"hs38DH.fa"
#target region
regions=${dir}/"xgen-exome-research-panel-v2-targets-hg38.bed"

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
        bamtools sort -in $sample_name.bam -out $sample_name.bam
        bamtools index -in ${sample_name}.bam
        rm $sample_name*_L*.bam
    fi
done