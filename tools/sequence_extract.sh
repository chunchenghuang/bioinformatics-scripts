#!bin/bash
dir="."
ref_dir="/home/public/hg38_1000G"
ref_file=${ref_dir}/"SMN1.fa"
for list_R1 in $(find $dir -name "*R1*.fastq.gz"|sort); do
    list_R1=${list_R1//.fastq.gz/}
    list_R2=${list_R1//_R1/_R2}
    sample_name=$(echo "$list_R1"|awk -F "/" '{print $NF}'|cut -d "_" -f 1)
    group_id=$(cat ${list_R1}_SMN.fastq | head -n 1 | tr -d "@" | cut -d ":" -f 3,4 | tr ":" ".")
    library=$(cat ${list_R1}_SMN.fastq | head -n 1 | awk -F ":" '{print $NF}')
    header=$(echo "@RG\tID:"$group_id"\tSM:"${sample_name}"\tPL:ILLUMINA\tLB:"$library)
    samtools view ${sample_name}.bam "chr5:70049612-70077592" "chr5:70925030-70953012" | awk '{print $1}' > ${sample_name}_SMN.list
    seqtk subseq $list_R1.fastq.gz ${sample_name}_SMN.list 1> ${list_R1}_SMN.fastq
    seqtk subseq $list_R2.fastq.gz ${sample_name}_SMN.list 1> ${list_R2}_SMN.fastq
    bwa mem -M -t 20 -R $header $ref_file ${list_R1}_SMN.fastq ${list_R2}_SMN.fastq |\
        samtools view -@ 20 -S -b - > ${sample_name}_SMN_realigned.bam
    bamtools sort -in ${sample_name}_SMN_realigned.bam -out ${sample_name}_SMN_realigned.bam
    bamtools index -in ${sample_name}_SMN_realigned.bam
done