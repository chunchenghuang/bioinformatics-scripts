#!/bin/bash

#environment setting
dir="."
ref_dir="/home/cch/hg38_reference/grch38_snp_tran"
cpu="7"
matrix_name="MPS"
#using downloaded index, input index prefix
index=$ref_dir/"genome_snp_tran"
gtf_file=$ref_dir/"Homo_sapiens.GRCh38.84.gtf"
gtf_file_name=${gtf_file//.gtf/}

for f in $ref_dir/*.ht*; do
    if [ ! -e "$f" ]; then
        #using raw .fa genome for indexing
        ref_file=$(find $ref_dir -name "*.fa")
        ref_file_name=${ref_file//.fa/}
        gtf_file=$(find $ref_dir -name "*.gtf")
        gtf_file_name=${gtf_file//.gtf/}
        index=${ref_file_name}_index
        hisat2-build -p $cpu $ref_file ${ref_file_name}_index
        #transcript option
        #extract_splice_sites.py $gtf_file > ${gtf_file_name}.ss
        #extract_exons.py $gtf_file > ${gtf_file_name}.exon
        #hisat2-build -p 20 --ss ${gtf_file_name}.ss --exon ${gtf_file_name}.exon $ref_file ${ref_file_name}_index
        #snp option (lack of RAM)
    fi
done

#loops
for list_R1 in $(find $dir -name "*R1*_paired.fastq.gz"|sort); do
    list_R2=${list_R1//_R1/_R2};
    output=${list_R1//_R1*/};
    sample_name=$(echo "$list_R1"|cut -d "_" -f 1)
    if [ ! -e "${sample_name}.sam" ] && [ ! -e "${sample_name}.bam" ]; then
        hisat2 -p $cpu -x ${index} -1 ${list_R1} -2 ${list_R2} -S ${output}.sam --summary-file ${output}.txt --new-summary &
        [ $(jobs|wc -l) -ge $(expr $(nproc) / $cpu) ] && wait
    fi
done
wait

#sam to bam
for sam in $(find $dir -name "*.sam"|sort); do
    bam=${sam//.sam/.bam};
    if [ ! -e "$bam" ]; then
        samtools view -S -b $sam > $bam
        rm $sam
    fi
done

#merge, sort bam files
for sample_name in $(find $dir -name "*.bam"|cut -d "_" -f 1 | sort | uniq); do
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
    fi
done

#featureCounts to produce read counts
featureCounts -T $cpu -t gene -g gene_id -a $gtf_file -o $matrix_name.out *.bam

#produce count matrix with gene name
#http://genomespot.blogspot.com/2015/01/generate-rna-seq-count-matrix-with.html
cut -f1,7- $matrix_name.out | sed 1d > $matrix_name.mx
grep -w gene $gtf_file | cut -d '"' -f2,6 | tr '"' '\t' | sort -k 1b,1 > ENS2genename.tmp
head -1 $matrix_name.mx > matrix_header.tmp
sed 1d $matrix_name.mx | sort -k 1b,1 |\
join -1 1 -2 1 ENS2genename.tmp - |\
tr ' ' '\t' | sed 's/\t/_/' |\
cat matrix_header.tmp - | sed 's/.bam//g'> ${matrix_name}_counts.mx

#TPMCalculator to produce TPM and FPKM
TPMCalculator -e -g $gtf_file -d $dir

for TPM_files in $(find $dir -name "*_genes.out");do
    sample_name=${TPM_files//_genes.out/}
    cut -f 1,7 $TPM_files | sed 1d | sort -k 1b,1 > ${sample_name}_TPM.tmp
    awk '{print $2}' ${sample_name}_TPM.tmp > ${sample_name}_TPM_only.tmp
done

first_sample=$(ls *_TPM.tmp|sort|head -n 1)
awk '{print $1}' $dir/$first_sample | tr "#" "\t" > table_TPM_row.header
join -1 1 -2 1 ENS2genename.tmp table_TPM_row.header | tr ' ' '\t' |\
awk '{print $1"#"$3"_"$2}' | sed 's/#_/_/' > table_TPM_row_genenames.header
sample_name_list=$(ls *_TPM_only.tmp|sort|cut -d "_" -f 1|tr "\n" "\t")
echo -e "Gene_Id\t$sample_name_list" > table_TPM_column.header
sample_list=$(ls *_TPM_only.tmp|sort|tr "\n" " ")
paste table_TPM_row_genenames.header $sample_list > ${matrix_name}_TPM_mx.tmp
cat table_TPM_column.header ${matrix_name}_TPM_mx.tmp > ${matrix_name}_TPM.mx

# ##estimate maximum depth
# column_number=$(awk '{print NF;exit}' ${matrix_name}_counts.mx)
# for i in $(seq 2 $column_number); do
#     mx_sample_header=$(awk -v i="$i" '{print $i}' ${matrix_name}_counts.mx | head -1)
#     awk -v i="$i" '{print $i}' ${matrix_name}_counts.mx | sed 1d| sort -nr |\
#     head -1 > $mx_sample_header.number.tmp
# done
# max_depth=$(cat *.number.tmp | sort -nr | head -1)

# ##variant calling from bcftools
# #https://www.biostars.org/p/335121/
# for bam_files in $(find $dir -name "*.bam" | sort); do
#     sample_name=${bam_files//.bam/}
#     bcftools mpileup -d $max_depth --threads $cpu -Ou -f $ref_file $sample_name.bam |\
#     bcftools call -mv -Ov -o $sample_name.vcf
# done

#clean-ups
rm $dir/*.tmp
rm $dir/*.header