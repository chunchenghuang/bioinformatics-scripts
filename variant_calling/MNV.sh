#!bin/bash
java -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar SelectVariants \
    -V CM103_markduplicates_bqsr_filtered_pass.vcf \
    -select-type SNP \
    -O snvs.vcf

grep -v "#" snvs.vcf | awk '{print $1"."$2"."$4"."$5}' > input_SNVs.txt

/opt/MAC/MAC_v1.2.pl -i input_SNVs.txt\
 -bam CM103.bam\
 -r /home/public/hg19_1000G/ucsc.hg19.fasta\
 -annotator annovar\
 -annovar_annotate_variation /home/hao/Downloads/annovar/annotate_variation.pl\
 -annovar_coding_change /home/hao/Downloads/annovar/coding_change.pl\
 -annovar_refgene /home/hao/Downloads/annovar/humandb/hg19_refGene.txt\
 -annovar_refmrna /home/hao/Downloads/annovar/humandb/hg19_refGeneMrna.fa