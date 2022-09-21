#!bin/bash
kmer="42"
dir="."
ref_dir="/home/public/hg38_1000G"
ref_file=${ref_dir}/"Homo_sapiens_assembly38.fasta"
for list_R1 in $(find $dir -name "*R1*.fastq.gz"|sort); do
    output=$(echo $list_R1 | cut -d "_" -f 1)
    output=$(echo ${output}_k${kmer})
    list_R1=${list_R1//.fastq.gz/}
    list_R2=${list_R1//_R1/_R2}
    ##assemble contigs and phase for de-brujin graph
    platanus_allee assemble -k $kmer -t 20 -f ${list_R1}.fastq.gz ${list_R2}.fastq.gz -o $output 2>assemble.log &
    platanus_allee phase \
        -c ${output}_contig.fa ${output}_junctionKmer.fa \
        -IP1 ${list_R1}.fastq.gz ${list_R2}.fastq.gz \
        -o $output \
        2>phase.log
    ##testing platanus allee consensus for scaffolds
    platanus_allee consensus -t 20 \
        -c ${output}_primaryBubble.fa ${output}_nonBubbleHomoCandidate.fa \
        -IP1 ${list_R1}.fastq.gz ${list_R2}.fastq.gz \
        -o $output \
        2>consensus.log
    ##testing platanus allee for gap closing
    platanus_allee gap_close -t 20 \
        -c ${output}_consensusScaffold.fa \
        -IP1 ${list_R1}.fastq.gz ${list_R2}.fastq.gz \
        2>gap_close.log &
    #map assembled contigs for visualisation
    bwa mem -M -t 20 -R '@RG\tID:foo\tSM:bar' $ref_file ${output}_allPhaseBlock.fa | samtools view -@ 20 -S -b - > ${output}_allPhaseBlock_1.bam
    bamtools sort -in ${output}_allPhaseBlock_1.bam -out ${output}_allPhaseBlock_1.bam
    bamtools index -in ${output}_allPhaseBlock_1.bam
    bwa mem -M -t 20 -R '@RG\tID:foo\tSM:bar' $ref_file ${output}_allPhaseBlock.fa | samtools view -@ 20 -S -b - > ${output}_allPhaseBlock_2.bam
    bamtools sort -in ${output}_allPhaseBlock_2.bam -out ${output}_allPhaseBlock_2.bam
    bamtools index -in ${output}_allPhaseBlock_2.bam
done
