#!bin/bash

dir="."

mkdir $dir/before_fastqc

fastqc -t 20 $dir/*.fastq.gz -o $dir/before_fastqc

for list_R1 in $(find $dir -name "*all_R1*.fastq.gz"|sort); do
list_R1=${list_R1//.fastq.gz/}
list_R1=${list_R1//$dir\//}
list_R2=${list_R1//_R1/_R2};
java -Xmx100g -jar /opt/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 20 -phred33 \
 $dir/${list_R1}.fastq.gz \
 $dir/${list_R2}.fastq.gz \
 $dir/${list_R1}_paired.fastq.gz \
 $dir/${list_R1}_unpaired.fastq.gz \
 $dir/${list_R2}_paired.fastq.gz \
 $dir/${list_R2}_unpaired.fastq.gz \
 ILLUMINACLIP:/opt/Trimmomatic-0.38/adapters/TruSeq.fa:2:30:10:1:true \
 LEADING:30 TRAILING:30 MINLEN:100 SLIDINGWINDOW:4:20
done

mkdir $dir/after_fastqc

fastqc -t 20 $dir/*_paired.fastq.gz -o $dir/after_fastqc