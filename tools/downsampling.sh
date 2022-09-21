#!/bin/bash
sample_size="100000"
seed=$(shuf -i 1-100 -n 1)
for list_R1 in $(find . -maxdepth 1 -name "*R1*.fastq.gz");do
list_R2=${list_R1//_R1/_R2}
output_R1=${list_R1//.fastq.gz/_sub.fastq}
output_R2=${list_R2//.fastq.gz/_sub.fastq}
seqtk sample -s $seed $list_R1 $sample_size > $output_R1
seqtk sample -s $seed $list_R2 $sample_size > $output_R2
seed=$(($seed+1))
done
gzip *.fastq