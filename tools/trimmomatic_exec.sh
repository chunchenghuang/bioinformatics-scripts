# #!bin/bash
dir="."
for list_R1 in $(find $dir -name "*R1*.fastq.gz"|sort); do
list_R1=${list_R1//.fastq.gz/}
list_R2=${list_R1//_R1/_R2};
#testing
#  echo ${list_R1}
#  echo ${list_R2}
#  echo ${list_R1}_25.fastq.gz
#  echo ${list_R1}_25_paired.fastq.gz
#  echo ${list_R1}_25_unpaired.fastq.gz
java -jar /opt/app/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 ${list_R1}.fastq.gz  ${list_R1}_25.fastq.gz  HEADCROP:25
java -jar /opt/app/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 ${list_R2}.fastq.gz  ${list_R2}_25.fastq.gz  HEADCROP:25
java -jar /opt/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 ${list_R1}_25.fastq.gz ${list_R2}_25.fastq.gz ${list_R1}_25_paired.fastq.gz ${list_R1}_25_unpaired.fastq.gz ${list_R2}_25_paired.fastq.gz ${list_R2}_25_unpaired.fastq.gz LEADING:30 TRAILING:30 MINLEN:100 SLIDINGWINDOW:4:20
done
rm "$dir"/*_25.fastq.gz