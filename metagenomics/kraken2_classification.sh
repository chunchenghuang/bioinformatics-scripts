##blast for cross validating kraken2's result
#kraken classification file
dir="."
db_dir="/home/cch/kraken_database"
cpu="24"
READ_LEN="250"

mkdir $dir/before_fastqc

fastqc -t 20 $dir/*.fastq.gz -o $dir/before_fastqc

for list_R1 in $(find $dir -name "*R1*.fastq.gz"|sort); do
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

for list_R1 in $(find $dir -maxdepth 1 -name "*_R1*_paired.fastq.gz"|sort); do
list_R2=${list_R1//_R1/_R2};
output=${list_R1//_R1*.fastq.gz/};
microbial_output_R1=${list_R1//.fastq.gz/};
microbial_output_R2=${microbial_output_R1//_R1/_R2};
kraken2 --report ${output}_report.txt --threads $cpu --use-names --db $db_dir --paired --gzip-compressed \
--classified-out $output#_C.fastq.gz --unclassified-out $output#_UC.fastq.gz \
$list_R1 $list_R2 1> ${output}_classification.out
#ktImportTaxonomy -m 3 -t 5 ${output}_report.txt -o ${output}_krona.html
done

# bracken-build -d $db_dir -t $cpu -l ${READ_LEN}

for list_R1 in $(find $dir -maxdepth 1 -name "*_report.txt"|sort); do
output=${list_R1//_report.txt/};
bracken -d $db_dir -i ${output}_report.txt -o ${output}.bracken -r ${READ_LEN}
ktImportTaxonomy -m 3 -t 5 ${output}_report_bracken_species.txt -o ${output}_krona.html
done
# /opt/Bracken/analysis_scripts/combine_bracken_outputs.py --files *_bracken_species.txt -o all_report.bracken
# combine_kreports.py --display-headers --no-headers --report-file *_report_bracken_species.txt -o all.braken_reports --sample-names L1-0 L1-12 L1-14 L1-1 L1-2 L1-3 L1-5 L1-7 L2-0 L2-12 L2-14 L2-1 L2-2 L2-3 L2-5 L2-7 L3-0 L3-12 L3-14 L3-1 L3-2 L3-3 L3-5 L3-7 L4-0 L4-12 L4-14 L4-1 L4-2 L4-3 L4-5 L4-7 L5-0 L5-1 L5-2 L5-3 L5-5 L6-0 L6-1 L6-2 L6-3 L6-5 L7-0 L7-1 L7-2 L7-3 L7-5 N1-0 N1-12 N1-14 N1-1 N1-2 N1-3 N1-5 N1-7 N2-0 N2-12 N2-14 N2-1 N2-2 N2-3 N2-5 N2-7 N3-0 N3-12 N3-14 N3-1 N3-2 N3-3 N3-5 N3-7
# nohup kraken-biom *_report.txt -o kraken.biom &
kraken-biom *_report_bracken_species.txt -o bracken.biom