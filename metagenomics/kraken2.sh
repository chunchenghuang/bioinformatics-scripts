dir="."
db_dir="/home/cch/kraken_database"
cpu="24"
for list_R1 in $(find $dir -maxdepth 1 -name "*_R1*_paired.fastq.gz"|sort); do
list_R2=${list_R1//_R1/_R2};
output=${list_R1//_R1*.fastq.gz/};
microbial_output_R1=${list_R1//.fastq.gz/};
microbial_output_R2=${microbial_output_R1//_R1/_R2};
kraken2 --report ${output}_report.txt --threads $cpu --use-names --db $db_dir --paired --gzip-compressed \
--classified-out $output#_C.fastq.gz --unclassified-out $output#_UC.fastq.gz \
$list_R1 $list_R2 1> ${output}_classification.out
# grep -v "Homo sapiens" ${output}_classification.out | awk '{print $2}' > ${output}_human_seq.list
# seqtk subseq $list_R1 ${output}_human_seq.list 1>${microbial_output_R1}_human.fastq
# seqtk subseq $list_R2 ${output}_human_seq.list 1>${microbial_output_R2}_human.fastq
# seqtk seq -A ${microbial_output_R1}_microbial.fastq > ${microbial_output_R1}_microbial.fasta
# seqtk seq -A ${microbial_output_R2}_microbial.fastq > ${microbial_output_R2}_microbial.fasta
# blastn -query ${microbial_output_R1}_microbial.fasta \
# -db /home/cch/kraken_database/library/all_for_blast \
# -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > ${output}_microbial_blast_R1.outfmt6
# blastn -query ${microbial_output_R2}_microbial.fasta \
# -db /home/cch/kraken_database/library/all_for_blast \
# -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > ${output}_microbial_blast_R2.outfmt6
# awk '$3>95{print $2}' ${output}_microbial_blast_R1.outfmt6 | sort | uniq -c | tr ' ' '\t' | sed 's/^[ \t]*//g' | cut -d ',' -f1 | sed -e 's/|/\t/3; s/|/ /g' |sed 1d > ${output}_blast_R1.list
# awk '$3>95{print $2}' ${output}_microbial_blast_R2.outfmt6 | sort | uniq -c | tr ' ' '\t' | sed 's/^[ \t]*//g' | cut -d ',' -f1 | sed -e 's/|/\t/3; s/|/ /g' |sed 1d > ${output}_blast_R2.list
#extract the Mycobacterium sequences out
# grep Myco ${output}_classification.out | grep -v Mycoplasma | awk '{print $2}' > ${output}_mycobacterium_seq.list
# seqtk subseq $list_R1 ${output}_mycobacterium_seq.list 1>${microbial_output_R1}_mycobacterium.fastq
# seqtk subseq $list_R2 ${output}_mycobacterium_seq.list 1>${microbial_output_R2}_mycobacterium.fastq
# seqtk seq -A ${microbial_output_R1}_mycobacterium.fastq > ${microbial_output_R1}_mycobacterium.fasta
# seqtk seq -A ${microbial_output_R2}_mycobacterium.fastq > ${microbial_output_R2}_mycobacterium.fasta
# blastn -query ${microbial_output_R1}_mycobacterium.fasta \
# -db /home/cch/kraken_database/library/all_for_blast \
# -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > ${output}_mycobacterium_blast_R1.outfmt6 &
# blastn -query ${microbial_output_R2}_mycobacterium.fasta \
# -db /home/cch/kraken_database/library/all_for_blast \
# -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > ${output}_mycobacterium_blast_R2.outfmt6 &
#ktImportTaxonomy -m 3 -t 5 ${output}_report.txt -o ${output}_krona.html
done
