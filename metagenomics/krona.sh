dir="."
for report in $(find $dir -maxdepth 1 -name "*_report.txt"|sort -V);do
    report=${report//_report.txt/};
    #ktImportTaxonomy -m 3 -t 5 ${report}_report.txt -o ${report}_krona.html
    kreport2mpa.py --percentages --display-header --no-intermediate-ranks -r ${report}_report.txt -o ${report}_mpa_p.txt
    sort -k 1b,1 ${report}_mpa_p.txt > ${report}_mpa_p_sorted.txt
    kreport2mpa.py --display-header --no-intermediate-ranks -r ${report}_report.txt -o ${report}_mpa.txt
    sort -k 1b,1 ${report}_mpa.txt > ${report}_mpa_sorted.txt
done

/home/cch/scripts/join_kraken2_report.sh *_mpa_p_sorted.txt > combined_mpa_p.txt
/home/cch/scripts/join_kraken2_report.sh *_mpa_sorted.txt > combined_mpa.txt

rm *_mpa_p.txt