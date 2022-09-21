for g in CYP2D6 CYP2A6 CYP2C19 CYP2C8 CYP2C9 CYP3A4 CYP3A5 CYP4F2 TPMT DPYD; do
    parallel --eta aldy genotype -p illumina {} -g $g -o {/}-${g}.log ::: *.bam;
done