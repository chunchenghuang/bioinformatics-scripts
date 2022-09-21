#!/bin/bash
dir="."
for f in *.vcf; do
#echo $f
bgzip $f
#it overwrites original file
tabix -p vcf ${f}.gz
done

vcf-merge -d *.vcf.gz > all.vcf