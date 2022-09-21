#!bin/bash
while IFS=$'\t' read -r -a line;
do
grep -P "\t${line}_" xgen-exome-research-panel-v2-targets-hg38.bed
done < target_gene.list > gene_list.bed