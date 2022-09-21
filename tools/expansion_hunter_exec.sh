#!bin/bash
nohup ExpansionHunter --reads NA04025_S1.bam \
    --reference /home/cch/hg19_reference/UCSC_hg19/hg38.fa \
    --variant-catalog /opt/ExpansionHunter-v3.2.2-linux_x86_64/variant_catalog/hg19/variant_catalog.json \
	--sex female \
    --output-prefix NA04025 &