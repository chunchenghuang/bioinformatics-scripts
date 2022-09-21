#!/bin/bash
present_dir="."
output_dir=$present_dir
for i in $(find $present_dir -name "*L001*.fastq.gz"|sort);do
o=${i//$present_dir/$output_dir};
o=${i//_L001_/_all_};
b=${i//_L001_/_L002_};
c=${i//_L001_/_L003_};
d=${i//_L001_/_L004_};
cat $i $b $c $d > ${o};
done