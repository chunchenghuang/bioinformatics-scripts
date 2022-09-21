#!bin/bash
for i in $(find . -name "*_sub*");do
d=${i//_sub/}
# echo $i
# echo $d
mv $i $d
done