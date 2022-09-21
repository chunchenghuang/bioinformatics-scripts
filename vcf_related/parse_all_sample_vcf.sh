#!/bin/bash
vcf_AF_ratio_table () {
all_vcf_file=${1//.vcf/};
sed '/##/d' ${all_vcf_file}.vcf | grep "#" | awk '{$3=$6=$7=$8=$9=""; print $0}' > ${all_vcf_file}.header
sed '/#/d' ${all_vcf_file}.vcf | awk '{$3=$6=$7=$8=$9=""; print $0}' > ${all_vcf_file}.stripped
column_number=$(awk '{print NF;exit}' ${all_vcf_file}.stripped)
for i in $(seq 5 $column_number); do
awk '{print $5}' ${all_vcf_file}.stripped | cut -d ":" -f 3,5 | awk -F ":" -v d="$DP" '{if($1 < d){$0=".:."}{print $2}}' > $i.DP
awk -F "," '{if($0 != "."){$0=$2/($2+$1)}{print $0}}' $i.DP > $i.ratio
paste ${all_vcf_file}.stripped $i.ratio | awk '{$5=""; print $0}' > $i.tmp
rm ${all_vcf_file}.stripped
mv $i.tmp ${all_vcf_file}.stripped
done
cat ${all_vcf_file}.header ${all_vcf_file}.stripped | tr -s " " "\t" > ${all_vcf_file}_ratio_table.tmp
awk '{if($(NF-5)!="."||$(NF-4)!="."||$(NF-3)!="."||$(NF-2)!="."||$(NF-1)!="."||$NF!="."){print $0}}' ${all_vcf_file}_ratio_table.tmp > ${all_vcf_file}_ratio_table.txt
rm *.DP
rm *.ratio
rm *.header
rm *.stripped
rm *.tmp
}

#usage description
usage="$(basename "$0") [options] <.vcf> - to calculate allele ratio from AD values in vcf file

options:
    -h/--help           show this help text
    -d/--depth <int>    depth filter [50];
                        Ratio smaller than this value will be dropped to none represents by '.'
                        If all ratios are dropped across all samples, the variant will be automatically dropped
    <.vcf>              input a vcf file"

#set getopt into options
options=$(getopt --options hd: --long help,depth: --name "$0" -- "$@")

#set if exit status does not equals to 0 (failure)
if [ $? != 0 ] ; then echo "Failed to parse options...exiting." >&2 ; exit 1 ; fi

#eval set together makees the $options work and set as script arguments while -- shifts all arguments so that the input would all start from $1
eval set -- "$options"

DP="50" #default

#option conditions setting
while true; do
  case "$1" in
    -h | --help) echo "$usage"
        exit
        ;;
    \?) echo "Unknown parameter, see usage:" >&2 #>&2 means send everything to stderr and exit 1 means failure
        echo "$usage" >&2
        exit 1
        ;;
    --) shift ## -- allows you to input strings similar to options like -h for your function later on to process by shift
        break ## break just exit the loop before the loop's ending and continue the script but exit exits whole script
        ;;
    -d | --depth) DP=$2
        shift 2
        ;;
  esac
done

#input condition setting
if [[ -z $1 || $1 == *[[:space:]]* ]];then #-z stands for if it's empty and *[[:space:]]* represents continuous spaces infront and after
    echo "Please input either file or parameter" >&2
    exit 1
# [[ $var ]] indicates a variable is set, [[ ! $var ]] indicates a variable is unste
elif [[ $2 ]];then
    echo "Please only input one file or parameters" >&2
    exit 1
# [[ $var == z* ]] means variable matches with a string starts with z, [[ $var == z* ]]
elif [[ $1 != *.vcf ]];then
    echo "Please input a vcf file" >&2
    exit 1
else 
    vcf_AF_ratio_table $1
fi