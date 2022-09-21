#!/bin/bash
bam_qc () {
declare DP="$1"
dir="$PWD"
echo "cpu=$cpu"
echo "depth=$DP"
echo "working directory=$dir"
echo "QC starting..."

if [ $(ls | grep .bam$ | wc -l) == 0 ]; then
    echo "ERROR: No BAM file(s) detected!!"
    echo "Must contain BAM file(s) and only one BED file in the directory!!"
    exit 1
fi
if [ $(ls | grep .bed$ | wc -l) != 1 ]; then
    echo "ERROR: The working directory must contain only one BED file!!"
    echo "Must contain BAM file(s) and only one BED file in the directory!!"
    exit 1
fi

bed=$(find $dir -name "*.bed")

for bam_file in $(find $dir -name "*.bam"); do
    bam_file=${bam_file//.bam/}
    bam_file=${bam_file//_mkdup*/}
    if [ ! -e "${bam_file}_mkdup.bam" ]; then
        java -Xmx100g -jar /opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar MarkDuplicates \
        -I ${bam_file}.bam -O ${bam_file}_mkdup.bam -M ${bam_file}_metrics.txt
    else
    echo "Marked duplicate BAM file exists!"
    fi
    if  [ ! -e "${bam_file}_QC.csv" ]; then
        samtools stats -@ $cpu ${bam_file}_mkdup.bam | grep ^SN | cut -f 2- > ${bam_file}_stats.tmp
        samtools stats -t $bed -@ $(expr $cpu - 1) ${bam_file}_mkdup.bam | grep ^SN | cut -f 2- > ${bam_file}_stats_region.tmp &
        samtools depth -aa -b $bed -d 0 -Q 10 ${bam_file}_mkdup.bam > ${bam_file}_DP.tmp
        wait
        samtools view -@ $cpu -L $bed -b ${bam_file}_mkdup.bam > ${bam_file}_mkdup_region.bam.tmp
        samtools fastq -@ $cpu ${bam_file}_mkdup_region.bam.tmp 1> ${bam_file}_mkdup.fq.tmp 2> fq.err.tmp
        fqtools qualtab ${bam_file}_mkdup.fq.tmp | awk '{print $2}'> ${bam_file}_qual_table.tmp
        Rscript ~/scripts/bam_qc_server.R $dir ${bam_file}_stats.tmp ${bam_file}_stats_region.tmp \
        ${bam_file}_qual_table.tmp ${bam_file}_DP.tmp $DP $bam_file
        #rm *.tmp
    else
    echo "QC report exists!"
    fi
done
}

usage="$(basename "$0") [options] <depth[200]> - bam file QC and output regions whose depth is under the specified value.

options:
    -h/--help           show this help text
    -c/--cpu <int>      number of cpu to be used (cpu number cannot be lower than 2) [4]"

#set getopt into options

if [[ "$1" =~ ^[-]([0-9]*[.])?([0-9]+)?$ ]] ; then echo "The input must be positive and numeric!!!" >&2 ; exit 1; fi

options=$(getopt --options hc: --long help,cpu: --name "$0" -- "$@")

#set if exit status does not equals to 0 (failure)
if [ $? != 0 ] ; then echo "Failed to parse options...exiting." >&2 ; exit 1 ; fi

#eval set together makees the $options work and set as script arguments while -- shifts all arguments so that the input would all start from $1
eval set -- "$options"

cpu="4" #default

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
    -c | --cpu) cpu="$2"
        shift 2
        ;;
  esac
done

if [[ $cpu < 2 ]]; then echo "The cpu number cannot be lower than 2, used 2 instead..."; fi

#input condition setting
if [[ -z $1 || $1 == *[[:space:]]* ]];then #-z stands for if it's empty and *[[:space:]]* represents continuous spaces infront and after
    bam_qc 200 #defult 200
# [[ $var ]] indicates a variable is set, [[ ! $var ]] indicates a variable is unset
elif [[ $2 ]];then
    echo "Please only input one positive numeric value." >&2
    exit 1
# [[ $var == z* ]] means variable matches with a string starts with z, [[ $var == z* ]]
elif ! [[ "$1" =~ ^[+]?([0-9]*[.])?([0-9]+)?$ ]] ;then
    echo "The input must be positive and numeric!!!" >&2
    exit 1
else 
    bam_qc $1
fi