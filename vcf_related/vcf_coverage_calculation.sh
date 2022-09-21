#!/bin/bash
usage="$(basename "$0") [-h] [<.genome.vcf/.gvcf>] - to extract depth from gvcf file

where:
    -h  show this help text
    <.genome.vcf/.gvcf> input a genome vcf file"

while getopts ':h' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    s) seed=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

parse_vcf_to_depth(){

genome_vcf_file=$1

##delete headers
stripped_genome_vcf=${genome_vcf_file//.genome.vcf/_stripped.genome.vcf}
sed '/#/d' ${genome_vcf_file}} > ${stripped_genome_vcf}

##parse by regions and bases
regions=$(grep "END=" ${stripped_genome_vcf})
bases=$(grep -v "END=" ${stripped_genome_vcf})

##BASES
##strip out all no call bases
low_MP_call_bases_sort=$(grep "\./." ${bases} | awk '{print $1"\t"$2"\t"$2"\t""0"}')
low_GQ_call_bases_sort=$(awk '{if($10=="."){print $0}}' ${bases} | awk '{print $1"\t"$2"\t"$2"\t""0"}')
base_calls=$(grep -v "\./." ${bases} | awk '{if($10!="."){print $0}}')

##separate variants and non-variants from base_calls
variants_sort=$(grep "DP=" ${base_calls} | awk '{print $1"\t"$2"\t"$2"\t"$8}'| sed -e 's/BaseQRankSum=.*;DP=//g' -e 's/DP=//g' -e 's/;.*//g')
non_variants_sort=$(grep -v "DP=" ${base_calls} | sed 's/\t/:/g'| awk -F ":" '{print $1"\t"$2"\t"$2"\t"$(NF-2)}')

##REGIONS
##strip out all no call bases regions
low_MP_call_regions_sort=$(grep "\./." ${regions} | awk '{print $1"\t"$2"\t"$8"\t""0"}' | sed 's/END=//g')
low_GQ_call_regions_sort=$(awk '{if($10=="."){print $0}}' ${regions} | awk '{print $1"\t"$2"\t"$8"\t""0"}' | sed 's/END=//g')
region_calls=$(grep -v "\./." ${regions} | awk '{if($10!="."){print $0}}')

##parse regions
regions_sort=$(sed -e 's/\t/:/g' -e 's/;LowMQ:/:/g' -e 's/;/:/g' ${region_calls} | awk -F ":" '{print $1"\t"$2"\t"$8"\t"$(NF-2)}' | sed 's/END=//g')

##concatenate all files, base and regions
total_gvcf_depth_unsorted=$(cat \
${low_MP_call_bases_sort} \
${low_GQ_call_bases_sort} \
${variants_sort} \
${non_variants_sort.txt} \
${low_MP_call_regions_sort} \
${low_GQ_call_regions_sort.txt} \
${regions_sort.txt} \
)


##check how to sort
##sort -V total_gvcf_depth_unsorted.txt | awk '{print $1}' | uniq
##sort -V total_gvcf_depth_unsorted.txt | awk '{print $0}' | head
##chromosomes and start position seemed fine using this sort -V, otherwise use sort -k 1,1 -k2,2n

##clean sorted output
total_gvcf_depth=${genome_vcf_file//.genome.vcf/_depth.bed}
sort -V total_gvcf_depth_unsorted.txt > ${total_gvcf_depth}.txt

}

parse_vcf_to_depth $1

##Validation
##calculate the base count that the coverage <20 by using start and stop position with conditions on depth column
#nohup awk 'BEGIN{sum=0}{if($NF < 20){sum+=($3-$2)+1}}END{print sum}' total_gvcf_depth.txt >  less_20_base_counts.txt &
#528499

##calculate total base count by using start and stop position with conditions on depth column
#nohup awk 'BEGIN{sum=0}{sum+=($3-$2)+1}END{print sum}' total_gvcf_depth.txt >  total_base_counts.txt &
#38852603

#echo "528499 / 38852603" | bc -l
#.01360266646741789732
##the percentage still not 0.005

##take out the low quality base calls
#nohup awk 'BEGIN{sum=0}{if($NF==0){sum+=($3-$2)+1}}END{print sum}' total_gvcf_depth.txt >  equal_0_base_counts.txt &
#375389
#echo "(528499-375389)/(38852603-375389)" | bc -l
#.00397923820575990766


##check <50X coverage percentage
#nohup awk 'BEGIN{sum=0}{if($NF >= 50){sum+=($3-$2)+1}}END{print sum}' total_gvcf_depth.txt >  more_50_base_counts &
#2523662
#echo "(2523662-375389)/(38852603-375389)" | bc -l
#echo "(1-0.05583234274706063697)*100" | bc -l
#94.41676572529393630300
#it's close to the BWA enrichment stat report