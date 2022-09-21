#!/bin/bash

# multijoin - join multiple files

join_rec() {
f1=$1; f2=$2; shift 2;
if [ $# -gt 0 ]; then
    join -t $'\t' -j 1 -a1 -a2 -o auto -e 0 "$f1" "$f2" | join_rec - "$@";
else
    join -t $'\t' -j 1 -a1 -a2 -o auto -e 0 "$f1" "$f2";
 fi
}

#usage description
usage="$(basename "$0") [options] <files> - to join multiple files 
options:
    -h/--help           show this help text
    <files>             need to input at least two files"

#set getopt into options
options=$(getopt --options hd: --long help: --name "$0" -- "$@")

#set if exit status does not equals to 0 (failure)
if [ $? != 0 ] ; then echo "Failed to parse options...exiting." >&2 ; exit 1 ; fi

#eval set together makees the $options work and set as script arguments while -- shifts all arguments so that the input would all start from $1
eval set -- "$options"

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
  esac
done

#input condition setting
if [[ -z $1 || $1 == *[[:space:]]* ]];then #-z stands for if it's empty and *[[:space:]]* represents continuous spaces infront and after
    echo "Please input either file or parameter" >&2
    exit 1
elif [[ $# < 2 ]]; then
    echo "Must input more than two files"
    exit 1
else 
    join_rec $@
fi