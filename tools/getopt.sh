#!/bin/bash

#define your function first

Hello () {
   echo $@
}

#set your help message

usage="$(basename "$0") [-h] <.genome.vcf/.gvcf> - to extract depth from gvcf file

where:
    -h  show this help text
    <.genome.vcf/.gvcf> input a genome vcf file"

#set your options with getopt with short ones and long ones
#add : behind your option name (parameters after --options) e.g. help: if you would need to add an argument after --help

options=$(getopt --options h --long help --name "$0" -- "$@")

#$@ expands to numerous positional parameters, starting from one, e.g. $1, $2, $3; $* is a parameter consisting all the parameters

#[ $? != 0 ] meanes if the exit satus do not equals to 0 

if [ $? != 0 ] ; then echo "The input options do not exist...exiting." >&2 ; exit 1 ; fi

#eval is when you defining a command into a variable and later on, if you want to use that command then you should use eval

#set sets the value of shell options and positional parameters for this script

#-- signifies the end of command options, after which only posistional parameters are accpeted

#-- this is useful because what we do is to sourcing a script inside of another script

eval set -- "$options"

#eval set together makees the $options work and set as script arguments while -- shifts 

#while true means loops forever, need to have commands like exit or break to break out the loop

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

if [[ -z $1 || $1 == *[[:space:]]* ]];then #-z stands for if it's empty and *[[:space:]]* represents continuous spaces infront and after
    echo "Please input file or parameters" >&2
    exit 1
else
    Hello $@
fi
