#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: sh <scriptname>.sh <orgname> <input>"
    echo "       <orgname> organism name (USCS version)"
    echo "       <input> .txt file containing prefexes of files"
    echo ""
    echo "Other Options: use sh -x so every command is printed"
    exit 0
fi

org=$1
orgfull=''
input=$2

if [ "$org" == 'mm10' ]; then
    orgfull='Mus_musculus'
elif [ "$org" == 'rn6' ]; then
    orgfull='Rattus_Norvegicus'
else
    echo 'ERROR: wrong organism name provided (USCS version)'
    exit 0
fi

refgen="/home/myrto/work/refs/${orgfull}/UCSC/${org}/"

echo $org
echo $orgfull
echo $input
echo $refgen
