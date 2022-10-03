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

if [[ "$org" == *'mm'* ]]; then
    orgfull='Mus_musculus'
elif [[ "$org" == *'rn'* ]]; then
    orgfull='Rattus_norvegicus'
elif [[ "$org" == *'hg'* ]]; then
    orgfull='Homo_sapiens'
elif [[ "$org" == *'ce'* ]]; then
    orgfull='Caenorhabditis_elegans'
elif [[ "$org" == *'sacCer'* ]]; then
    orgfull='Saccharomyces_cerevisiae'
elif [[ "$org" == *'dm'* ]]; then
    orgfull='Drosophila_melanogaster'
elif [[ "$org" == *'danRer'* ]]; then
    orgfull='Danio_rerio'
else
    echo 'ERROR: wrong organism name provided (USCS version)'
    exit 1
fi

refgen="/home/myrto/work/refs/${orgfull}/UCSC/${org}/"

echo $org
echo $orgfull
echo $input
echo $refgen
