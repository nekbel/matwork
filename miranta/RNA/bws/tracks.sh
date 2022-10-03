#!/bin/bash

touch tracks.txt
color=''
n=''
url=$(pwd | sed "s|/home/myrto/work|https://raw.githubusercontent.com/nekbel/matwork/main|")
for fn in ./*.bw; do
    eval "n=$(echo $fn | cut -c3-)"

    if [[ S"$n" = S*'WT'* ]]; then
        eval "color=255,0,0"
    else
        eval "color=0,255,0"
    fi
    echo "track type=bigWig name=${n%.*} description=${n%.*} visibility=full bigDataUrl=${url}/${n} color=${color}" >> tracks.txt
    done
