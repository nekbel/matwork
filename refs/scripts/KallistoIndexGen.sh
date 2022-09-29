#!/bin/bash

if [ -z $1 ] ; then
    echo 'give path to file'
    exit 1
else
    echo 'running for file ${1}...'
fi

ff=$1
kallisto index -i ~work/genomes/kallisto/${ff}.index $ff
