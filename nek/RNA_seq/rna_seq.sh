#!/bin/bash

refgen=$1
input=$2
mkdir -p fqs trimmed_fqs STAR bais bws stranded_bws

while IFS="" read -r line || [ -n "$line" ]
        do

        date +'%m/%d/%Y %r'
        echo "! starting on ${line} !"

        #curl -X GET \
        #"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${line}&fields=fastq_ftp&result=read_run" |
        #cut -f 2 |
        #sed -n 2p |
        #tr \; "\n" |
        #wget -P ./fqs -i - &&
        #gunzip ./fqs/*.gz

        # trim off primers
        #trim_galore -q 20 -o ./trimmed_fqs --paired --nextera --fastqc ./fqs/${line}*

        # allign using STAR
        STAR    --genomeDir ~/work/genomes/STAR/${refgen} \
                --readFilesIn ./trimmed_fqs/${line}*.fq \
                --outFileNamePrefix ./STAR/${line}_ \
                --outSAMtype BAM SortedByCoordinate

        # uscs viewing stuff

        # indexing, what does bai stand for ?
        samtools index "./STAR/${line}_Aligned.sortedByCoord.out.bam"

        # creating bigwig file
        bamCoverage -b "./STAR/${line}_Aligned.sortedByCoord.out.bam" \
		-o "./bws/aligned_${line}.bw"

        bamCoverage -b "./STAR/${line}_Aligned.sortedByCoord.out.bam" \
                -o "./stranded_bws/aligned_fw_${line}.bw" \
                -p 24 \
                --binSize 10 \
                --normalizeUsing RPKM \
                --maxFragmentLength 0 \
                --filterRNAstrand=forward

        bamCoverage -b "./STAR/${line}_Aligned.sortedByCoord.out.bam" \
                -o "./stranded_bws/aligned_rv_${line}.bw" \
                -p 24 \
                --binSize 10 \
                --normalizeUsing RPKM \
                --maxFragmentLength 0 \
                --filterRNAstrand=reverse

        done < $input
