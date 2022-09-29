#!/bin/bash

mkdir -p fqs trimmed_fqs bowtie macs macs/sorted macs/merged macs/bed6 bws stranded_bws
input=$1

while IFS="" read -r line || [ -n "$line" ]; do

        date +'%m/%d/%Y %r'
        echo "starting on \"${line}\"... "

        curl -X GET \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${line}&fields=fastq_ftp&result=read_run" |
        cut -f 2 |
        sed -n 2p |
        tr \; "\n" |
        wget -P ./fqs -i - &&
        gunzip ./fqs/*.gz

        # trim off primers
        trim_galore -o ./trimmed_fqs \
                --nextera \
                --fastqc \
                --paired \
                -q 20 \
                ./fqs/${line}*

        # alignment, requires bowtie2-build
        # samtools sort requires the @SQ headers
        bowtie2 --local \
                --very-sensitive-local \
                -X 3000 \
                -x ~/work/genomes/bowtie2_indexes/rn6 \
                -1 ./trimmed_fqs/${line}_1_val_1.fq \
	        -2 ./trimmed_fqs/${line}_2_val_2.fq |
               	samtools sort -O BAM - > ./bowtie/${line}.sorted.bam

        # peak calling
        macs2 callpeak \
                --treatment ./bowtie/${line}.sorted.bam \
                --qvalue 0.001 \
                --gsize 2729860805 \
                --format BAMPE \
                --outdir ./macs \
                --name ${line}

        # sort and merge peaks
        sort -k1,1 -k2,2n ./macs/${line}_peaks.narrowPeak > ./macs/sorted/${line}_sortedPeaks.bed
        bedtools merge -i ./macs/sorted/${line}_sortedPeaks.bed -d 1000 > ./macs/merged/${line}_sortedMergedPeaks.bed

        # BED6+4 -> BED
        cut -f 1-6 ./macs/merged/${line}_sortedMergedPeaks.bed  > ./macs/bed6/${line}_peaks.6.bed

        done < $input
