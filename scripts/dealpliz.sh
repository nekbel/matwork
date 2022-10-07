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



#master.sh -> RNA
#!/bin/bash

input="./test.txt"

grep -v '^#' "$input" | while IFS= read -r line
        do

        date +'%m/%d/%Y %r'
        echo "! starting on ${line} !"

        # download SRR files from NCBI
        prefetch "$line"

        # convert SRA to fastq files
        fastq-dump --split-files --outdir ./fqs "/home/myrto/ncbi/public/sra/${line}.sra"

        # trim off primers
        trim_galore -q 20 -o ./trimmed_fqs --paired --nextera --fastqc ./fqs/${line}*

        # allign using STAR
        STAR    --genomeDir ./ref_gen/gd \
                --readFilesIn ./trimmed_fqs/${line}*.fq \
                --outFileNamePrefix ./STAR/${line}_ \
                --outSAMtype BAM SortedByCoordinate

        # indexing
        samtools index "./STAR/${line}_Aligned.sortedByCoord.out.bam"
        mv ${line}_Aligned.sortedByCoord.out.bam ./bais

        # creating bigwig file
        bamCoverage -b "./STAR/${line}_Aligned.sortedByCoord.out.bam" -of bigwig -o "./bws/aligned_${line}.bw"

        bamCoverage -b "./STAR/${line}_Aligned.sortedByCoord.out.bam" \
                                -o "./stranded_bws/aligned_fw_${line}.bw" \
                                -p 24 \
                                --binSize 10 \
                                --normalizeUsingRPKM \
                                --maxFragmentLength 0 \
                                --filterRNAstrand=forward

        bamCoverage -b "./STAR/${line}_Aligned.sortedByCoord.out.bam" \
                                -o "./stranded_bws/aligned_rv_${line}.bw" \
                                -p 24 \
                                --binSize 10 \
                                --normalizeUsingRPKM \
                                --maxFragmentLength 0 \
                                --filterRNAstrand=reverse

        done

# kallisto.sh
#!/bin/bash

input=$1

mkdir -p kallisto

while IFS="" read -r line || [ -n "$line" ] ; do


        date +'%m/%d/%Y %r'
        echo "! starting on ${line} !"

        # kallisto
        kallisto quant \
                --rf-stranded \
                --index=/home/myrto/work/genomes/kallisto/rn6_cdna_kallisto_index \
                --output-dir=./kallisto/${line} \
                --threads=4 \
                ./trimmed_fqs/*${line}*.fq

        done < $input