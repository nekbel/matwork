#!/bin/bash
# SCRIPTS USED:
#       orgUCSC
source ~/.nek_bash_commands.sh
echo "Imported source commands"

ref_build="${1:-}"
input="${2:-}"
gsize="${3:-}"
# rn need gsize 2729860805
# else $(echo "${ref_build}" | head -c2)

# check number of arguments given
if [ "$#" -ne 2 ]; then
    echo "Usage: bash <name>.sh <ref_build> <input>"
    echo "      <ref_build> organism reference build name (USCS version)"
    echo "      <input> .txt list of SRR codes"
    echo ""
    echo "Options: run script with -x option"
    exit 1
fi

gen_path=`orgUCSC ${ref_build}`
echo "Path to genome used: $gen_path"

echo "Creating directories (if they don't already exist)"
mkdir -p fqs trimmed_fqs bowtie2 macs macs/sorted macs/merged macs/bed6 bws stranded_bws

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
            -x ${gen_path}/Sequence/Bowtie2Index/genome \
            -1 ./trimmed_fqs/${line}.fq \
            -2 ./trimmed_fqs/${line}.fq |
            samtools sort -O BAM - > ./bowtie2/${line}.sorted.bam

    # peak calling
    macs2 callpeak \
            --treatment ./bowtie/${line}.sorted.bam \
            --qvalue 0.001 \
            --gsize ${gsize} \
            --format BAMPE \
            --outdir ./macs \
            --name ${line}

    # sort and merge peaks
    sort -k1,1 -k2,2n ./macs/${line}_peaks.narrowPeak > ./macs/sorted/${line}_sortedPeaks.bed
    bedtools merge -i ./macs/sorted/${line}_sortedPeaks.bed -d 1000 > ./macs/merged/${line}_sortedMergedPeaks.bed

    # BED6+4 -> BED
    cut -f 1-6 ./macs/merged/${line}_sortedMergedPeaks.bed  > ./macs/bed6/${line}_peaks.6.bed

    done < $input
