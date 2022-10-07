#!/bin/bash

##########################
# Main rnaseq pipeline
##########################

# Miranta library type -> ISR Inward Stranded Reverse
# Carullo library type -> ISR

# SCRIPTS USED:
#       orgUCSC
source ~/.nek_bash_commands.sh
echo "Imported source commands"

ref_build="${1:-}"
input="${2:-}"

# check number of arguments given
if [ "$#" -ne 2 ]; then
    echo "Usage: bash rna_seq_miranta.sh <ref_build> <input>"
    echo "      <ref_build> organism reference build name (USCS version)"
    echo "      <input> .txt list of SRR codes"
    echo ""
    echo "Options: run script with -x option"
    exit 1
fi

gen_path=`orgUCSC ${ref_build}`
echo "Path to genome used: $gen_path"

echo "Creating directories (if they don't already exist)"
mkdir -p fqs trimmed_fqs STAR bais bws stranded_bws

while IFS="" read -r line || [ -n "$line" ]; do

    date +'%m/%d/%Y %r'
    echo "! starting on ${line} !"

    curl -X GET \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${line}&fields=fastq_ftp&result=read_run" |
    cut -f 2 |
    sed -n 2p |
    tr \; "\n" |
    wget -P ./fqs -i - &&
    gunzip ./fqs/*.gz

    # trim off primers
    trim_galore --paired --nextera --fastqc -q 20 -o ./trimmed_fqs ./fqs/${line}*

    # NOTE: “STAR’s default parameters are optimized for mammalian genomes. 
    # Other species may require significant modifications of some alignment 
    # parameters; in particular, the maximum and minimum intron sizes have 
    # to be reduced for organisms with smaller introns”.

    # allign using STAR
    STAR    --genomeDir ${gen_path}/Sequence/STARIndex \
		--readFilesIn ./trimmed_fqs/${line}*.fq \
		--outFileNamePrefix ./STAR/${line}_ \
		--outSAMtype BAM SortedByCoordinate

    # UCSC viewing stuff

    # indexing
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


##########################
# Manolo htseq count (gff)
##########################

echo $'\n'"Counting reads per feature..."

REFERENCE_ANNOTATION="/data/refs/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"

for sample in fastq/*.fastq;
do
	FILE_ID=${sample:6:4}
	echo $FILE_ID

# 3' QuantSeq
	htseq-count -f bam -s yes -i gene_id bam/$FILE_ID.bam $REFERENCE_ANNOTATION > $COUNTS_DIR/NGS$FILE_ID &

# Full RNA
	htseq-count -f bam -s no -i gene_id bam/$FILE_ID.bam $REFERENCE_ANNOTATION > $COUNTS_DIR/NGS$FILE_ID &

done;

htseq-count -f bam -s yes -i gene_id bam/$FILE_ID.bam $REFERENCE_ANNOTATION > $COUNTS_DIR/NGS$FILE_ID
htseq-count -f bam -s no -i gene_id SRR11785645_Aligned.sortedByCoord.out.bam rn6.gtf > SRR11785645_count


##########################
# Merge aligns per origin
##########################

# WT merged
samtools merge ./merged_bams/RNA_MIRANTA_WT.sorted.merged.bam ./STAR/*WT*.bam
samtools index ./merged_bams/RNA_MIRANTA_WT.sorted.merged.bam

bamCoverage --binSize 10 \
        --normalizeUsing RPKM \
        --maxFragmentLength 0 \
        -b ./merged_bams/RNA_MIRANTA_WT.sorted.merged.bam \
        -o ./bws/RNA_MIRANTA_WT.sorted.merged.bw \
        -p 24


# KO merged
samtools merge ./merged_bams/RNA_MIRANTA_KO.sorted.merged.bam ./STAR/*KO*.bam
samtools index ./merged_bams/RNA_MIRANTA_KO.sorted.merged.bam

bamCoverage --binSize 10 \
        --normalizeUsing RPKM \
        --maxFragmentLength 0 \
        -b ./merged_bams/RNA_MIRANTA_KO.sorted.merged.bam \
        -o ./bws/RNA_MIRANTA_KO.sorted.merged.bw \
        -p 24
