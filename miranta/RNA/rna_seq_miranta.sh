#!/bin/bash

refgen="${1:-}"
input="${2:-}"

mkdir -p fqs trimmed_fqs STAR bws

while IFS="" read -r line || [ -n "$line" ] ; do

        date +'%m/%d/%Y %r'
        echo "! starting on ${line} !"

	trim_galore -q 20 \
		-o ./trimmed_fqs \
		--paired \
		--fastqc \
		./fqs/*${line}*.txt

        # allign using STAR / mm10
        STAR    --genomeDir ~/work/genomes/STAR/${refgen} \
                --readFilesIn ./trimmed_fqs/*${line}*.fq \
                --outFileNamePrefix ./STAR/${line}_ \
                --outSAMtype BAM SortedByCoordinate

        # uscs viewing stuff

        # indexing
        samtools index "./STAR/${line}_Aligned.sortedByCoord.out.bam"

        # creating bigwig file (bai files have to be in the same directory)
        bamCoverage -b "./STAR/${line}_Aligned.sortedByCoord.out.bam" \
                -o "./bws/aligned_${line}.bw"

        bamCoverage -b "./STAR/${line}_Aligned.sortedByCoord.out.bam" \
                -o "./bws/aligned_fw_${line}.bw" \
                -p 24 \
                --binSize 10 \
                --normalizeUsing RPKM \
                --maxFragmentLength 0 \
                --filterRNAstrand=forward

        bamCoverage -b "./STAR/${line}_Aligned.sortedByCoord.out.bam" \
                -o "./bws/aligned_rv_${line}.bw" \
                -p 24 \
                --binSize 10 \
                --normalizeUsing RPKM \
                --maxFragmentLength 0 \
                --filterRNAstrand=reverse
        # do not, bam files have to be in the same directory, maybe later for storage
        # mv ./STAR/${line}_Aligned.sortedByCoord.out.bam.bai ./bais/aligned_${line}.bai

	#kallisto
        kallisto quant \
		--rf-stranded \
		--index=/home/myrto/work/genomes/kallisto/mmm10_cdna_kallisto_index \
		--output-dir=kallisto \
		--threads=4 \
		--plaintext \
		./trimmed_fqs/*${line}*.fq

        cd kallisto; mv abundance.tsv ${line}_abundance.tsv; cd ..
        done < $input

# WT merged
samtools merge RNA_MIRANTA_WT.sorted.merged.bam ./STAR/*WT*.bam
samtools index ./STAR/RNA_MIRANTA_WT.sorted.merged.bam

bamCoverage --binSize 10 \
	--normalizeUsing RPKM \
	--maxFragmentLength 0 \
	-b RNA_MIRANTA_WT.sorted.merged.bam \
	-o RNA_MIRANTA_WT.sorted.merged.bw \
	-p 24


# KO merged
samtools merge RNA_MIRANTA_KO.sorted.merged.bam ./*KO*.bam
samtools index RNA_MIRANTA_KO.sorted.merged.bam

bamCoverage --binSize 10 \
	--normalizeUsing RPKM \
	--maxFragmentLength 0 \
	-b RNA_MIRANTA_KO.sorted.merged.bam \
	-o RNA_MIRANTA_KO.sorted.merged.bw \
	-p 24
