#!/bin/bash
# SRR11785646.sra
# add conda activate portokali

# list required packages (prob on conda enviroment)
# create folders (-p -> if they don't already exist)
mkdir -p fqs trimmed_fqs bowtie macs macs/sorted macs/merged bws stranded_bws
input = "./srrList.txt"

while IFS="" read -r line || [ -n "$line" ]

		do

        date +'%m/%d/%Y %r'
        echo "starting on \"${line}\"... "

        # download SRR files from NCBI (prob need an alternative) -> download from EBI!
        # prefetch $line

        # convert SRA to fastq files
        # fastq-dump --split-files --outdir ./fqs "~/ncbi/public/sra/${line}.sra"

        # need to time these
        # get time on /usr/bin/time on the server, check on local compu, or maya!

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

        # alignment, requires bowtie2-build for your genome beforehand
        # bowtie2-build rn6.fa rn6
        echo 'export GENOME_DB=~/work/genomes' >> ~/.profile
		echo 'export BOWTIE2_INDEXES=~/work/genomes/bowtie2_indexes' >> ~/.profile
		source ~/.profile
        # maybe can be swapped out with STAR
        bowtie2 --local \
        		--very-sensitive-local \
        		-X 3000 \
        		-x rn6 \
        		-1 ./trimmed_fqs/${line}_1_val_1.fq \
        		-2 ./trimmed_fqs/${line}_2_val_2.fq |
       	grep -v '@SQ' |
        samtools sort -O BAM - > ./bowtie/aligned_${line}_sorted_clean.bam 

        # -X: fragment size
        # remove @SQ wste 165_266_405 -> 54_933_470
        # samtools sort can take SAM file -> sorted BAM in one go (instead of view-> sort)
        # the (-) on samtools sort is for taking in piped input 

        # peak calling

        # mm10 (mouse) effective genome size = 
        # non-N bases = 
        # oliko (2,730,871,774) - N (78,088,274) = 
        # 2,652,783,500

        # rn6 (rat) egs = 2,729,860,805

        # --format BAMPE = bam paired end

        # macs2 can take bam and sam input, (un)sorted. Bam is faster
		macs2 callpeak \
			--treatment ./bowtie/aligned_${line}_sorted_clean.bam \
       		--qvalue 0.001 \
      		--gsize 2729860805 \
       		--format BAMPE \
       		--outdir ./macs \
       		--name ${line}
       		# bw option, 

       	sort -k1,1 -k2,2n ./macs/${line}_peaks.bed > ./macs/sorted/${line}_sortedPeaks.bed
		bedtools merge -i ./macs/sorted/${line}_sortedPeaks.bed -d 1000 > ./macs/${line}_sortedMergedPeaks.bed 
 		# BED6+4 ---> BED to be uploaded to ucsc (not sure we need to do that)
		cut -f 1-6 ./macs/${line}_peaks.narrowPeak > ./macs/${line}_peaks.bed
		
		# sort and merge peaks
		

		# extras (bws viewing)

		bamCoverage -b "./bowtie/aligned_${line}_sorted_clean.bam" -of bigwig -o "./bws/aligned_${line}.bw"

        bamCoverage -b "./bowtie/aligned_${line}_sorted_clean.bam" \
                -o "./stranded_bws/aligned_fw_${line}.bw" \
                -p 24 \
                --binSize 10 \
                --normalizeUsing RPKM \
                --maxFragmentLength 0 \
                --filterRNAstrand=forward

        bamCoverage -b "./bowtie/aligned_${line}_sorted_clean.bam" \
                -o "./stranded_bws/aligned_rv_${line}.bw" \
                -p 24 \
                --binSize 10 \
                --normalizeUsing RPKM \
                --maxFragmentLength 0 \
                --filterRNAstrand=reverse

        done < $input




#!/bin/bash

mkdir -p fqs trimmed_fqs bowtie macs macs/sorted macs/merged macs/bed6 bws stranded_bws
input="./test_atac.txt"

while IFS="" read -r line || [ -n "$line" ]

        do

        date +'%m/%d/%Y %r'
        echo "starting on \"${line}\"... "

        #curl -X GET \
        #"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${line}&fields=fastq_ftp&result=read_run" |
        #cut -f 2 |
        #sed -n 2p |
        #tr \; "\n" |
        #wget -P ./fqs -i - &&
        #gunzip ./fqs/*.gz

        # trim off primers
        trim_galore -o ./trimmed_fqs \
                                --nextera \
                                --fastqc \
                                --paired \
                                -q 20 \
                                ./fqs/${line}*

        # alignment, requires bowtie2-build
        # samtools sort requires the @SQ headers, so we could removre afterwards
        bowtie2 --local \
                --very-sensitive-local \
                -X 3000 \
                -x rn6 \
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
        sort -k1,1 -k2,2n ./macs/${line}_peaks.bed > ./macs/sorted/${line}_sortedPeaks.bed
        bedtools merge -i ./macs/sorted/${line}_sortedPeaks.bed -d 1000 > ./macs/${line}_sortedMergedPeaks.bed

        # BED6+4 -> BED
        cut -f 1-6 ./macs/${line}_sortedMergedPeaks.bed  > ./macs/bed6/${line}_peaks.6.bed

        # extras (bws viewing)
        bamCoverage -b "./bowtie/${line}.sorted.bam" \
                -o "./bws/${line}.bw"

        bamCoverage -b "./bowtie/${line}.sorted.bam" \
                -o "./stranded_bws/fw_${line}.bw" \
                -p 24 \
                --binSize 10 \
                --normalizeUsing RPKM \
                --maxFragmentLength 0 \
                --filterRNAstrand=forward

        bamCoverage -b "./bowtie/${line}.sorted.bam" \
                -o "./stranded_bws/rv_${line}.bw" \
                -p 24 \
                --binSize 10 \
                --normalizeUsing RPKM \
                --maxFragmentLength 0 \
                --filterRNAstrand=reverse

        done < $input



*_random

samtools view -o RNA1MIRANTAWT1_chr1-Y.bam RNA1MIRANTAWT1_Aligned.sortedByCoord.out.bam `seq 1 22 | sed 's/^/chr/'; echo chrX; echo chrY`

# mm10
hisat
samtools view -H lane1RNA1MIRANTAWT1.bam | perl -lne '/SN:(\S+)/ and print $1'

1-19 chrX chrY chrM

samtools view RNA1MIRANTAWT1_Aligned.sortedByCoord.out.bam chrM >> RNA1MIRANTAWT1_chr1-Y.bam



samtools view -o .chrs.bam  `seq 1 19 | sed 's/^/chr/'; echo chrX; echo chrY; echo chrM`

bamCompare -b1  -b2  -o log2ratio.bw 
--effectiveGenomeSize 2652783500
--scaleFactorsMethod None
--normalizeUsing RPKM

bamCompare --effectiveGenomeSize 2652783500 --scaleFactorsMethod None --normalizeUsing RPKM -b1 ../hisat2/lane1RNA1MIRANTAWT1.chrs.bam -b2 ../STAR/RNA1MIRANTAWT1.chrs.bam -o RNA1MIRANTAWT1.log2ratio.bw


bwa mem -R "@RG\tID:lane1ATAC6MIRANTAKO3\tSM:lane1ATAC6MIRANTAKO3\tPL:illumina\tCN:EMBL" \
		-t 8 \
		/home/myrto/work/genomes/test/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa \
		./*ATAC6MIRANTAKO3* > ATAC6MIRANTAKO3.sam


bwa mem -R "@RG\tID:lane1ATAC6MIRANTAKO3\tSM:lane1ATAC6MIRANTAKO3\tPL:illumina\tCN:EMBL" -t 8 /home/myrto/work/refs/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa ./*ATAC6MIRANTAKO3* > ATAC6MIRANTAKO3.sam




#!/bin/bash

mkdir -p fqs trimmed_fqs bowtie macs macs/sorted macs/merged macs/bed6 bws stranded_bws
input=$1

while IFS="" read -r line || [ -n "$line" ]

        do

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


#        bamCoverage -b "./bowtie/${line}.sorted.bam" \
#                -o "./bws/${line}.bw"
#
#        bamCoverage -b "./bowtie/${line}.sorted.bam" \
#                	-o "./stranded_bws/fw_${line}.bw" \
#                	-p 24 \
#                	--binSize 10 \
#                	--normalizeUsing RPKM \
#                	--maxFragmentLength 0 \
#                	--filterRNAstrand=forward
#
#        bamCoverage -b "./bowtie/${line}.sorted.bam" \
#        	        -o "./stranded_bws/rv_${line}.bw" \
#        	        -p 24 \
#        	        --binSize 10 \
#        	        --normalizeUsing RPKM \
#        	        --maxFragmentLength 0 \
#        	        --filterRNAstrand=reverse




##################################################################
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