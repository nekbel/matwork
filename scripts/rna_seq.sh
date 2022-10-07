#!/bin/bash
# SRR11785654

# check ! SRR11785688 (RNA)
# add conda activate portokali
mkdir -p fqs trimmed_fqs STAR bais bws stranded_bws
input="./srrList.txt"

while IFS="" read -r line || [ -n "$line" ]
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

        # uscs viewing stuff

        # indexing
        # samtools index "./STAR/${line}_Aligned.sortedByCoord.out.bam"
        # mv ${line}_Aligned.sortedByCoord.out.bam ./bais

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

        done < $input


echo $'\n'"Counting reads per feature..."

#/data/refs/Arabidopsis_thaliana/TAIR10/Arabidopsis_thaliana.TAIR10.47.gtf
#REFERENCE_ANNOTATION="/data/refs/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
#/data/refs/Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.gtf
REFERENCE_ANNOTATION="/data/refs/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
#/data/refs/Phaeodactylum_tricornutum/Phatr2_chromosomes_GeneModels_AllModels_20070523.gff
#/data/refs/Nicotiana_benthamiana/Niben101_annotation.allfeatures.gff
#REFERENCE_ANNOTATION="/data/refs/Botrytis_cinerea/ASM83294v1/genes.gff3"
#/data/refs/Cucumis_sativus/Ensembl/ASM407v2/Cucumis_sativus.ASM407v2.47.gtf

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



SUBSTRING=$(echo $INPUT| cut -d'_' -f 2)


# (syn75) 61 -> 5
# cut -d'-' -f 3
RNA1MIRANTAWT1_1
RNA4MIRANTAKO2_1
RNA1MIRANTAWT1_2
RNA4MIRANTAKO2_2
RNA2MIRANTAWT2_1
RNA5MIRANTAKO3_1
RNA2MIRANTAWT2_2
RNA5MIRANTAKO3_2
RNA3MIRANTAWT3_1
RNA6MIRANTAKO4_1
RNA3MIRANTAWT3_2
RNA6MIRANTAKO4_2

ATAC1MIRANTAWT1
ATAC2MIRANTAWT2
ATAC3MIRANTAWT3
ATAC4MIRANTAKO2
ATAC5MIRANTAKO3
ATAC6MIRANTAKO4


RNA1MIRANTAWT1_sorted.bam

AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA1MIRANTAWT1_1_sequence.txt
AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA4MIRANTAKO2_1_sequence.txt
AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA1MIRANTAWT1_2_sequence.txt
AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA4MIRANTAKO2_2_sequence.txt
AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA2MIRANTAWT2_1_sequence.txt
AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA5MIRANTAKO3_1_sequence.txt
AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA2MIRANTAWT2_2_sequence.txt
AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA5MIRANTAKO3_2_sequence.txt
AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA3MIRANTAWT3_1_sequence.txt
AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA6MIRANTAKO4_1_sequence.txt
AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA3MIRANTAWT3_2_sequence.txt
AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA6MIRANTAKO4_2_sequence.txt

#!/bin/bash

mkdir -p fqs trimmed_fqs STAR bais bws stranded_bws

input="./srr_lists/RNA_cortex.txt"

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
        trim_galore -q 20 -o ./trimmed_fqs --paired --nextera --fastqc ./fqs/${line}*

        # allign using STAR
        STAR    --genomeDir ~/work/genomes/STAR \
                --readFilesIn ./trimmed_fqs/${line}*.fq \
                --outFileNamePrefix ./STAR/${line}_ \
                --outSAMtype BAM SortedByCoordinate
                
        # uscs viewing stuff

        # indexing, what does bai stand for ?
        samtools index "./STAR/${line}_Aligned.sortedByCoord.out.bam" > ./bais/aligned_${line}.bai

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

#!/bin/bash

if [ -z $1 ] ; then
    echo 'give path to file'
    exit 1
else
    echo 'running for file ${1}...'
fi

mkdir -p fqs trimmed_fqs STAR bais bws stranded_bws
input=$1

while IFS="" read -r line || [ -n "$line" ] ; do

        date +'%m/%d/%Y %r'
        echo "! starting on ${line} !"

        # allign using STAR
        STAR    --genomeDir ~/work/genomes/STAR \
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
        # do not, bam files have to be in the same directory, maybe later for storage
        # mv ./STAR/${line}_Aligned.sortedByCoord.out.bam.bai ./bais/aligned_${line}.bai

        # kallisto rf-stranded ???
        kallisto quant \
        		--rf-stranded \
        		--index=/home/myrto/work/genomes/kallisto/rn6_cdna_kallisto_index \
        		--output-dir=kallisto/${line} \
        		--threads=4 \
        		./trimmed_fqs/*${line}*.fq

        done < $input 
# does grep -v '^#' work here ?


trim_galore -q 20 \
			-o ./trimmed_fqs \
			--paired \
			--fastqc \
			./fqs/*${line}*.txt




# STAR OPTION STRAND BW
# --outWigStrand
# default: Stranded

# string: strandedness of wiggle/bedGraph output

# Stranded
# separate strands, str1 and str2
# Unstranded
# collapsed strands

################################################################
# If you use STAR's --quantMode GeneCounts option, 
# it will output counts for 3 types of strandedness: 
# unstranded, forward, reverse, so you would have to 
# select the right orientation before sending it to DEseq. 
# If you count with featureCounts, you would have to specify 
# the option to count alignments in a certain strand orientation.
###############################################################


# --sjdbGTFfile

RNA1MIRANTAWT1

track type=bigWig name=RNA1MIRANTAWT1 description=RNA1MIRANTAWT1_both_strands visibility=full bigDataUrl=http://139.91.171.49/Lavigne/Nektarios/tracks/miranta_PV/aligned_RNA1MIRANTAWT1.bw
track type=bigWig name=RNA1MIRANTAWT1_rv description=RNA1MIRANTAWT1_rv_strand visibility=full bigDataUrl=http://139.91.171.49/Lavigne/Nektarios/tracks/miranta_PV/aligned_rv_RNA1MIRANTAWT1.bw
track type=bigWig name=RNA1MIRANTAWT1_fw description=RNA1MIRANTAWT1_fw_strand visibility=full bigDataUrl=http://139.91.171.49/Lavigne/Nektarios/tracks/miranta_PV/aligned_fw_RNA1MIRANTAWT1.bw

input='RNA_hippo.txt'

# bamtools merge -list ../srr_lists/RNA_cortex.txtbamfn.txt -out RNA_CARULLO_CORTEX.sorted.merged.bam

sed -e "s/$/_Aligned.sortedByCoord.out.bam/" ~/work/nek/RNA_seq/srr_lists/${input} > ~/work/nek/RNA_seq/srr_lists/${input}bamfn.txt

bamtools merge -list ~/work/nek/RNA_seq/srr_lists/${input}bamfn.txt -out RNA_{tissue}.sorted.merged.bam

samtools index RNA_{tissue}.sorted.merged.bam

bamCoverage --binSize 10 \
        --normalizeUsing RPKM \
        --maxFragmentLength 0 \
        -b RNA_{tissue}.sorted.merged.bam \
        -o RNA_{tissue}.sorted.merged.bw \
        -p 24

input='RNA_hippo.txt'; sed -e "s/$/_Aligned.sortedByCoord.out.bam/" ~/work/nek/RNA_seq/srr_lists/${input} > ~/work/nek/RNA_seq/srr_lists/${input}bamfn.txt
input='RNA_hippo.txt'; bamtools merge -list ~/work/nek/RNA_seq/srr_lists/${input}bamfn.txt -out RNA_hippo.sorted.merged.bam






#kallisto run
#give WHOLE path to index, not relative, kallisto is stupid 
line='SRR11785645'; kallisto quant --rf-stranded --index=/home/myrto/work/genomes/kallisto/rn6_cdna_kallisto_index --output-dir=kallisto/${line} --threads=4 ./trimmed_fqs/${line}*.fq




#!/bin/bash

input=$1

mkdir -p kallisto

while IFS="" read -r line || [ -n "$line" ] ; do


        date +'%m/%d/%Y %r'
        echo "! starting on ${line} !"

        # kallisto
        kallisto quant \
                --rf-stranded \
                --index=/home/myrto/work/genomes/kallisto/mmm10_cdna_kallisto_index \
                --output-dir=./kallisto/${line} \
                --threads=4 \
                ./trimmed_fqs/*${line}*.fq

        done < $input

# check Carullo Gene Metadata



# MIRANTA

#!/bin/bash
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
mkdir -p fqs trimmed_fqs STAR bws merged_bams

while IFS="" read -r line || [ -n "$line" ] ; do

    date +'%m/%d/%Y %r'
    echo "! starting on ${line} !"

    trim_galore \
    	-q 20 \
		-o ./trimmed_fqs \
		--paired \
		--fastqc \
		./fqs/*${line}*.txt

    # allign using STAR / mm10
    STAR \
		--genomeDir ${gen_path}/Sequence/STARIndex \
		--readFilesIn ./trimmed_fqs/*${line}*.fq \
		--outFileNamePrefix ./STAR/${line}_ \
		--outSAMtype BAM SortedByCoordinate

    # uscs viewing stuff

    # indexing
    samtools index "./STAR/${line}_Aligned.sortedByCoord.out.bam"

    # creating bigwig file (bai files have to be in the same directory)
    bamCoverage -b "./STAR/${line}_Aligned.sortedByCoord.out.bam" \
            -o "./bws/aligned_${line}.bw"

    bamCoverage \
    	-b "./STAR/${line}_Aligned.sortedByCoord.out.bam" \
		-o "./bws/aligned_fw_${line}.bw" \
		-p 24 \
		--binSize 10 \
		--normalizeUsing RPKM \
		--maxFragmentLength 0 \
		--filterRNAstrand=forward

    bamCoverage 
    	-b "./STAR/${line}_Aligned.sortedByCoord.out.bam" \
		-o "./bws/aligned_rv_${line}.bw" \
		-p 24 \
		--binSize 10 \
		--normalizeUsing RPKM \
		--maxFragmentLength 0 \
		--filterRNAstrand=reverse
    done < $input


# Merge aligns per origin

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

# kallisto 
    # kallisto quant \
    #         --rf-stranded \
    #         --index=${gen_path}/KallistoIndex/${ref_build}_cdna_kallisto_index \
    #         --output-dir=kallisto \
    #         --threads=4 \
    #         --plaintext \
    #         ./trimmed_fqs/*${line}*.fq

    # cd kallisto; mv abundance.tsv ${line}_abundance.tsv; cd ..



#####################################################################################################################

#!/bin/bash
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


grep
-B print n lines before match
-A print n lines after match
--no-group-separator dont print -- separator between matches
-n print line number of match
-v return everything that DOES NOT match
-c return the number of matches rather than the actual matching lines
-H report the filename

cut -f 1,4 -d \t -

sort 
-k which column to sort by
-u only return unique values
 
scp [OPTION] [user@]SRC_HOST:]file1 [user@]DEST_HOST:]file2

# rsync works similarly to scp.
# v is for verbosity - means that it prints on the screen what is being copied
# -e ssh is for encryption - means that we want to use the ssh protocol for encryption of the file transfer

basename -> retrieve filename without leading directories (prefix). Also remove specified suffix

variables -> there cannot be any spaces in their creation

.bashrc and .bash_profile located on ~
# You only log in once, and that's when ~/.bash_profile or ~/.profile is read and executed. 
# environment variables (LESS, PATH, MANPATH, LC_*, ...) e.g. 
# On the next runned shells, your environment variables will be passed along, so your non-login shells don't need to load them anymore
# Non-login shells only execute ~/.bashrc, not /.profile or ~/.bash_profile, for this exact reason, 
# define everything that only applies to bash (functions, aliases, bash-only variables like HISTSIZE 
# (this is not an environment variable, shell options with set and shopt
# You'll need to add in my ~/.profile too: source ~/.bashrc.

alias rm='rm -i'

.bashrc versus .bash_profile
# .bash_profile is executed for login shells, while .bashrc is executed for interactive non-login shells. 
# When you login (type username and password) to O2 the .bash_profile is executed. So if you want the alias 
# available only when you login, you will want to put it in your .bash_profile.

# download files and dirs
wget
-r/--recursive 
-np/--no-parent & -nd no parent directories 
-R "index.html*" reject this file

# download individual file (not dirs)
cURL
-O will save the file in the current working directory with the same file name as remote. 
-o lets you specify a different file name or location
-l list contents (directory -> directories and files)
				 (file -> print content in stdout)

https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/

curl -l <ftp.link> | grep href= 
returs list of files and directories contained

| grep href= | head -n -1 | sed 's/<a href="//g' | egrep -o '^.+"' | sed 's/.$//' | tail -n +2

fastq-dump  --split-files vs --split-3

CollectMultipleMetrics

traditional NGS analyses (e.g., QC, read mapping, variant calling, assembly etc.) .

• Transcriptome gives us the complete set of transcripts
It is in the form of a GTF (gene transfer format) file

Transcriptome reference file

Mapping results are only as good as the quality of the
reference transcriptome

• If one does not exist, it can be created using coordinates
from a GTF file and the genome sequence file 

salmon_in_serial
Request 6 cores 
Request 12 hours of runtime
Request 8G of memory


for fq in ~/rnaseq/raw_data/*.fq

do

# create a prefix for the output file
samplename=`basename $fq .fq`

# run salmon
salmon quant -i /n/holylfs05/LABS/hsph_bioinfo/Everyone/Workshops/Intro_to_rnaseq/indicies/salmon_index \
 -l A \
 -r $fq \
 -o ${samplename}.salmon \
 -p 6 \
 --numBootstraps 30 \ # 
 --seqBias \
 --useVBOpt \
 --validateMappings

done

# Bootstrapping is any test or metric that uses random sampling with replacement and falls under the broader 
# class of resampling methods. It uses sampling with replacement to estimate the sampling distribution for 
# the desired estimator. This approach is used to assess the reliability of sequence-based phylogeny.

# Bootstrapping is a statistical procedure that resamples a single dataset to create many simulated samples. 
# This process allows for the calculation of standard errors, confidence intervals, and hypothesis testing

# NOTE: --numBootstraps is necessary if performing isoform-level differential expression analysis with Sleuth, 
# but not for gene-level differential expression analysis. Due to the statistical procedure required to assign 
# reads to gene isoforms, in addition to the random processes underlying RNA-Seq, there will be technical 
# variability in the abundance estimates output from the pseudo-alignment tool [2, 3] for the isoform level 
# abundance estimates (not necessary for gene-level estimates). Therefore, we would need technical replicates 
# to distinguish technical variability from the biological variability for gene isoforms.

# The bootstraps estimate technical variation per gene by calculating the abundance estimates for all genes 
# using a different sub-sample of reads during each round of bootstrapping. The variation in the abundance 
# estimates output from each round of bootstrapping is used for the estimation of the technical variance for 
# each gene.

GRCh38, contains information about alternative alleles for various locations on the genome, use the HISAT2 aligner
as it is able to utilize this information during the alignment. 
There is a version of GRCh38 that does not have alleles represented, which is the appropriate version to use with STAR. 

STAR
--runThreadN: number of threads / cores
--outSAMtype: output filetype (SAM default)
--outSAMunmapped: what to do with unmapped reads (Within)
--outSAMattributes Standard 

NOTE: Note that “STAR’s default parameters are optimized for mammalian genomes. Other species may require significant modifications of some alignment parameters; in particular, the maximum and minimum intron sizes have to be reduced for organisms with smaller introns” [1].

qualimap rnaseq \
-outdir results/qualimap/Mov10_oe_1 \
-a proportional \
-bam results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam \
-p strand-specific-reverse \
-gtf /n/holylfs05/LABS/hsph_bioinfo/Everyone/Workshops/Intro_to_rnaseq/reference_data/Homo_sapiens.GRCh38.92.gtf \
--java-mem-size=8G

RNASEQ Metrics
number of raw reads
percentage of reads aligned to genome
percentage of reads associated with genes

multiqc -n multiqc_report_rnaseq \
full_dataset_results/fastqc/*zip \
full_dataset_results/STAR/*Log.final.out \
full_dataset_results/qualimap/* \
full_dataset_results/salmon/*salmon

always use full path for parameters of toools

# * DO NOT RUN *
sh  run_rnaseq.sh  input.fq  input.gtf  12
$0 => run_rnaseq.sh

$1 => input.fq

$2 => input.gtf

$3 => 12



#!/bin/bash
# This script takes a fastq file of RNA-seq data, runs FastQC, STAR, Qualimap and Salmon.
# USAGE: sh rnaseq_analysis_on_input_file.sh <name of fastq file>

set -x
fq=$1

# grab base of filename for naming outputs

samplename=`basename $fq .subset.fq`
echo "Sample name is $samplename"         

# specify the number of cores to use

cores=6

# directory with the genome and transcriptome index files + name of the gene annotation file

genome=/n/holylfs05/LABS/hsph_bioinfo/Everyone/Workshops/Intro_to_rnaseq/indicies/ensembl38_STAR_index
transcriptome=/n/holylfs05/LABS/hsph_bioinfo/Everyone/Workshops/Intro_to_rnaseq/indicies/salmon_index
gtf=/n/holylfs05/LABS/hsph_bioinfo/Everyone/Workshops/Intro_to_rnaseq/reference_data/Homo_sapiens.GRCh38.92.gtf

# make all of the output directories
# The -p option means mkdir will create the whole path if it 
# does not exist and refrain from complaining if it does exist

mkdir -p results/fastqc/
mkdir -p results/STAR/
mkdir -p results/qualimap/
mkdir -p results/salmon/

# set up output filenames and locations

fastqc_out=results/fastqc/
align_out=results/STAR/${samplename}_
align_out_bam=results/STAR/${samplename}_Aligned.sortedByCoord.out.bam
qualimap_out=results/qualimap/${samplename}.qualimap
salmon_out=results/salmon/${samplename}.salmon
salmon_mappings=results/salmon/${samplename}_salmon.out

# set up the software environment (use version numbers)

module load fastqc/0.11.8-fasrc01
module load STAR/2.7.0e-fasrc01
unset DISPLAY

echo "Processing file $fq"

echo "Starting QC for $samplename"

# Run FastQC and move output to the appropriate folder
fastqc -o $fastqc_out $fq


# Run STAR
STAR --runThreadN $cores \
	--genomeDir $genome \
	--readFilesIn $fq \
	--outFileNamePrefix $align_out \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes Standard

# Run Qualimap
qualimap rnaseq \
-outdir $qualimap_out \
-a proportional \
-bam $align_out_bam \
-p strand-specific-reverse \
-gtf $gtf \
--java-mem-size=8G

# Run salmon

module load salmon/1.5.2-fasrc01

echo "Starting Salmon run for $samplename"

salmon quant -i $transcriptome \
-p $cores \
-l A \
-r $fq \
-o $salmon_out \
--seqBias \
--useVBOpt








# this `for` loop, will take the fastq files as input and run the script for all of them one after the other. 
for fq in ~/rnaseq/raw_data/*.fq
do
  echo "running analysis on $fq"
  sh ~/rnaseq/scripts/rnaseq_analysis_on_input_file.sh $fq
done



wrappedProgram "$@"
# ^^^ this is correct and will hand over all arguments in the way
#     we received them, i. e. as several arguments, each of them
#     containing all the spaces and other uglinesses they have.
wrappedProgram "$*"
# ^^^ this will hand over exactly one argument, containing all
#     original arguments, separated by single spaces.
wrappedProgram $*
# ^^^ this will join all arguments by single spaces as well and
#     will then split the string as the shell does on the command
#     line, thus it will split an argument containing spaces into
#     several arguments.



for f in ./*.txt
do
  echo $f | sed 's/txt/fq/g'
done


#!/bin/bash

input=$1

while IFS="" read -r line || [ -n "$line" ]; do

    date +'%m/%d/%Y %r'
    echo "! starting on ${line} !"
	salmon quant -i /home/myrto/work/refs/Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/SalmonIndex \
		-p 8 \
		-l A \
		-1 ../fqs/*${line}_1* \
		-2 ../fqs/*${line}_2* \
		-o ${line} \
		--seqBias \
		--useVBOpt

    done < $input


miranta library type -> ISR Inward Stranded Reverse
Carullo library type -> ISR