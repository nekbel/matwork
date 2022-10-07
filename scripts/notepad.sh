##########################
# PROTOTYPE SCRIPT
##########################

#!/bin/bash
# SCRIPTS USED:
#       orgUCSC
source ~/.nek_bash_commands.sh #better to source just the function 
echo "Imported source commands"

ref_build="${1:-}"
input="${2:-}"

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
mkdir -p


# index file first before bw creation
# do not mv .bai, bam files have to be in the same directory

##########################
# default value :-
##########################
refgen="${1:-hg38}"


##########################
# wget directory, accept/reject files
##########################
wget -r -nH -nd -np -R "*_random.fa.gz, chrUn_*" https://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/
# -r : recursive
# -nH : Disable generation of host-prefixed directories
# -nd : all files will get saved to the current directory
# -np : Do not ever ascend to the parent directory when retrieving recursively. 
# -R : don't download files with this files pattern
# -A : get only *.html files (for this case)


##########################
# LAVIGNE RNA PIPELINE LOOP
##########################
FILE="lane1RNA2MIRANTAWT2"

# alignment
hisat2 -p24 \
-x $REFERENCE_GENOME \
-1 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA2MIRANTAWT2_1_sequence.txt \
-2 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA2MIRANTAWT2_2_sequence.txt \
-S sam/$FILE.sam \
--rna-strandness FR \
--score-min L,0,-0.5 \
-k 3

samtools sort -@ 32 -o sam-bam/$FILE.bam sam/$FILE.sam
rm sam/$FILE.sam
samtools index sam-bam/$FILE.bam

mv sam-bam/$FILE.bam /data/www/Lavigne/tracks/miranta_bams/
mv sam-bam/$FILE.bam.bai /data/www/Lavigne/tracks/miranta_bams/

echo 'track type=bam pairEndsByName=. name=lane1RNA2MIRANTAWT2_bam description=lane1RNA2MIRANTAWT2 visibility=pack bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA2MIRANTAWT2.bam'

# WT & KO MERGE 
samtools merge -@ 32 /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.bam /data/www/Lavigne/tracks/miranta_bams/*RNA*WT*.bam
samtools sort -@ 32 -o /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.sorted.bam /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.bam
samtools index /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.sorted.bam

# WT & KO BW
bamCoverage \
--binSize 10 \
--normalizeUsing RPKM \
--maxFragmentLength 0 \
-b RNA_MIRANTA_WT.sorted.merged.bam \
-o RNA_MIRANTA_WT.sorted.merged.bw \
-p 24

echo 'track type=bigWig name=RNA_MIRANTA_WT.merged_bw description=RNA_MIRANTA_WT.merged_bw visibility=full maxHeightPixels=64:64:11 color=50,50,150 bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.sorted.bw'

# Process others
bamCoverage \
--binSize 10 \
--normalizeUsingRPKM \
--maxFragmentLength 0 \
-b /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1.bam \
-o /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1.bw \
-p 24
echo 'track type=bigWig name=lane1RNA1MIRANTAWT1_bw description=lane1RNA1MIRANTAWT1_fwd visibility=full bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1.bw'

bamCoverage \
--binSize 10 \
--normalizeUsingRPKM \
--maxFragmentLength 0 \
--filterRNAstrand=forward \
-b /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1.bam \
-o /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1_fwd.bw \
-p 24
echo 'track type=bigWig name=lane1RNA1MIRANTAWT1_fwd_bw description=lane1RNA1MIRANTAWT1_fwd visibility=full bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1_fwd.bw'

bamCoverage \
--binSize 10 \
--normalizeUsingRPKM \
--maxFragmentLength 0 \
--filterRNAstrand=reverse \
-b /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1.bam \
-o /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1_rev.bw \
-p 24
echo 'track type=bigWig name=lane1RNA1MIRANTAWT1_rev_bw description=lane1RNA1MIRANTAWT1_rev visibility=full bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1_rev.bw'


##########################
# disk used
##########################
du -hs * | sort -h


##########################
# find string contained in file 
##########################
grep -rnwl ./ -e 'SRR11785646'


##########################
# GIT
# Matwork Repository
##########################
git add <files>
git commit -m "commit message"
git push -u origin main
token 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

git log 

git checkout <commit-hash> (chech past commit)
git chekout master (return to last commit)

git branch (list branches)
git branch <new-branch-name>
git checkout <branch-name> (check branch)

git merge <branch-name> # merge branch to master


##########################
# kati 8a kanei kai auto mwre
##########################
cat <file.fa> | grep ">" | perl -ne '$_ =~ s/\>//; print $_' | sort | uniq > tidlist.txt


##########################
# GeneCounts STAR
##########################
STAR \
--runThreadN 6  \
--genomeDir /home/myrto/work/refs/Rattus_norvegicus/UCSC/rn6/Sequence/STARIndex \
--sjdbGTFfile /home/myrto/work/refs/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf \
--sjdbGTFtagExonParentTranscript transcript_id \
--readFilesIn sub* \
--quantMode GeneCounts


##########################
# GET FASTQ from SRR
##########################
curl -X GET \
"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR12080417&fields=fastq_ftp&result=read_run" |
cut -f 2 |
sed -n 2p |
tr \; "\n" |
wget -P ./fqs -i - &&
gunzip *.gz


##########################
# CHIP SEQ ALIGN BOWTIE2 (SINGLE-END)
##########################
trim_galore -o ./trimmed_fqs -q 20 ./fqs/SRR12080416; trim_galore -o ./trimmed_fqs -q 20 ./fqs/SRR12080417

bowtie2 --local -X 3000 -p 25 -x /home/myrto/work/refs/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -U ./trimmed_fqs/SRR12080416_trimmed.fq | samtools sort -O BAM - > ./bowtie/SRR12080416.sorted.bam
bowtie2 --local -X 3000 -p 25 -x /home/myrto/work/refs/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -U ./trimmed_fqs/SRR12080417_trimmed.fq | samtools sort -O BAM - > ./bowtie/SRR12080417.sorted.bam