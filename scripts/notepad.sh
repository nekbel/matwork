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


##########################
# track comp hub 
##########################
track type=bigWig name=RNA_MIRANTA_KO.sorted.merged description=RNA_MIRANTA_KO.sorted.merged visibility=full bigDataUrl=http://139.91.171.49/Lavigne/Nektarios/tracks/miranta_PV/RNA_MIRANTA_KO.sorted.merged.bw color=0,255,0
track type=bigWig name=RNA_MIRANTA_WT.sorted.merged description=RNA_MIRANTA_WT.sorted.merged visibility=full bigDataUrl=http://139.91.171.49/Lavigne/Nektarios/tracks/miranta_PV/RNA_MIRANTA_WT.sorted.merged.bw color=255,0,0

hub MirantaHub
shortLabel Miranta Hub
longLabel RNA and ATAC runs on the Miranta files
useOneFile on
email nek.belmezos@gmail.com

genome mm10

track ComparisonKOWT
type bigWig
container multiWig
visibility full
aggregate transparentOverlay
showSubtrackColorOnUi on
maxHeightPixels 500:100:8

        track KO
        type bigWig
        bigDataUrl http://139.91.171.49/Lavigne/Nektarios/tracks/miranta_PV/RNA_MIRANTA_KO.sorted.merged.bw
        parent ComparisonKOWT
        color 70,250,150

        track WT
        type bigWig
        bigDataUrl http://139.91.171.49/Lavigne/Nektarios/tracks/miranta_PV/RNA_MIRANTA_WT.sorted.merged.bw
        parent ComparisonKOWT
        color 207,21,62


##########################
# Hub track example
##########################
track multiWig1
type bigWig
container multiWig
shortLabel Ex. multiWig container
longLabel This multiWig overlay track graphs points from three bigWig files.
visibility full
aggregate transparentOverlay
showSubtrackColorOnUi on
maxHeightPixels 500:100:8
viewLimits 1:20
priority 1
html examplePage

				track wig1
				bigDataUrl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep1.bigWig
				shortLabel Overlay bigWig1
				longLabel This is an example bigWig1 displaying Raw Signal from the ENCODE RNA-seq CSHL track, graphing just points as default.
				parent multiWig1
				graphTypeDefault points
				type bigWig
				color 255,0,0

				track wig2
				bigDataUrl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep1.bigWig
				shortLabel Overlay bigWig2
				longLabel This is an example bigWig2 displaying Raw Signal from the ENCODE RNA-seq CSHL track, graphing just points as default.
				graphTypeDefault points
				parent multiWig1
				type bigWig
				color 0,255,0

				track wig3
				bigDataUrl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/wgEncodeCshlLongRnaSeqAg04450CellLongnonpolyaPlusRawSigRep1.bigWig
				shortLabel Overlay bigWig3
				longLabel This is an example bigWig3 displaying Raw Signal from the ENCODE RNA-seq CSHL track, graphing just points as default. 
				graphTypeDefault points
				parent multiWig1
				type bigWig
				color 95,158,160